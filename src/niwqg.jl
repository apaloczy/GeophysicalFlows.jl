module NIWQG

export
  Problem,
  set_q!,
  set_phi!,
  set_planewave!,
  updatevars!,

  waveaction,
  qgke,
  wavepe,
  coupledenergy
  
using 
  Reexport,
  FFTW

@reexport using FourierFlows

using LinearAlgebra: mul!, ldiv!
using FourierFlows: getfieldspecs, varsexpression, parsevalsum, parsevalsum2


nothingfunction(args...) = nothing
onefunction(args...) = 1

# --
# Problem
# --

"""
    Problem(; parameters...)

Construct a VerticallyFourierBoussinesq initial value problem.
"""
function Problem(;
  # Numerical parameters
          nx = 64,
          Lx = 2π,
          ny = nx,
          Ly = Lx,
          dt = 1e-9,
  # Drag and/or hyper-/hypo-viscosity
          nu = 0, # wave (hyper-)viscosity
         nnu = 1, # wave (hyper-)viscous order (1=Laplacian)
         kap = 0, # PV (hyper-)visosity
        nkap = 1, # PV (hyper-)viscous order
         muq = 0, # drag/arbitrary-order dissipation for q
        nmuq = 0, # order of 2nd dissipation term for q
         muw = 0, # drag/arbitrary-order dissipation for w
        nmuw = 0, # order of 2nd dissipation term for w
  # Physical parameters
         eta = 0, # dispersivity: eta = N^2 / f m^2
           f = 1, # inertial frequency
  # Optional uniform and steady background flow
          Ub = 0,
          Vb = 0,
  # Timestepper and eqn options
     stepper = "RK4",
       calcF = nothingfunction,
  stochastic = false,
    nthreads = Sys.CPU_THREADS,
           T = Float64,
     dealias = false
  )

  grid = TwoDGrid(nx, Lx, ny, Ly; T=T)
  params = Params{T}(nu, nnu, kap, nkap, muw, nmuw, muq, nmuq, eta, f, 1/f, Ub, Vb, calcF)
  vars = (calcF == nothingfunction ? 
            Vars(grid) : 
            (stochastic ? StochasticForcedVars(grid) : ForcedVars(grid)))
  eqn = Equation(params, grid)

  FourierFlows.Problem(eqn, stepper, dt, grid, vars, params)
end


# --
# Params
# --

"""
    Params(nu, nnu, kap, nkap, muq, nmuq, muw, nmuw, eta, f, invf, Ub, Vb, calcF!)

Construct parameters for an NIWQG problem.
"""
struct Params{T} <: AbstractParams
  nu::T             # Wave viscosity
  nnu::Int          # Wave hyperviscous order
  kap::T            # PV 'diffusivity'
  nkap::Int         # PV hyperdiffusive order
  muw::T            # Wave drag / hypoviscosityty
  nmuw::Int         # Wave drag / hypoviscous orderer
  muq::T            # PV drag / hypoviscosity
  nmuq::Int         # PV drag / hypoviscous order
  eta::T            # Dispersivity
  f::T              # Planetary vorticity
  invf::T           # 1/f
  Ub::T             # Steady barotropic background ('mean') x-velocity
  Vb::T             # Steady barotropic background ('mean') y-velocity
  calcF!::Function  # Forcing
end

# ---------
# Equations
# ---------

function Equation(p, g)
  Lq = @. - p.kap * g.Krsq ^ p.nkap  -  p.muq * g.Krsq ^ p.nmuq
  Lw = @. - p.nu  * g.Ksq  ^ p.nnu   -  p.muw * g.Ksq  ^ p.nmuw  -  0.5*im*p.eta*g.Ksq
  Lq[1, 1] = 0
  Lw[1 ,1] = 0
  L = [Lq, Lw]
  FourierFlows.Equation(L, calcN!, g)
end


# --
# Vars
# --

physicalvarsr = [:q, :U, :V, :zeta, :psi, :Uq, :Vq, :phijac, :phisq]
physicalvarsc = [:phi, :Uphi, :Vphi, :zetaphi, :phix, :phiy] 
transformvarsr = [ Symbol(var, :h) for var in physicalvarsr ]
transformvarsc = [ Symbol(var, :h) for var in physicalvarsc ]
forcedvars = [:F]
stochforcedvars = [:prevsol]

varspecs = cat(getfieldspecs(physicalvarsr, :Tr),
               getfieldspecs(physicalvarsc, :Tc),
               getfieldspecs(transformvarsr, :Tc),
               getfieldspecs(transformvarsc, :Tc),
               dims=1)

forcedvarspecs = cat(varspecs, getfieldspecs(forcedvars, :Tsol), dims=1)
stochforcedvarspecs = cat(forcedvarspecs, getfieldspecs(stochforcedvars, :Tsol), dims=1)

eval(varsexpression(:Vars, varspecs, typeparams=[:Tr, :Tc]))
eval(varsexpression(:ForcedVars, forcedvarspecs, typeparams=[:Tr, :Tc, :Tsol]))
eval(varsexpression(:StochasticForcedVars, stochforcedvarspecs, typeparams=[:Tr, :Tc, :Tsol]))

function Vars(g::AbstractGrid{T}) where T
  @zeros T (g.nx, g.ny) q U V zeta psi Uq Vq phijac phisq
  @zeros Complex{T} (g.nkr, g.nl) qh Uh Vh zetah psih Uqh Vqh phijach phisqh
  @zeros Complex{T} (g.nx, g.ny) phi Uphi Vphi zetaphi phix phiy
  @zeros Complex{T} (g.nk, g.nl) phih Uphih Vphih zetaphih phixh phiyh
  Vars(
    q, U, V, zeta, psi, Uq, Vq, phijac, phisq,
    phi, Uphi, Vphi, zetaphi, phix, phiy,
    qh, Uh, Vh, zetah, psih, Uqh, Vqh, phijach, phisqh,
    phih, Uphih, Vphih, zetaphih, phixh, phiyh
  )
end

function ForcedVars(g::AbstractGrid{T}) where T
  v = Vars(g)
  F = [zeros(Complex{T}, (g.nkr, g.nl)), zeros(Complex{T}, (g.nk, g.nl))]
  ForcedVars(getfield.(Ref(v), fieldnames(typeof(v)))..., F)
end

function StochasticForcedVars(g::AbstractGrid{T}) where T
  v = ForcedVars(g)
  prevsol = [zeros(Complex{T}, (g.nkr, g.nl)), zeros(Complex{T}, (g.nk, g.nl))]
  StochasticForcedVars(getfield.(Ref(v), fieldnames(typeof(v)))..., prevsol)
end


# --
# Solver routines
# --

function calczetah!(zetah, qh, phih, phi, v, p, g)
  # Calc qw
  @. v.phixh = im*g.k*phih
  @. v.phiyh = im*g.l*phih

  ldiv!(v.phix, g.fftplan, v.phixh)
  ldiv!(v.phiy, g.fftplan, v.phiyh)

  @. v.phijac = real(im*(conj(v.phix)*v.phiy - conj(v.phiy)*v.phix))
  @. v.phisq = abs2(phi)

  mul!(v.phijach, g.rfftplan, v.phijac)
  mul!(v.phisqh, g.rfftplan, v.phisq)

  #   zeta = q  -               *** q^w ***
  @. zetah = qh - p.invf*(-0.25*g.Krsq*v.phisqh + 0.5*v.phijach)
  zetah[1, 1] = 0

  nothing
end

function calcUhVh!(Uh, Vh, psih, Ub, Vb, g)
  @. Uh = - im * g.l  * psih # U = - ∂y ψ
  @. Vh =   im * g.kr * psih # V = + ∂x ψ
  Uh[1, 1] += Ub*g.nx*g.ny
  Vh[1, 1] += Vb*g.nx*g.ny
  nothing
end

function calcadvectionterms!(Uq, Vq, Uphi, Vphi, U, V, q, phi)
  @. Uq = U*q
  @. Vq = V*q
  @. Uphi = U*phi
  @. Vphi = V*phi
  nothing
end

function calcN!(N, sol, t, cl, v, p, g)
  @views @. v.qh = sol[1] # inverse transforms may destroy sol[1]
  @views @. v.phih = sol[2] # this copy may not be necessary, but clarifies code

  ldiv!(v.phi, g.fftplan, sol[2])
  @. v.psih = -g.invKrsq*v.zetah
  calcUhVh!(v.Uh, v.Vh, v.psih, p.Ub, p.Vb, g)

  # inverse transforms may destroy v.Uh, v.Vh, v.zetah, v.qh
  ldiv!(v.U, g.rfftplan, v.Uh)
  ldiv!(v.V, g.rfftplan, v.Vh)
  ldiv!(v.zeta, g.rfftplan, v.zetah)
  ldiv!(v.q, g.rfftplan, v.qh)

  calcadvectionterms!(v.Uq, v.Vq, v.Uphi, v.Vphi, v.U, v.V, v.q, v.phi)
  @. v.zetaphi = v.zeta*v.phi

  mul!(v.Uphih, g.fftplan, v.Uphi)
  mul!(v.Vphih, g.fftplan, v.Vphi)
  mul!(v.Uqh, g.rfftplan, v.Uq)
  mul!(v.Vqh, g.rfftplan, v.Vq)
  mul!(v.zetaphih, g.fftplan, v.zetaphi)

  @views @. N[1] = - im*g.kr*v.Uqh - im*g.l*v.Vqh
  @views @. N[2] = - im*g.k*v.Uphih - im*g.l*v.Vphih - 0.5*im*v.zetaphih

  addforcing!(N, sol, t, cl, v, p, g)

  nothing
end

addforcing!(N, sol, t, cl, v::Vars, p, g) = nothing

function addforcing!(N, sol, t, cl, v::ForcedVars, p, g)
  p.calcF!(v.Fh, sol, t, cl, v, p, g)
  @views @. N[1] += v.Fh[1]
  @views @. N[2] += v.Fh[2]
  nothing
end

function addforcing!(N, sol, t, cl, v::StochasticForcedVars, p, g)
  if t == cl.t # not a substep
    @views @. v.prevsol[1] = sol[1]
    @views @. v.prevsol[2] = sol[2]
    p.calcF!(v.Fh, sol, t, cl, v, p, g)
  end
  @views @. N[1] += v.Fh[1]
  @views @. N[2] += v.Fh[2]
  nothing
end

# --
# Helper functions
# --

"""
    updatevars!(prob)

Update variables in `prob.vars` using the solution in `prob.sol`.
"""
function updatevars!(v, sol, p, g)
  @views @. v.qh .= sol[1]
  @views @. v.phih .= sol[2]

  ldiv!(v.phi, g.fftplan, v.phih)
  @views calczetah!(v.zetah, sol[1], sol[2], v.phi, v, p, g)

  @. v.psih = -g.invKrsq*v.zetah
  calcUhVh!(v.Uh, v.Vh, v.psih, p.Ub, p.Vb, g)

  # Copy variables that are destroyed by inverse transforms
  ldiv!(v.psi, g.rfftplan, deepcopy(v.psih))
  ldiv!(v.q, g.rfftplan, deepcopy(v.qh))
  ldiv!(v.U, g.rfftplan, deepcopy(v.Uh))
  ldiv!(v.V, g.rfftplan, deepcopy(v.Vh))

  ldiv!(v.phi, g.fftplan, v.phih)

  nothing
end
updatevars!(prob) = updatevars!(prob.vars, prob.sol, prob.params, prob.grid)

"""
    set_q!(prob, q)

Set potential vorticity and update variables.
"""
function set_q!(prob, q)
  @. prob.vars.q = q
  @views mul!(prob.sol[1], prob.grid.rfftplan, prob.vars.q)
  updatevars!(prob)
  nothing
end

"""
    set_phi!(prob, phi)

Set the wave field amplitude, phi.
"""
function set_phi!(prob, phi)
  @. prob.vars.phi = phi
  @views mul!(prob.sol[2], prob.grid.fftplan, prob.vars.phi)
  updatevars!(prob)
  nothing
end

"""
    set_planewave!(prob, uw, nkw, θ=0; kwargs...)

Set plane wave solution in `prob` with initial speed `uw` and non-dimensional wave
number `nkw`. Keyword argument `envelope=env(x, y)` multiplies the plane wave by `env(x, y)`
"""
function set_planewave!(prob, uw, nkw, θ=0; envelope=onefunction)
  k = 2π/g.Lx*round(Int, nkw*cos(θ))
  l = 2π/g.Lx*round(Int, nkw*sin(θ))
  x, y = gridpoints(prob.grid)

  Φ = @. k*x + l*y
  phi = @. uw * exp(im*Φ) * envelope(x, y)
  set_phi!(prob, phi)
  nothing
end

"""
    waveaction(prob)

Returns the domain-averaged near-inertial action.
"""
waveaction(prob) = waveaction(prob.sol, prob.params, prob.grid)

waveaction(sol, p, g) = @views 0.5*p.invf*parsevalsum2(sol[2], g) / (g.Lx*g.Ly)

"""
    qgke(prob)

Returns the QG kinetic energy.
"""
qgke(prob) = qgke(prob.sol, prob.vars, prob.params, prob.grid)

function qgke(sol, v, p, g)
  @views @. v.qh = sol[1]
  @views ldiv!(v.phi, g.fftplan, sol[2])
  @views calczetah!(v.zetah, v.qh, sol[2], v.phi, v, p, g)
  @. v.Uh = g.invKrsq * abs2(v.zetah)
  1/(g.Lx*g.Ly) * 0.5 * parsevalsum(v.Uh, g)
end

"""
    wavepe(prob)

Returns the potential energ of the near-inertial waves.
"""
wavepe(prob) = wavepe(prob.sol, prob.vars, prob.params, prob.grid)

function wavepe(phih, v, p, g)
  @. v.phixh = g.k*phih
  @. v.phiyh = g.l*phih
  1/(g.Lx*g.Ly) * 0.25*p.eta*p.invf * (parsevalsum2(v.phixh, g) + parsevalsum2(v.phiyh, g))
end

"""
    coupledenergy(prob)

Returns the 'coupled energy' in the NIW-QG flow.
"""
coupledenergy(prob) = qgke(prob) + wavepe(prob)

end # module
