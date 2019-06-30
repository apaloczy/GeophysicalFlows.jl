var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#GeophysicalFlows.jl-Documentation-1",
    "page": "Home",
    "title": "GeophysicalFlows.jl Documentation",
    "category": "section",
    "text": ""
},

{
    "location": "#Overview-1",
    "page": "Home",
    "title": "Overview",
    "category": "section",
    "text": "GeophysicalFlows.jl is a collection of modules which leverage the  FourierFlows.jl framework to provide solvers for problems in Geophysical Fluid Dynamics, on periodic domains and using Fourier-based pseudospectral methods."
},

{
    "location": "#Developers-1",
    "page": "Home",
    "title": "Developers",
    "category": "section",
    "text": "FourierFlows is currently being developed by Navid C. Constantinou and Gregory L. Wagner."
},

{
    "location": "#Cite-1",
    "page": "Home",
    "title": "Cite",
    "category": "section",
    "text": "The code is citable via zenodo."
},

{
    "location": "modules/twodturb/#",
    "page": "TwoDTurb Module",
    "title": "TwoDTurb Module",
    "category": "page",
    "text": ""
},

{
    "location": "modules/twodturb/#TwoDTurb-Module-1",
    "page": "TwoDTurb Module",
    "title": "TwoDTurb Module",
    "category": "section",
    "text": "newcommandJmathsfJ"
},

{
    "location": "modules/twodturb/#Basic-Equations-1",
    "page": "TwoDTurb Module",
    "title": "Basic Equations",
    "category": "section",
    "text": "This module solves two-dimensional incompressible turbulence. The flow is given through a streamfunction psi as (uupsilon) = (-partial_ypsi partial_xpsi). The dynamical variable used here is the component of the vorticity of the flow normal to the plane of motion, zeta=partial_x upsilon- partial_y u = nabla^2psi. The equation solved by the module is:partial_t zeta + J(psi zeta) = underbrace-leftmu(-1)^n_mu nabla^2n_mu\n+nu(-1)^n_nu nabla^2n_nuright zeta_textrmdissipation + f where J(a b) = (partial_x a)(partial_y b)-(partial_y a)(partial_x b). On the right hand side, f(xyt) is forcing, mu is hypoviscosity, and nu is hyperviscosity. Plain old linear drag corresponds to n_mu=0, while normal viscosity corresponds to n_nu=1."
},

{
    "location": "modules/twodturb/#Implementation-1",
    "page": "TwoDTurb Module",
    "title": "Implementation",
    "category": "section",
    "text": "The equation is time-stepped forward in Fourier space:partial_t widehatzeta = - widehatJ(psi zeta) -left(mu k^2n_mu\n+nu k^2n_nuright) widehatzeta  + widehatf In doing so the Jacobian is computed in the conservative form: J(ab) = partial_y  (partial_x a) b -partial_x (partial_y a) b.Thus:mathcalL = -mu k^-2n_mu - nu k^2n_nu mathcalN(widehatzeta) = - mathrmik_x mathrmFFT(u zeta)-\n	mathrmik_y mathrmFFT(upsilon zeta) + widehatf "
},

{
    "location": "modules/twodturb/#AbstractTypes-and-Functions-1",
    "page": "TwoDTurb Module",
    "title": "AbstractTypes and Functions",
    "category": "section",
    "text": "ParamsFor the unforced case (f=0) parameters AbstractType is build with Params and it includes:nu:   Float; viscosity or hyperviscosity coefficient.\nnnu: Integer0; the order of viscosity n_nu. Case n_nu=1 give normal viscosity.\nmu: Float; bottom drag or hypoviscosity coefficient.\nnmu: Integerge 0; the order of hypodrag n_mu. Case n_mu=0 give plain linear drag mu.For the forced case (fne 0) parameters AbstractType is build with ForcedParams. It includes all parameters in Params and additionally:calcF!: Function that calculates the forcing widehatfVarsFor the unforced case (f=0) variables AbstractType is build with Vars and it includes:zeta: Array of Floats; relative vorticity.\nu: Array of Floats; x-velocity, u.\nv: Array of Floats; y-velocity, upsilon.\nsol: Array of Complex; the solution, widehatzeta.\nzetah: Array of Complex; the Fourier transform widehatzeta.\nuh: Array of Complex; the Fourier transform widehatu.\nvh: Array of Complex; the Fourier transform widehatupsilon.For the forced case (fne 0) variables AbstractType is build with ForcedVars. It includes all variables in Vars and additionally:Fh: Array of Complex; the Fourier transform widehatf.\nprevsol: Array of Complex; the values of the solution sol at the previous time-step (useful for calculating the work done by the forcing).calcN! functionThe nonlinear term mathcalN(widehatzeta) is computed via functions:calcN_advection!: computes - widehatJ(psi zeta) and stores it in array N.\ncalcN_forced!: computes - widehatJ(psi zeta) via calcN_advection! and then adds to it the forcing widehatf computed via calcF! function. Also saves the solution widehatzeta of the previous time-step in array prevsol.\nupdatevars!: uses sol to compute zeta, u, upsilon, widehatu, and widehatupsilon and stores them into corresponding arrays of Vars/ForcedVars."
},

{
    "location": "modules/twodturb/#Examples-1",
    "page": "TwoDTurb Module",
    "title": "Examples",
    "category": "section",
    "text": "examples/twodturb_mcwilliams1984.jl: A script that simulates decaying two-dimensional turbulence reproducing the results of the paper by\nMcWilliams, J. C. (1984). The emergence of isolated coherent vortices in turbulent flow. J. Fluid Mech., 146, 21-43.\nexamples/twodturb_randomdecay.jl: A script that simulates decaying two-dimensional turbulence starting from random initial conditions.\nexamples/twodturb_stochasticforcing.jl: A script that simulates forced-dissipative two-dimensional turbulence with isotropic temporally delta-correlated stochastic forcing."
},

{
    "location": "modules/barotropicqg/#",
    "page": "BarotropicQG Module",
    "title": "BarotropicQG Module",
    "category": "page",
    "text": ""
},

{
    "location": "modules/barotropicqg/#BarotropicQG-Module-1",
    "page": "BarotropicQG Module",
    "title": "BarotropicQG Module",
    "category": "section",
    "text": "newcommandJmathsfJ"
},

{
    "location": "modules/barotropicqg/#Basic-Equations-1",
    "page": "BarotropicQG Module",
    "title": "Basic Equations",
    "category": "section",
    "text": "This module solves the quasi-geostrophic barotropic vorticity equation on a beta-plane of variable fluid depth H-h(xy). The flow is obtained through a streamfunction psi as (u upsilon) = (-partial_ypsi partial_xpsi). All flow fields can be obtained from the quasi-geostrophic potential vorticity (QGPV). Here the QGPV isunderbracef_0 + beta y_textplanetary PV + underbrace(partial_x upsilon\n	- partial_y u)_textrelative vorticity +\n	underbracefracf_0 hH_texttopographic PVThe dynamical variable is the component of the vorticity of the flow normal to the plane of motion, zetaequiv partial_x upsilon- partial_y u = nabla^2psi. Also, we denote the topographic PV with etaequiv f_0 hH. Thus, the equation solved by the module is:partial_t zeta + J(psi underbracezeta + eta_equiv q) +\nbetapartial_xpsi = underbrace-leftmu + nu(-1)^n_nu nabla^2n_nu\nright zeta _textrmdissipation + f where J(a b) = (partial_x a)(partial_y b)-(partial_y a)(partial_x b). On the right hand side, f(xyt) is forcing, mu is linear drag, and nu is hyperviscosity. Plain old viscosity corresponds to n_nu=1. The sum of relative vorticity and topographic PV is denoted with qequivzeta+eta."
},

{
    "location": "modules/barotropicqg/#Implementation-1",
    "page": "BarotropicQG Module",
    "title": "Implementation",
    "category": "section",
    "text": "The equation is time-stepped forward in Fourier space:partial_t widehatzeta = - widehatJ(psi q) +betafracmathrmik_xk^2widehatzeta -left(mu\n+nu k^2n_nuright) widehatzeta  + widehatf In doing so the Jacobian is computed in the conservative form: J(fg) = partial_y  (partial_x f) g -partial_x (partial_y f) g.Thus:mathcalL = betafracmathrmik_xk^2 - mu - nu k^2n_nu mathcalN(widehatzeta) = - mathrmik_x mathrmFFT(u q)-\n	mathrmik_y mathrmFFT(upsilon q) "
},

{
    "location": "modules/barotropicqg/#Examples-1",
    "page": "BarotropicQG Module",
    "title": "Examples",
    "category": "section",
    "text": "examples/barotropicqg/decayingbetaturb.jl: An script that simulates decaying quasi-geostrophic flow on a beta-plane demonstrating zonation.\nexamples/barotropicqg/forcedbetaturb.jl: An script that simulates forced-dissipative quasi-geostrophic flow on a beta-plane demonstrating zonation. The forcing is temporally delta-corraleted and its spatial structure is isotropic with power in a narrow annulus of total radius kf in wavenumber space.\nexamples/barotropicqg/ACConelayer.jl: A script that simulates barotropic quasi-geostrophic flow above topography reproducing the results of the paper by\nConstantinou, N. C. (2018). A barotropic model of eddy saturation. J. Phys. Oceanogr., 48 (2), 397-411."
},

{
    "location": "modules/multilayerqg/#",
    "page": "MultilayerQG Module",
    "title": "MultilayerQG Module",
    "category": "page",
    "text": ""
},

{
    "location": "modules/multilayerqg/#MultilayerQG-Module-1",
    "page": "MultilayerQG Module",
    "title": "MultilayerQG Module",
    "category": "section",
    "text": "newcommandJmathsfJ"
},

{
    "location": "modules/multilayerqg/#Basic-Equations-1",
    "page": "MultilayerQG Module",
    "title": "Basic Equations",
    "category": "section",
    "text": "This module solves the layered quasi-geostrophic equations on a beta-plane of variable fluid depth H-h(xy). The flow in each layer is obtained through a streamfunction psi_j as (u_j upsilon_j) = (-partial_ypsi_j partial_xpsi_j), j=1n, where n is the number of fluid layers.The QGPV in each layer ismathrmQGPV_j = q_j  + underbracef_0+beta y_textrmplanetary PV + delta_jnunderbracefracf_0 hH_n_textrmtopographic PVquad j=1nwhereq_1 = nabla^2psi_1 + F_32 1 (psi_2-psi_1)\nq_j = nabla^2psi_j + F_j-12 j (psi_j-1-psi_j) + F_j+12 j (psi_j+1-psi_j)quad j=2dotsn-1\nq_n = nabla^2psi_n + F_n-12 n (psi_n-1-psi_n)withF_j+12 k = fracf_0^2g_j+12 H_kquadtextandquad\ng_j+12 = gfracrho_j+1-rho_jrho_j+1 Therefore, in Fourier space the q\'s and psi\'s are related throughbeginpmatrix widehatq_boldsymbolk1vdotswidehatq_boldsymbolkn endpmatrix =\nunderbraceleft(-boldsymbolk^2mathbb1 + mathbbF right)_equiv mathbbS_boldsymbolk\nbeginpmatrix widehatpsi_boldsymbolk1vdotswidehatpsi_boldsymbolkn endpmatrixwheremathbbF equiv beginpmatrix\n -F_32 1               F_32 1     0     cdots     0\n  F_32 2  -(F_32 2+F_52 2)  F_52 2         vdots\n 0                             ddots     ddots    ddots  \n vdots                                                      0 \n 0                  cdots                0    F_n-12 n  -F_n-12 n\nendpmatrixIncluding an imposed zonal flow U_j(y) in each layer the equations of motion are:partial_t q_j + J(psi_j q_j ) + (U_j - partial_ypsi_j) partial_x Q_j +  U_j partial_x q_j  + (partial_y Q_j)(partial_xpsi_j) = -delta_jnmunabla^2psi_n - nu(-1)^n_nu nabla^2n_nu q_jwithpartial_y Q_j equiv beta - partial_y^2 U_j - (1-delta_j1)F_j-12 j (U_j-1-U_j) - (1-delta_jn)F_j+12 j (U_j+1-U_j) + delta_jnpartial_yeta \npartial_x Q_j equiv delta_jnpartial_xetaThe eddy kinetic energy in each layer is:textrmKE_j = dfracH_jH int dfrac12 boldsymbolnablapsi_j^2 fracmathrmd^2boldsymbolxL_x L_yquad j=1dotsnwhile the eddy potential energy related to each of fluid interface istextrmPE_j+12 = int dfrac12 dfracf_0^2g_j+12 (psi_j-psi_j+1)^2 fracmathrmd^2boldsymbolxL_x L_yquad j=1dotsn-1The lateral eddy fluxes in each layer are:textrmlateralfluxes_j = dfracH_jH int dfrac12 U_jupsilon_j partial_y u_j fracmathrmd^2boldsymbolxL_x L_yquad j=1dotsnwhile the vertical fluxes accros fluid interfaces are:textrmverticalfluxes_j+12 = int dfracf_0^2g_j+12 H (U_j-U_j+1)upsilon_j+1psi_j fracmathrmd^2boldsymbolxL_x L_yquad j=1dotsn-1\n"
},

{
    "location": "modules/multilayerqg/#Implementation-1",
    "page": "MultilayerQG Module",
    "title": "Implementation",
    "category": "section",
    "text": "Matrices mathbbS_boldsymbolk as well as mathbbS^-1_boldsymbolk are included in params as params.S and params.invS respectively.You can get widehatpsi_j from widehatq_j with streamfunctionfrompv!(psih, qh, invS, grid), while to go from widehatpsi_j back to widehatq_j pvfromstreamfunction!(qh, psih, S, grid).The equations are time-stepped forward in Fourier space:partial_t widehatq_j = - widehatJ(psi_j q_j)  - widehatU_j partial_x Q_j - widehatU_j partial_x q_j\n+ widehat(partial_ypsi_j) partial_x Q_j  - widehat(partial_xpsi_j)(partial_y Q_j) + delta_jnmu k^2 widehatpsi_n - nu k^2n_nu widehatq_jIn doing so the Jacobian is computed in the conservative form: J(fg) = partial_y  (partial_x f) g -partial_x (partial_y f) g.Thus:mathcalL = - nu k^2n_nu mathcalN(widehatq_j) = - widehatJ(psi_j q_j) - widehatU_j partial_x Q_j - widehatU_j partial_x q_j\n + widehat(partial_ypsi_j)(partial_x Q_j) - widehat(partial_xpsi_j)(partial_y Q_j) + delta_jnmu k^2 widehatpsi_n "
},

{
    "location": "man/types/#",
    "page": "Private types",
    "title": "Private types",
    "category": "page",
    "text": ""
},

{
    "location": "man/types/#Private-types-1",
    "page": "Private types",
    "title": "Private types",
    "category": "section",
    "text": ""
},

{
    "location": "man/types/#Private-types-in-module-GeophysicalFlows:-1",
    "page": "Private types",
    "title": "Private types in module GeophysicalFlows:",
    "category": "section",
    "text": "Modules = [GeophysicalFlows]\nPublic = false\nOrder = [:type]"
},

{
    "location": "man/types/#GeophysicalFlows.TwoDTurb.ForcedVars-Tuple{Any}",
    "page": "Private types",
    "title": "GeophysicalFlows.TwoDTurb.ForcedVars",
    "category": "method",
    "text": "ForcedVars(g)\n\nReturns the vars for forced two-dimensional turbulence with grid g.\n\n\n\n\n\n"
},

{
    "location": "man/types/#GeophysicalFlows.TwoDTurb.Params",
    "page": "Private types",
    "title": "GeophysicalFlows.TwoDTurb.Params",
    "category": "type",
    "text": "Params(nu, nnu, mu, nmu, calcF!)\n\nReturns the params for two-dimensional turbulence.\n\n\n\n\n\n"
},

{
    "location": "man/types/#GeophysicalFlows.TwoDTurb.StochasticForcedVars-Tuple{Any}",
    "page": "Private types",
    "title": "GeophysicalFlows.TwoDTurb.StochasticForcedVars",
    "category": "method",
    "text": "StochasticForcedVars(g; T)\n\nReturns the vars for stochastically forced two-dimensional turbulence with grid g.\n\n\n\n\n\n"
},

{
    "location": "man/types/#GeophysicalFlows.TwoDTurb.Vars-Tuple{Any}",
    "page": "Private types",
    "title": "GeophysicalFlows.TwoDTurb.Vars",
    "category": "method",
    "text": "Vars(g)\n\nReturns the vars for unforced two-dimensional turbulence with grid g.\n\n\n\n\n\n"
},

{
    "location": "man/types/#Private-types-in-module-TwoDTurb:-1",
    "page": "Private types",
    "title": "Private types in module TwoDTurb:",
    "category": "section",
    "text": "Modules = [GeophysicalFlows.TwoDTurb]\nPublic = false\nOrder = [:type]"
},

{
    "location": "man/types/#GeophysicalFlows.BarotropicQG.ForcedVars-Tuple{Any}",
    "page": "Private types",
    "title": "GeophysicalFlows.BarotropicQG.ForcedVars",
    "category": "method",
    "text": "ForcedVars(g)\n\nReturns the vars for forced two-dimensional barotropic QG problem with grid g.\n\n\n\n\n\n"
},

{
    "location": "man/types/#GeophysicalFlows.BarotropicQG.Params",
    "page": "Private types",
    "title": "GeophysicalFlows.BarotropicQG.Params",
    "category": "type",
    "text": "Params(g::TwoDGrid, f0, beta, FU, eta, mu, nu, nnu, calcFU, calcFq)\n\nReturns the params for an unforced two-dimensional barotropic QG problem.\n\n\n\n\n\n"
},

{
    "location": "man/types/#GeophysicalFlows.BarotropicQG.Params-Tuple{Any,Any,Any,Function,Any,Any,Any,Any,Any}",
    "page": "Private types",
    "title": "GeophysicalFlows.BarotropicQG.Params",
    "category": "method",
    "text": "Params(g::TwoDGrid, f0, beta, eta::Function, mu, nu, nnu, calcFU, calcFq)\n\nConstructor for Params that accepts a generating function for the topographic PV.\n\n\n\n\n\n"
},

{
    "location": "man/types/#GeophysicalFlows.BarotropicQG.StochasticForcedVars-Tuple{Any}",
    "page": "Private types",
    "title": "GeophysicalFlows.BarotropicQG.StochasticForcedVars",
    "category": "method",
    "text": "StochasticForcedVars(g)\n\nReturns the vars for stochastically forced two-dimensional barotropic QG problem with grid g.\n\n\n\n\n\n"
},

{
    "location": "man/types/#GeophysicalFlows.BarotropicQG.Vars-Tuple{Any}",
    "page": "Private types",
    "title": "GeophysicalFlows.BarotropicQG.Vars",
    "category": "method",
    "text": "Vars(g)\n\nReturns the vars for unforced two-dimensional barotropic QG problem with grid g.\n\n\n\n\n\n"
},

{
    "location": "man/types/#Private-types-in-module-BarotropicQG:-1",
    "page": "Private types",
    "title": "Private types in module BarotropicQG:",
    "category": "section",
    "text": "Modules = [GeophysicalFlows.BarotropicQG]\nPublic = false\nOrder = [:type]"
},

{
    "location": "man/functions/#",
    "page": "Functions",
    "title": "Functions",
    "category": "page",
    "text": ""
},

{
    "location": "man/functions/#Functions-1",
    "page": "Functions",
    "title": "Functions",
    "category": "section",
    "text": ""
},

{
    "location": "man/functions/#Functions-exported-from-GeophysicalFlows:-1",
    "page": "Functions",
    "title": "Functions exported from GeophysicalFlows:",
    "category": "section",
    "text": "Modules = [GeophysicalFlows]\nPrivate = false\nOrder = [:function]"
},

{
    "location": "man/functions/#GeophysicalFlows.TwoDTurb.Problem-Tuple{}",
    "page": "Functions",
    "title": "GeophysicalFlows.TwoDTurb.Problem",
    "category": "method",
    "text": "Problem(; parameters...)\n\nConstruct a 2D turbulence problem.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#GeophysicalFlows.TwoDTurb.dissipation-Tuple{Any}",
    "page": "Functions",
    "title": "GeophysicalFlows.TwoDTurb.dissipation",
    "category": "method",
    "text": "dissipation(prob)\n\nReturns the domain-averaged dissipation rate. nnu must be >= 1.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#GeophysicalFlows.TwoDTurb.drag-Tuple{Any}",
    "page": "Functions",
    "title": "GeophysicalFlows.TwoDTurb.drag",
    "category": "method",
    "text": "drag(prob)\n\nReturns the extraction of domain-averaged energy by drag/hypodrag mu.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#GeophysicalFlows.TwoDTurb.energy-Tuple{Any}",
    "page": "Functions",
    "title": "GeophysicalFlows.TwoDTurb.energy",
    "category": "method",
    "text": "energy(prob)\n\nReturns the domain-averaged kinetic energy in the Fourier-transformed vorticity solution sol.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#GeophysicalFlows.TwoDTurb.enstrophy-Tuple{Any}",
    "page": "Functions",
    "title": "GeophysicalFlows.TwoDTurb.enstrophy",
    "category": "method",
    "text": "enstrophy(prob)\n\nReturns the domain-averaged enstrophy in the Fourier-transformed vorticity solution sol.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#GeophysicalFlows.TwoDTurb.set_zeta!-Tuple{Any,Any}",
    "page": "Functions",
    "title": "GeophysicalFlows.TwoDTurb.set_zeta!",
    "category": "method",
    "text": "set_zeta!(prob, zeta)\n\nSet the solution sol as the transform of zeta and update variables v on the grid g.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#GeophysicalFlows.TwoDTurb.updatevars!-Tuple{Any}",
    "page": "Functions",
    "title": "GeophysicalFlows.TwoDTurb.updatevars!",
    "category": "method",
    "text": "updatevars!(prob)\n\nUpdate the vars in v on the grid g with the solution in sol.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#GeophysicalFlows.TwoDTurb.work-Tuple{Any,GeophysicalFlows.TwoDTurb.ForcedVars,Any}",
    "page": "Functions",
    "title": "GeophysicalFlows.TwoDTurb.work",
    "category": "method",
    "text": "work(prob)\nwork(sol, v, g)\n\nReturns the domain-averaged rate of work of energy by the forcing Fh.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#Functions-exported-from-TwoDTurb:-1",
    "page": "Functions",
    "title": "Functions exported from TwoDTurb:",
    "category": "section",
    "text": "Modules = [GeophysicalFlows.TwoDTurb]\nPrivate = false\nOrder = [:function]"
},

{
    "location": "man/functions/#GeophysicalFlows.BarotropicQG.dissipation-NTuple{4,Any}",
    "page": "Functions",
    "title": "GeophysicalFlows.BarotropicQG.dissipation",
    "category": "method",
    "text": "dissipation(prob)\ndissipation(s, v, p, g)\n\nReturns the domain-averaged dissipation rate. nnu must be >= 1.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#GeophysicalFlows.BarotropicQG.drag-NTuple{4,Any}",
    "page": "Functions",
    "title": "GeophysicalFlows.BarotropicQG.drag",
    "category": "method",
    "text": "drag(prob)\ndrag(s, v, p, g)\n\nReturns the extraction of domain-averaged energy by drag mu.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#GeophysicalFlows.BarotropicQG.energy-Tuple{Any}",
    "page": "Functions",
    "title": "GeophysicalFlows.BarotropicQG.energy",
    "category": "method",
    "text": "Calculate the domain-averaged kinetic energy.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#GeophysicalFlows.BarotropicQG.enstrophy-Tuple{Any}",
    "page": "Functions",
    "title": "GeophysicalFlows.BarotropicQG.enstrophy",
    "category": "method",
    "text": "Returns the domain-averaged enstrophy.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#GeophysicalFlows.BarotropicQG.meanenergy-Tuple{Any}",
    "page": "Functions",
    "title": "GeophysicalFlows.BarotropicQG.meanenergy",
    "category": "method",
    "text": "Returns the energy of the domain-averaged U.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#GeophysicalFlows.BarotropicQG.meanenstrophy-Tuple{Any}",
    "page": "Functions",
    "title": "GeophysicalFlows.BarotropicQG.meanenstrophy",
    "category": "method",
    "text": "Returns the enstrophy of the domain-averaged U.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#GeophysicalFlows.BarotropicQG.set_zeta!-Tuple{Any,GeophysicalFlows.BarotropicQG.Vars,Any,Any,Any}",
    "page": "Functions",
    "title": "GeophysicalFlows.BarotropicQG.set_zeta!",
    "category": "method",
    "text": "set_zeta!(prob, zeta)\nset_zeta!(s, v, g, zeta)\n\nSet the solution sol as the transform of zeta and update variables v on the grid g.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#GeophysicalFlows.BarotropicQG.updatevars!-NTuple{4,Any}",
    "page": "Functions",
    "title": "GeophysicalFlows.BarotropicQG.updatevars!",
    "category": "method",
    "text": "updatevars!(v, s, g)\n\nUpdate the vars in v on the grid g with the solution in sol.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#GeophysicalFlows.BarotropicQG.work-Tuple{Any,GeophysicalFlows.BarotropicQG.ForcedVars,Any}",
    "page": "Functions",
    "title": "GeophysicalFlows.BarotropicQG.work",
    "category": "method",
    "text": "work(prob)\nwork(s, v, p, g)\n\nReturns the domain-averaged rate of work of energy by the forcing Fqh.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#Functions-exported-from-BarotropicQG:-1",
    "page": "Functions",
    "title": "Functions exported from BarotropicQG:",
    "category": "section",
    "text": "Modules = [GeophysicalFlows.BarotropicQG]\nPrivate = false\nOrder = [:function]"
},

]}
