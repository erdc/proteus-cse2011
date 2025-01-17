from pyadh import *
from pyadh.default_p import *
from sloshbox import *
"""
Two-phase incompressible Navier-Stokes flow of air and water in a perturbed box
"""

## \page Tests Test Problems
#\ref twp_navier_stokes_ls_so_sloshbox_2d_p.py "Two-phase incompressible Navier-Stokes flow of air and water in a perturbed box"
# 

##\ingroup test
# \file twp_navier_stokes_ls_so_sloshbox_2d_p.py
# @{
#
# \brief Two-phase incompressible Navier-Stokes flow of air and water in a perturbed box
#
# The governing equations are described by the
# pyadh::TransportCoefficients::TwophaseNavierStokes_LS_SO class. The
# domain is the unit square. The boundary conditions are slip and no
# flow on the whole boundary. The pressure is set (arbitrarily) to
# zero at the node in the top right corner. The free surface is
# initially described by the line
#\f[
#(\frac{1}{2}-x_0)\sin(\theta)+(y-y_0)\cos(\theta) = 0
#\f]
#where \f$\tan(\theta)\f$ is the slope of the free surface and and \f$(x_0,y_0)\f$ is a point on the free surface
# 
#\image html sloshbox_1.jpg
#\image html sloshbox_2.jpg
#
# <A href="https://adh.usace.army.mil/pyadh-images/sloshbox.avi"> AVI Animation of velocity and pressure </A>
#
# <A href="https://adh.usace.army.mil/pyadh-images/sloshbox_adh.avi"> AVI Animation of Richard Stockstill's moving mesh simulation with ADH</A>
#

                     
analyticalSolution = None

coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=epsFact_viscosity,
                                             sigma=sigma_01,
                                             rho_0 = rho_0,
                                             nu_0 = nu_0,
                                             rho_1 = rho_1,
                                             nu_1 = nu_1,
                                             g=g,
                                             nd=nd,
                                             LS_model=2,
                                             KN_model=0,
                                             epsFact_density=epsFact_density,
                                             stokes=useStokes)
# coefficients = TwophaseNavierStokes_LS_SO(epsFact=epsFact_viscosity,
#                                              rho_0 = rho_0,
#                                              nu_0 = nu_0,
#                                              rho_1 = rho_1,
#                                              nu_1 = nu_1,
#                                              g=g,
#                                              nd=nd,
#                                              LS_model=2)
# if VOF:
#     coefficients = TwophaseNavierStokes_VOF_SO(epsFact=epsFact,
#                                                rho_0 = rho_0,
#                                                nu_0 = nu_0,
#                                                rho_1 = rho_1,
#                                                nu_1 = nu_1,
#                                                g=g,
#                                                nd=nd)
# else:
#     coefficients = TwophaseNavierStokes_LS_SO(epsFact=epsFact,
#                                               rho_0 = rho_0,
#                                               nu_0 = nu_0,
#                                               rho_1 = rho_1,
#                                               nu_1 = nu_1,
#                                               g=g,
#                                               nd=nd,
#                                               LS_model=1)
#coefficients = TwophaseStokes_LS_SO(g=[0.0,9.8],nd=nd,steady=True)
#coefficients.rho_1 = coefficients.rho_0
#coefficients.mu_1 = coefficients.mu_0

