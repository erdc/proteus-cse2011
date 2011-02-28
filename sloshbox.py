from math import *
import pyadh.MeshTools
applyCorrection=True
applyRedistancing=True
#rdtimeIntegration='newton'
#freezeLevelSet=False
rdtimeIntegration='osher'
freezeLevelSet=True#False
sloshbox_quad_order = 3
nodalPartitioning=False#True
nd = 2
useBackwardEuler=True
useBackwardEuler_ls=True

dt_init=1.0e-4
T=20.0
#nDTout = None
nDTout=200
runCFL = 0.33

lag_ns_subgridError=True
lag_ls_shockCapturing=True

ns_shockCapturingFactor=0.9
ls_shockCapturingFactor=0.9
vof_shockCapturingFactor=0.9
rd_shockCapturingFactor=0.9

#epsilons for Heaviside/Dirac/etc smoothing
epsFact_density = 1.5
epsFact_viscosity = 1.5
epsFact_redistance = 0.33
epsFact_curvature=3.0
epsFact_consrv_heaviside=1.5
epsFact_consrv_dirac=1.5
epsFact_consrv_diffusion=10.0
epsFact_vof=1.5

usePETSc=True#False

spaceOrder=1

L = (10.0,2.0,1.0)

useStokes=False

#nnx=201
#nny=41
nnx=21
nny=5
he = L[0]/float(nnx-1)
nLevels = 1#

useShock=True
shock_x = 0.5*L[0]
shock_y = 0.5*L[1]

waterLevel = 0.5*L[1]
slopeAngle = 0.5*(pi/2.0)

#water
rho_0=998.2
nu_0=1.004e-6
#air
rho_1=1.205
nu_1= 1.500e-5 #* 1000.0
sigma_01=72.8e-3
#gravity
g=[0.0,-9.8,0.0]
VOF = False
closedTop = False#True
sidesNoSlip=False

def getDBC_p_sloshbox(x,flag):
    if closedTop:#set pressure  in top right corner
       if x[1] == L[1]:
           if x[0] == L[0]:
               return lambda x,t: 0.0
    else:
        if x[1] == L[1]:
            return lambda x,t: 0.0

def getDBC_u_sloshbox(x,flag):
    return None

def getDBC_v_sloshbox(x,flag):
    return None

def getDBC_w_sloshbox(x,flag):
    return None

dirichletConditions = {0:getDBC_p_sloshbox,
                       1:getDBC_u_sloshbox,
                       2:getDBC_v_sloshbox}

def getAFBC_p_sloshbox(x,flag):
    if closedTop:
        if (x[0] == 0.0 or
            x[0] == L[0] or
            x[1] == 0.0 or
            x[1] == L[1]):
            return lambda x,t: 0.0
    else:
        if (x[0] == 0.0 or
            x[0] == L[0] or
            x[1] == 0.0):
            return lambda x,t: 0.0

def getAFBC_u_sloshbox(x,flag):
    if closedTop:
        if (x[0] == 0.0 or
            x[0] == L[0] or
            x[1] == 0.0 or
            x[1] == L[1]):
            return lambda x,t: 0.0
    else:
        if (x[0] == 0.0 or
            x[0] == L[0] or
            x[1] == 0.0):
            return lambda x,t: 0.0
        
def getAFBC_v_sloshbox(x,flag):
    if closedTop:
        if (x[0] == 0.0 or
            x[0] == L[0] or
            x[1] == 0.0 or
            x[1] == L[1]):
            return lambda x,t: 0.0
    else:
        if (x[0] == 0.0 or
            x[0] == L[0] or
            x[1] == 0.0):
            return lambda x,t: 0.0

def getDFBC_u_sloshbox(x,flag):
    if (x[0] == 0.0 or
        x[0] == L[0] or
        x[1] == 0.0 or
        x[1] == L[1]):
        return lambda x,t: 0.0
    
def getDFBC_v_sloshbox(x,flag):
    if (x[0] == 0.0 or
        x[0] == L[0] or
        x[1] == 0.0 or
        x[1] == L[1]):
        return lambda x,t: 0.0

fluxBoundaryConditions = {0:'mixedFlow',
                          1:'mixedFlow',
                          2:'mixedFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_sloshbox,
                                    1:getAFBC_u_sloshbox,
                                    2:getAFBC_v_sloshbox}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u_sloshbox,2:getDFBC_u_sloshbox},
                                   2:{1:getDFBC_v_sloshbox,2:getDFBC_v_sloshbox}}

def shockSignedDistance(x):
    phi_x = shock_x - x[0]
    phi_y = x[1] - shock_y
    phi2 = sqrt((x[0]-shock_x)**2 + (x[1] - shock_y)**2)
    if x[0] < shock_x:
        if x[1] < shock_y:
            return phi_x
        else:
            return phi2
    else:
        if x[1] > shock_y:
            return phi_y
        elif fabs(phi_y) < fabs(phi_x):
            return phi_y
        else:
            return phi_x

class PerturbedSurface_p:
    def __init__(self,waterLevel,slopeAngle,shock_x,shock_y):
        self.waterLevel=waterLevel
        self.slopeAngle=slopeAngle
        self.shock_x=shock_x
        self.shock_y=shock_y
    def uOfXT(self,x,t):
        if useShock:
            if x[0] < shock_x:
                return -rho_1*g[1]*(L[1] - x[1])
            else:
                if x[1] > shock_y:
                    return -rho_1*g[1]*(L[1] - x[1])
                else:
                    return -rho_1*g[1]*(L[1] - shock_y)-rho_0*g[1]*(shock_y - x[1])
        else:
            surfaceNormal = [-sin(self.slopeAngle),cos(self.slopeAngle)]
            topRight_z= min(L[1],-(L[0]-0.5)*surfaceNormal[0] + self.waterLevel*surfaceNormal[1])
            signedDistance = (x[0] - 0.5)*surfaceNormal[0]+(x[1] - self.waterLevel)*surfaceNormal[1]
            signedDistanceCorner = (0.0 - 0.5)*surfaceNormal[0]+(L[1] - self.waterLevel)*surfaceNormal[1]
#         if signedDistance > 0.0:
#             p = -(L[1]-x[1])*rho_1*g[1]
#         else:
#             p = -(L[1]-topRight_z)*rho_1*g[1] - (topRight_z-x[1])*rho_0*g[1]
            if signedDistance > 0.0:
                p = -(signedDistanceCorner-signedDistance)*rho_1*g[1]
            else:
                p = -signedDistanceCorner*rho_1*g[1] + signedDistance*rho_0*g[1]
        return p

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:PerturbedSurface_p(waterLevel,slopeAngle,shock_x,shock_y),
                     1:AtRest(),
                     2:AtRest()}

if nodalPartitioning:
    restrictFineSolutionToAllMeshes=False
    parallelPartitioningType = pyadh.MeshTools.MeshParallelPartitioningTypes.node
    nLayersOfOverlapForParallel = 1
else:
    restrictFineSolutionToAllMeshes=False
    parallelPartitioningType = pyadh.MeshTools.MeshParallelPartitioningTypes.element
    nLayersOfOverlapForParallel = 2
restrictFineSolutionToAllMeshes=False
parallelPartitioningType = pyadh.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 2
