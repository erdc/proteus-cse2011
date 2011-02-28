from pyadh.default_so import *
import sloshbox
pnList = [("curvature_sloshbox_2d_p" , 
           "curvature_sloshbox_2d_c0p1_n"),
          ("twp_navier_stokes_sloshbox_2d_p" , 
           "twp_navier_stokes_sloshbox_2d_c0p1c0p1_n"),
          ("ls_sloshbox_2d_p" , "ls_sloshbox_2d_c0p1_n"),
          ("vof_sloshbox_2d_p" , "vof_sloshbox_2d_c0p1_n"),
          ("redist_sloshbox_2d_p" , "redist_sloshbox_2d_c0p1_n"),
          ("ls_consrv_sloshbox_2d_p" , 
           "ls_consrv_sloshbox_2d_c0p1_n")]
name = "twp_navier_stokes_sloshbox_2d"
systemStepControllerType = Sequential_MinAdaptiveModelStep
needEBQ_GLOBAL = True
needEBQ = True
useOneArchive = False
archiveFlag = ArchiveFlags.EVERY_USER_STEP
tnList = [0.0,sloshbox.dt_init]+[sloshbox.dt_init+\
i*(sloshbox.T-sloshbox.dt_init)/(float(sloshbox.nDTout)-1) \
for i in range(1,sloshbox.nDTout)]
