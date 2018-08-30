!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This module wraps the MPI include statement because on certain platforms
!! just writing "USE MPI" does not work.
MODULE NTMPIModule
#if USE_MPIH
  include "mpif.h"
#else
  USE MPI
#endif
END MODULE
