# Metis
METISDIR       = $(HOME)/work/ParMetis-3.2.0
METISLIBDIR    = $(METISDIR)
#METISINCDIR    = $(METISDIR)/METISLib
#METISLIBS      = -L$(METISLIBDIR) -lmetis

# ParMetis
PARMETISDIR    = $(HOME)/work/ParMetis-3.2.0
PARMETISLIBDIR = $(PARMETISDIR)
#PARMETISINCDIR = $(PARMETISDIR)/ParMETISLib
PARMETISLIBS   = -L$(PARMETISLIBDIR) -lparmetis

# Scotch
SCOTCHDIR      = $(HOME)/work/scotch_5.1.12_esmumps
#SCOTCHINCDIR   = $(SCOTCHDIR)/include
SCOTCHLIBDIR   = $(SCOTCHDIR)/lib
SCOTCHLIBS     = -L$(SCOTCHLIBDIR) -lptesmumps -lptscotch -lptscotcherr

# BLAS and ScaLAPACK
MKLDIR         = /opt/intel/composerxe/mkl
MKLLIBDIR      = $(MKLDIR)/lib/intel64
BLAS           = -L$(MKLLIBDIR) -Wl,--start-group \
               $(MKLLIBDIR)/libmkl_intel_lp64.a \
               $(MKLLIBDIR)/libmkl_sequential.a \
               $(MKLLIBDIR)/libmkl_core.a -Wl,--end-group
#BLAS           = -L$(MKLLIBDIR) -Wl,--start-group \
#               $(MKLLIBDIR)/libmkl_intel_lp64.a \
#               $(MKLLIBDIR)/libmkl_intel_thread.a \
#               $(MKLLIBDIR)/libmkl_core.a -Wl,--end-group
#BLAS           = -L$(MKLLIBDIR) -Wl,--start-group \
#               $(MKLLIBDIR)/libmkl_gf_lp64.a \
#               $(MKLLIBDIR)/libmkl_sequential.a \
#               $(MKLLIBDIR)/libmkl_core.a -Wl,--end-group
#BLAS           = -L$(MKLLIBDIR) -Wl,--start-group \
#               $(MKLLIBDIR)/libmkl_gf_lp64.a \
#               $(MKLLIBDIR)/libmkl_gnu_thread.a \
#               $(MKLLIBDIR)/libmkl_core.a -Wl,--end-group
#BLAS           = -L$(MKLLIBDIR) \
#               -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
#BLAS           = -L$(MKLLIBDIR) \
#               -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core
SCALAP         = -L$(MKLLIBDIR) -Wl,--start-group \
               $(MKLLIBDIR)/libmkl_scalapack_lp64.a \
               $(MKLLIBDIR)/libmkl_blacs_openmpi_lp64.a \
               -Wl,--end-group
#SCALAP         = -L$(MKLLIBDIR) \
#               -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64

# MUMPS
MUMPSDIR       = $(HOME)/work/MUMPS_4.10.0
MUMPSINC       = -I$(MUMPSDIR)/include
PORDLIBS       = -L$(MUMPSDIR)/lib -lpord
MUMPSLIBS      = -L$(MUMPSDIR)/lib -ldmumps -lmumps_common \
               $(PARMETISLIBS) $(METISLIBS) $(PORDLIBS) $(SCOTCHLIBS) \
	       $(SCALAP) $(BLAS)

# Fortran compiler settings
F90FLAGS_EXTRA = $(MUMPSINC)
#F90LDFLAGS_EXTRA = $(MUMPSLIBS) -lgfortran -L/opt/intel/lib/intel64 -liomp5 -lpthread -lm
F90LDFLAGS_EXTRA = $(MUMPSLIBS)
