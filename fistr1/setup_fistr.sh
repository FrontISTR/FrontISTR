#!/bin/sh

SERIAL=1
DEBUGMODE=0
REMOVEMAKEFILES=0
GATHERMAKEFILES=0
WITHMETIS=0
WITHPARMETIS=0
WITHTOOLS=0
WITHRCAP=0
WITHREFINER=0
WITHMKL=0
WITHMUMPS=0
WITHML=0

BUILDTARGET="build-default"
NOBUILDTARGET="no-build"
BUILDTARGET_SERIAL="build-serial"
TOOLSTARGET="build-tools"
ALLBUILDTARGET=""

SETUPFILE="setup_fistr.sh"
USER_CONFIGFILE="Makefile.conf"
HECMW_CONFIGFILE="Makefile.dev"
INTERMED_CONFIGFILE="Makefile.mid"
MAKEFILE_NAME="Makefile"
MAKEFILE_SETUPFILE="Makefile.am"
MAKEFILEARCHIVE="FSTR-hecmw-makefiles.tar"

LIBSRCDIRS="\
    src \
    src/common \
    src/lib \
    src/lib/utilities \
    src/lib/element \
    src/lib/physics \
    src/lib/contact \
    src/lib/user \
    src/analysis \
    src/analysis/static \
    src/analysis/heat \
    src/analysis/dynamic \
    src/analysis/dynamic/mode \
    src/analysis/dynamic/transit \
    src/analysis/dynamic/freq \
    src/main \
    tools \
    tools/neu2fstr \
    tools/neu2fstr/converter \
    tools/neu2fstr/HECD \
    tools/neu2fstr/NFD"

BUILDDIRS="${LIBSRCDIRS} ."

for i in $*
do
	if [ "\"$i\"" = "\"-g\"" -o "\"$i\"" = "\"-debug\"" -o "\"$i\"" = "\"--debug\"" ]; then
		DEBUGMODE=1
	elif [ "\"$i\"" = "\"-p\"" -o "\"$i\"" = "\"-parallel\"" -o "\"$i\"" = "\"--parallel\"" ]; then
		SERIAL=0
	elif [ "\"$i\"" = "\"-with-metis\"" -o "\"$i\"" = "\"--with-metis\"" ]; then
		WITHMETIS=1
	elif [ "\"$i\"" = "\"-with-parmetis\"" -o "\"$i\"" = "\"--with-parmetis\"" ]; then
		WITHPARMETIS=1
	elif [ "\"$i\"" = "\"-with-tools\"" -o "\"$i\"" = "\"--with-tools\"" ]; then
		WITHTOOLS=1
	elif [ "\"$i\"" = "\"-with-revocap\"" -o "\"$i\"" = "\"--with-revocap\"" ]; then
		WITHRCAP=1
	elif [ "\"$i\"" = "\"-with-refiner\"" -o "\"$i\"" = "\"--with-refiner\"" ]; then
		WITHREFINER=1
	elif [ "\"$i\"" = "\"-with-mkl\"" -o "\"$i\"" = "\"--with-mkl\"" ]; then
		WITHMKL=1
	elif [ "\"$i\"" = "\"-with-mumps\"" -o "\"$i\"" = "\"--with-mumps\"" ]; then
		WITHMUMPS=1
	elif [ "\"$i\"" = "\"-with-ml\"" -o "\"$i\"" = "\"--with-ml\"" ]; then
		WITHML=1
	elif [ "\"$i\"" = "\"-remove-makefiles\"" -o "\"$i\"" = "\"--remove-makefiles\"" ]; then
		REMOVEMAKEFILES=1
	elif [ "\"$i\"" = "\"-gather-makefiles\"" -o "\"$i\"" = "\"--gather-makefiles\"" ]; then
		GATHERMAKEFILES=1
	elif [ "\"$i\"" = "\"-show-all-options\"" -o "\"$i\"" = "\"--show-all-options\"" ]; then
		cat 1>&2 <<- EOF
			Usage: setup_fstr.sh [-options]
			-g, --debug             debug mode
			-p, --parallel          for parallel environment with MPI
			--with-tools            compile tools
			--with-metis            compile with METIS
			--with-parmetis         compile with ParMETIS
			--with-revocap          link revocap
			--with-refiner          link refiner
			--with-mkl              link mkl
			--with-mumps            link mumps
			--with-ml               link ml
			--remove-makefiles      remove all MAKEFILEs
			--gather-makefiles      archive all MAKEFILEs
			--show-all-options      print all options (show this message)
		EOF
		exit 1
	elif [ "\"$i\"" = "\"-h\"" -o "\"$i\"" = "\"-help\"" -o "\"$i\"" = "\"--help\"" ]; then
		cat 1>&2 <<- EOF
			Usage: setup_fstr.sh [-options]
			-g, --debug             debug mode
			-p, --parallel          for parallel environment with MPI
			--with-tools            compile tools
			--with-metis            compile with METIS
			--with-parmetis         compile with ParMETIS
			--with-revocap          link revocap
			--with-refiner          link refiner
			--with-mkl              link mkl
			--with-mumps            link mumps
			--with-ml               link ml
			-h, --help              show help(this message)
		EOF
		exit 1
	#else
	#	echo "Unknown parameter: " $i " (ignored, -h:help)"
	fi
done

#------------------------------------------------------------------------------#
#
# create intermediate config file
#
sed -e "s!\([[:alnum:]_]\)[[:blank:]]*=[[:blank:]]*\(.*\)!\1='\2'!g" \
	${HECMW_CONFIGFILE} > ${INTERMED_CONFIGFILE}
sed -e "s!\([[:alnum:]_]\)[[:blank:]]*=[[:blank:]]*\(.*\)!\1='\2'!g" \
	${USER_CONFIGFILE} >> ${INTERMED_CONFIGFILE}

. ./${INTERMED_CONFIGFILE}

#------------------------------------------------------------------------------#
#
# gather Makefiles
#
if [ ${GATHERMAKEFILES} -eq 1 ]; then
	tar cvf ${MAKEFILEARCHIVE} ./${USER_CONFIGFILE} ./${HECMW_CONFIGFILE} ./${SETUPFILE}
	for i in ${BUILDDIRS}
	do
		tar rvf ${MAKEFILEARCHIVE} $i/${MAKEFILE_SETUPFILE}
	done
	${RM} ${INTERMED_CONFIGFILE}
	exit 0
fi

#
# remove Makefiles
#
if [ ${REMOVEMAKEFILES} -eq 1 ]; then
	for i in ${BUILDDIRS}
	do
		${RM} $i/${MAKEFILE_NAME}
	done
	${RM} ${INTERMED_CONFIGFILE}
	exit 0
fi

#
# debug mode
#
if [ ${DEBUGMODE} -eq 0 ]; then
	OPTFLAGS="${OPTFLAGS} ${NODEBUGOPTFLAGS}"
	F90OPTFLAGS="${F90OPTFLAGS}"
else
	OPTFLAGS="${OPTFLAGS} ${DEBUGFLAGS} ${DEBUGOPTFLAGS}"
	F90OPTFLAGS="${F90OPTFLAGS} ${F90DEBUGFLAGS}"
fi

#
# serial
#
if [ ${SERIAL} -eq 1 ];  then
	CFLAGS="${CFLAGS} -DHECMW_SERIAL"
	MPI_CFLAGS=""
	MPI_LDFLAGS=""
	MPI_F90FLAGS=""
	MPI_F90LDFLAGS=""
else
	if [ -z ${MPIDIR} ]; then
		MPI_CFLAGS=""
		MPI_LDFLAGS=""
		MPI_F90FLAGS=""
		MPI_F90LDFLAGS=""
	fi
fi

#
# with metis
#
if [ ${WITHMETIS} -eq 0 ]; then
	METIS_CFLAGS=""
	METIS_LDFLAGS=""
	METIS_F90FLAGS=""
	METIS_F90LDFLAGS=""
	HECMW_METIS_VER="0"
fi

#
# with parmetis
#
if [ ${WITHPARMETIS} -eq 0 ]; then
	PARMETIS_CFLAGS=""
	PARMETIS_LDFLAGS=""
	PARMETIS_F90FLAGS=""
	PARMETIS_F90LDFLAGS=""
fi

#
# with revocap
#
if [ ${WITHRCAP} -eq 1 ]; then
	F90FLAGS="${F90FLAGS} -DWITH_REVOCAP"
	ALLBUILDTARGET=${BUILDTARGET}
else
	ALLBUILDTARGET=${BUILDTARGET}
	REVOCAP_F90FLAGS=""
	REVOCAP_F90LDFLAGS=""
fi

#
# with tools
#
if [ ${WITHTOOLS} -eq 1 ]; then
	ALLBUILDTARGET="${ALLBUILDTARGET} ${TOOLSTARGET}"
fi

#
# with refiner
#
if [ ${WITHREFINER} -eq 1 ]; then
	HECMW_LDFLAGS="${HECMW_LDFLAGS} ${REFINER_LDFLAGS}"
	HECMW_F90LDFLAGS="${HECMW_F90LDFLAGS} ${REFINER_LDFLAGS}"
fi

#
# with mkl
#
if [ ${WITHMKL} -eq 1 ]; then
	F90FLAGS="${F90FLAGS} -DWITH_MKL"
fi

#
# with mumps
#
if [ ${WITHMUMPS} -eq 0 ]; then
	MUMPS_CFLAGS=""
	MUMPS_LDFLAGS=""
	MUMPS_F90FLAGS=""
	MUMPS_F90LDFLAGS=""
fi

#
# with ml
#
if [ ${WITHML} -ne 1 ]; then
	ML_LDFLAGS=""
	ML_F90LDFLAGS=""
fi

#
# F90 linker
#
if [ "${F90LINKER}" = "" ]; then
	F90LINKER="${F90}"
fi

#------------------------------------------------------------------------------#
#
# create Makefile
#
for i in ${BUILDDIRS}
do
	sed -e "s!@mpidir@!${MPIDIR}!" \
		-e "s!@mpibindir@!${MPIBINDIR}!" \
		-e "s!@mpilibdir@!${MPILIBDIR}!" \
		-e "s!@mpiincdir@!${MPIINCDIR}!" \
		-e "s!@mpilibs@!${MPILIBS}!" \
		-e "s!@prefix@!${PREFIX}!" \
		-e "s!@bindir@!${BINDIR}!" \
		-e "s!@libdir@!${LIBDIR}!" \
		-e "s!@includedir@!${INCLUDEDIR}!" \
		-e "s!@revocapdir@!${REVOCAPDIR}!" \
		-e "s!@revocapincdir@!${REVOCAPINCDIR}!" \
		-e "s!@revocaplibdir@!${REVOCAPLIBDIR}!" \
		-e "s!@revocaplibs@!${REVOCAPLIBS}!" \
		-e "s!@refinerdir@!${REFINERDIR}!" \
		-e "s!@refinerincdir@!${REFINERINCDIR}!" \
		-e "s!@refinerlibdir@!${REFINERLIBDIR}!" \
		-e "s!@refinerlibs@!${REFINERLIBS}!" \
		-e "s!@cc@!${CC}!" \
		-e "s!@cflags@!${CFLAGS}!" \
		-e "s!@ldflags@!${LDFLAGS}!" \
		-e "s!@optflags@!${OPTFLAGS}!" \
		-e "s!@cpp@!${CPP}!" \
		-e "s!@cppflags@!${CPPFLAGS}!" \
		-e "s!@cppldflags@!${CPPLDFLAGS}!" \
		-e "s!@cppoptflags@!${CPPOPTFLAGS}!" \
		-e "s!@f90@!${F90}!" \
		-e "s!@f90flags@!${F90FLAGS}!" \
		-e "s!@f90ldflags@!${F90LDFLAGS}!" \
		-e "s!@f90fpp@!${F90FPP}!" \
		-e "s!@f90optflags@!${F90OPTFLAGS}!" \
		-e "s!@f90linker@!${F90LINKER}!" \
		-e "s!@make@!${MAKE}!" \
		-e "s!@ar@!${AR}!" \
		-e "s!@cp@!${CP}!" \
		-e "s!@rm@!${RM}!" \
		-e "s!@ranlib@!${RANLIB}!" \
		-e "s!@mkdir@!${MKDIR}!" \
		-e "s!@fstrexec_targetfile@!${FSTREXEC_TARGETFILE}!" \
		-e "s!@fstrlib_targetfile@!${FSTRLIB_TARGETFILE}!" \
		-e "s!@fstrlib_f90targetfile@!${FSTRLIB_F90TARGETFILE}!" \
		-e "s!@f90modulepostfix@!${F90MODULEPOSTFIX}!" \
		-e "s!@cobjfilepostfix@!${COBJFILEPOSTFIX}!" \
		-e "s!@cppobjfilepostfix@!${CPPOBJFILEPOSTFIX}!" \
		-e "s!@f90objfilepostfix@!${F90OBJFILEPOSTFIX}!" \
		-e "s!@base_cflags@!${BASE_CFLAGS}!" \
		-e "s!@base_cppflags@!${BASE_CPPFLAGS}!" \
		-e "s!@base_f90flags@!${BASE_F90FLAGS}!" \
		-e "s!@debug_optflags@!${DEBUG_OPTFLAGS}!" \
		-e "s!@nodebug_optflags@!${NODEBUG_OPTFLAGS}!" \
		-e "s!@mpi_cflags@!${MPI_CFLAGS}!" \
		-e "s!@mpi_ldflags@!${MPI_LDFLAGS}!" \
		-e "s!@mpi_f90flags@!${MPI_F90FLAGS}!" \
		-e "s!@mpi_f90ldflags@!${MPI_F90LDFLAGS}!" \
		-e "s!@hecmwlibs@!${HECMWLIBS}!" \
		-e "s!@hecmw_cflags@!${HECMW_CFLAGS}!" \
		-e "s!@hecmw_ldflags@!${HECMW_LDFLAGS}!" \
		-e "s!@hecmw_cppflags@!${HECMW_CPPFLAGS}!" \
		-e "s!@hecmw_cppldflags@!${HECMW_CPPLDFLAGS}!" \
		-e "s!@hecmw_f90flags@!${HECMW_F90FLAGS}!" \
		-e "s!@hecmw_f90ldflags@!${HECMW_F90LDFLAGS}!" \
		-e "s!@fstrlibs@!${FSTRLIBS}!" \
		-e "s!@fstr_cflags@!${FSTR_CFLAGS}!" \
		-e "s!@fstr_ldflags@!${FSTR_LDFLAGS}!" \
		-e "s!@fstr_cppflags@!${FSTR_CPPFLAGS}!" \
		-e "s!@fstr_cppldflags@!${FSTR_CPPLDFLAGS}!" \
		-e "s!@fstr_f90flags@!${FSTR_F90FLAGS}!" \
		-e "s!@fstr_f90ldflags@!${FSTR_F90LDFLAGS}!" \
		-e "s!@metisdir@!${METISDIR}!" \
		-e "s!@metislibdir@!${METISLIBDIR}!" \
		-e "s!@metisincdir@!${METISINCDIR}!" \
		-e "s!@metislibs@!${METISLIBS}!" \
		-e "s!@metis_cflags@!${METIS_CFLAGS}!" \
		-e "s!@metis_ldflags@!${METIS_LDFLAGS}!" \
		-e "s!@metis_f90flags@!${METIS_F90FLAGS}!" \
		-e "s!@metis_f90ldflags@!${METIS_F90LDFLAGS}!" \
		-e "s!@hecmw_metis_ver@!${HECMW_METIS_VER}!" \
		-e "s!@parmetisdir@!${PARMETISDIR}!" \
		-e "s!@parmetislibdir@!${PARMETISLIBDIR}!" \
		-e "s!@parmetisincdir@!${PARMETISINCDIR}!" \
		-e "s!@parmetislibs@!${PARMETISLIBS}!" \
		-e "s!@parmetis_cflags@!${PARMETIS_CFLAGS}!" \
		-e "s!@parmetis_ldflags@!${PARMETIS_LDFLAGS}!" \
		-e "s!@parmetis_f90flags@!${PARMETIS_F90FLAGS}!" \
		-e "s!@parmetis_f90ldflags@!${PARMETIS_F90LDFLAGS}!" \
		-e "s!@mumpsdir@!${MUMPSDIR}!" \
		-e "s!@mumpslibdir@!${MUMPSLIBDIR}!" \
		-e "s!@mumpsincdir@!${MUMPSINCDIR}!" \
		-e "s!@mumpslibs@!${MUMPSLIBS}!" \
		-e "s!@mumps_cflags@!${MUMPS_CFLAGS}!" \
		-e "s!@mumps_ldflags@!${MUMPS_LDFLAGS}!" \
		-e "s!@mumps_f90flags@!${MUMPS_F90FLAGS}!" \
		-e "s!@mumps_f90ldflags@!${MUMPS_F90LDFLAGS}!" \
		-e "s!@mldir@!${MLDIR}!" \
		-e "s!@mllibdir@!${MLLIBDIR}!" \
		-e "s!@mllibs@!${MLLIBS}!" \
		-e "s!@ml_ldflags@!${ML_LDFLAGS}!" \
		-e "s!@ml_f90ldflags@!${ML_F90LDFLAGS}!" \
		-e "s!@mkldir@!${MKLDIR}!" \
		-e "s!@mkllibdir@!${MKLLIBDIR}!" \
		-e "s!@mklincdir@!${MKLINCDIR}!" \
		-e "s!@mkllibs@!${MKLLIBS}!" \
		-e "s!@mkl_cflags@!${MKL_CFLAGS}!" \
		-e "s!@mkl_ldflags@!${MKL_LDFLAGS}!" \
		-e "s!@mkl_f90flags@!${MKL_F90FLAGS}!" \
		-e "s!@mkl_f90ldflags@!${MKL_F90LDFLAGS}!" \
		-e "s!@build_target@!${BUILDTARGET}!" \
		-e "s!@all_build_target@!${ALLBUILDTARGET}!" \
		-e "s!@revocap_f90flags@!${REVOCAP_F90FLAGS}!" \
		-e "s!@revocap_f90ldflags@!${REVOCAP_F90LDFLAGS}!" \
		-e "s!@refiner_cflags@!${REFINER_CFLAGS}!" \
		-e "s!@refiner_ldflags@!${REFINER_LDFLAGS}!" \
		$i/${MAKEFILE_SETUPFILE} > $i/${MAKEFILE_NAME}
done

${RM} ${INTERMED_CONFIGFILE}
