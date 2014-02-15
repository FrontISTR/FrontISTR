#!/bin/sh

SERIAL=1
DEBUGMODE=0
REMOVEMAKEFILES=0
GATHERMAKEFILES=0
WITHTOOLS=0
WITHRCAP=0
WITHREFINER=0

BUILDTARGET="build-default"
NOBUILDTARGET="no-build"
BUILDTARGET_SERIAL="build-serial"
TOOLSTARGET="build-tools"
BUILDTARGET_RCAP="build-with-rcap"
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
    src/hecmw2 \
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
    src/main"

BUILDDIRS="${LIBSRCDIRS} ."

for i in $*
do
	if [ "\"$i\"" = "\"-g\"" -o "\"$i\"" = "\"-debug\"" -o "\"$i\"" = "\"--debug\"" ]; then
		DEBUGMODE=1
	elif [ "\"$i\"" = "\"-p\"" -o "\"$i\"" = "\"-parallel\"" -o "\"$i\"" = "\"--parallel\"" ]; then
		SERIAL=0
	elif [ "\"$i\"" = "\"-with-tools\"" -o "\"$i\"" = "\"--with-tools\"" ]; then
		WITHTOOLS=1
	elif [ "\"$i\"" = "\"-with-revocap\"" -o "\"$i\"" = "\"--with-revocap\"" ]; then
		WITHRCAP=1
	elif [ "\"$i\"" = "\"-with-refiner\"" -o "\"$i\"" = "\"--with-refiner\"" ]; then
		WITHREFINER=1
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
			--with-revocap          link revocap
			--with-refiner          link refiner
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
			--with-revocap          link revocap
			--with-refiner          link refiner
			-h, --help              show help(this message)
		EOF
		exit 1
	#else
	#	echo "Unknown paramer: " $i " (ignored, -h:help)"
	fi
done

#------------------------------------------------------------------------------#
#
# create intermidiate config file
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
# with revocap
#
if [ ${WITHRCAP} -eq 1 ]; then
	BUILDTARGET=${BUILDTARGET_RCAP}
	ALLBUILDTARGET=${BUILDTARGET_RCAP}
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
		-e "s!@hecmwlibs2@!${HECMWLIBS2}!" \
		-e "s!@hecmw_cflags@!${HECMW_CFLAGS}!" \
		-e "s!@hecmw_ldflags@!${HECMW_LDFLAGS}!" \
		-e "s!@hecmw_cppflags@!${HECMW_CPPFLAGS}!" \
		-e "s!@hecmw_cppldflags@!${HECMW_CPPLDFLAGS}!" \
		-e "s!@hecmw_f90flags@!${HECMW_F90FLAGS}!" \
		-e "s!@hecmw_f90ldflags@!${HECMW_F90LDFLAGS}!" \
		-e "s!@hecmw_cflags2@!${HECMW_CFLAGS2}!" \
		-e "s!@hecmw_ldflags2@!${HECMW_LDFLAGS2}!" \
		-e "s!@hecmw_cppflags2@!${HECMW_CPPFLAGS2}!" \
		-e "s!@hecmw_cppldflags2@!${HECMW_CPPLDFLAGS2}!" \
		-e "s!@hecmw_f90flags2@!${HECMW_F90FLAGS2}!" \
		-e "s!@hecmw_f90ldflags2@!${HECMW_F90LDFLAGS2}!" \
		-e "s!@fstrlibs@!${FSTRLIBS}!" \
		-e "s!@fstr_cflags@!${FSTR_CFLAGS}!" \
		-e "s!@fstr_ldflags@!${FSTR_LDFLAGS}!" \
		-e "s!@fstr_cppflags@!${FSTR_CPPFLAGS}!" \
		-e "s!@fstr_cppldflags@!${FSTR_CPPLDFLAGS}!" \
		-e "s!@fstr_f90flags@!${FSTR_F90FLAGS}!" \
		-e "s!@fstr_f90ldflags@!${FSTR_F90LDFLAGS}!" \
		-e "s!@build_target@!${BUILDTARGET}!" \
		-e "s!@all_build_target@!${ALLBUILDTARGET}!" \
		-e "s!@revocap_f90flags@!${REVOCAP_F90FLAGS}!" \
		-e "s!@revocap_f90ldflags@!${REVOCAP_F90LDFLAGS}!" \
		-e "s!@refiner_cflags@!${REFINER_CFLAGS}!" \
		-e "s!@refiner_ldflags@!${REFINER_LDFLAGS}!" \
		$i/${MAKEFILE_SETUPFILE} > $i/${MAKEFILE_NAME}
done

${RM} ${INTERMED_CONFIGFILE}
