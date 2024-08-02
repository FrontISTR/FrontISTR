#!/bin/sh

#
# Flags
#
WITHFORTRAN=1
DEBUGMODE=0
WITHLEX=0
WITHMESSAGE=0
WITHMETIS=0
WITHPARMETIS=0
WITHTOOLS=0
REMOVEMAKEFILES=0
GATHERMAKEFILES=0
MESSAGEONLY=0
LEXONLY=0
SERIAL=1
WITHREFINER=0
WITHPARACON=0
WITHMUMPS=0
WITHMKL=0
WITHML=0
WITHLAPACK=0
OLDRESFORMAT=0

#
# Targets
#
ALLBUILDTARGET=""
BUILDTARGET="build-default"
NOBUILDTARGET="no-build"
CREATEINCLUDEDIRTARGET="create-include-dir"
CREATELIBDIRTARGET="create-lib-dir"
CREATEBINDIRTARGET="create-bin-dir"
LIBSTARGET="build-libs"
TOOLSTARGET="build-tools"
MESSAGETARGET="setup-msg"
LEXTARGET="setup-lex"

#
# Files
#
SETUPFILE="setup_hecmw.sh"
USER_CONFIGFILE="Makefile.conf"
HECMW_CONFIGFILE="Makefile.dev"
INTERMED_CONFIGFILE="Makefile.mid"
MAKEFILE_NAME="Makefile"
MAKEFILE_SETUPFILE="Makefile.am"
MAKEFILE_ARCHIVE="hecmw-makefiles.tar"

#
# Environment Variables
#
METIS_CFLAGS=""
METIS_LDFLAGS=""
PARMETIS_CFLAGS=""
PARMETIS_LDFLAGS=""

#
# Directories
#
UTILDIRS="\
	util"

LIBSRCDIRS="\
	src \
	src/common \
	src/operations \
	src/operations/adaptation \
	src/operations/dynamic_load_balancing \
	src/operations/element_smoothing \
	src/operations/jacobian \
	src/couple \
	src/solver \
	src/solver/matrix \
	src/solver/las \
	src/solver/iterative \
	src/solver/precond/bilu \
	src/solver/precond/diag \
	src/solver/precond/ml \
	src/solver/precond/rif \
	src/solver/precond/sainv \
	src/solver/precond/ssor \
	src/solver/precond \
	src/solver/solver_direct \
	src/solver/solver_direct_parallel \
	src/solver/solver_direct_lag \
	src/solver/sparse_matrix \
	src/solver/mumps \
	src/solver/mkl \
	src/solver/clustermkl \
	src/solver/communication \
	src/solver/init \
	src/solver/mpc \
	src/solver/main \
	src/solver/contact \
	src/visualizer \
	src/hecmw \
	src/etc"

TOOLSDIRS="\
	tools \
	tools/partitioner \
	tools/visualizer \
	tools/hec2rcap \
	tools/result_type_converter \
	tools/result_file_merger"

BUILDDIRS="${UTILDIRS} ${LIBSRCDIRS} ${TOOLSDIRS} ."

#------------------------------------------------------------------------------#
###
### Setup Options
###
for i in $*
do
	if [ "\"$i\"" = "\"-g\"" -o "\"$i\"" = "\"-debug\"" -o "\"$i\"" = "\"--debug\"" ]; then
		DEBUGMODE=1
	elif [ "\"$i\"" = "\"-p\"" -o "\"$i\"" = "\"-parallel\"" -o "\"$i\"" = "\"--parallel\"" ]; then
		SERIAL=0
	elif [ "\"$i\"" = "\"-l\"" -o "\"$i\"" = "\"-with-lex\"" -o "\"$i\"" = "\"--with-lex\"" ]; then
		WITHLEX=1
	elif [ "\"$i\"" = "\"-m\"" -o "\"$i\"" = "\"-with-message\"" -o "\"$i\"" = "\"--with-message\"" ]; then
		WITHMESSAGE=1
	elif [ "\"$i\"" = "\"-with-metis\"" -o "\"$i\"" = "\"--with-metis\"" ]; then
		WITHMETIS=1
	elif [ "\"$i\"" = "\"-with-parmetis\"" -o "\"$i\"" = "\"--with-parmetis\"" ]; then
		WITHPARMETIS=1
	elif [ "\"$i\"" = "\"-with-tools\"" -o "\"$i\"" = "\"--with-tools\"" ]; then
		WITHTOOLS=1
	elif [ "\"$i\"" = "\"-with-refiner\"" -o "\"$i\"" = "\"--with-refiner\"" ]; then
		WITHREFINER=1
	elif [ "\"$i\"" = "\"-with-paracon\"" -o "\"$i\"" = "\"--with-paracon\"" ]; then
		WITHPARACON=1
	elif [ "\"$i\"" = "\"-with-mumps\"" -o "\"$i\"" = "\"--with-mumps\"" ]; then
		WITHMUMPS=1
	elif [ "\"$i\"" = "\"-with-mkl\"" -o "\"$i\"" = "\"--with-mkl\"" ]; then
		WITHMKL=1
	elif [ "\"$i\"" = "\"-with-ml\"" -o "\"$i\"" = "\"--with-ml\"" ]; then
		WITHML=1
	elif [ "\"$i\"" = "\"-with-lapack\"" -o "\"$i\"" = "\"--with-lapack\"" ]; then
		WITHLAPACK=1
	elif [ "\"$i\"" = "\"-old-res-format\"" -o "\"$i\"" = "\"--old-res-format\"" ]; then
		OLDRESFORMAT=1
	elif [ "\"$i\"" = "\"-remove-makefiles\"" -o "\"$i\"" = "\"--remove-makefiles\"" ]; then
		REMOVEMAKEFILES=1
		GATHERMAKEFILES=0
		MESSAGEONLY=0
		LEXONLY=0
	elif [ "\"$i\"" = "\"-gather-makefiles\"" -o "\"$i\"" = "\"--gather-makefiles\"" ]; then
		REMOVEMAKEFILES=0
		GATHERMAKEFILES=1
		MESSAGEONLY=0
		LEXONLY=0
	elif [ "\"$i\"" = "\"-only-message\"" -o "\"$i\"" = "\"--only-message\"" ]; then
		REMOVEMAKEFILES=0
		GATHERMAKEFILES=0
		MESSAGEONLY=1
		LEXONLY=0
	elif [ "\"$i\"" = "\"-only-lex\"" -o "\"$i\"" = "\"--only-lex\"" ]; then
		REMOVEMAKEFILES=0
		GATHERMAKEFILES=0
		MESSAGEONLY=0
		LEXONLY=1
	elif [ "\"$i\"" = "\"-show-all-options\"" -o "\"$i\"" = "\"--show-all-options\"" ]; then
		cat 1>&2 <<- EOF
			Usage: setup_hecmw.sh [-options]
			-g, --debug             debug mode
			-p, --parallel          for parallel environment with MPI
			-l, --with-lex          perform lexical analysis
			-m, --with-message      build message file by perl
			--with-tools            compile tools
			--with-metis            compile with METIS
			--with-parmetis         compile with ParMETIS
			--with-refiner          compile with REVOCAP_Refiner
			--with-paracon          for parallel contact
			--with-mumps            compile with MUMPS
			--with-mkl              compile with MKL PARDISO
			--with-ml               compile with ML
			--with-lapack           compile with LAPACK
			--only-message          only create error message files
			--only-lex              only perform lexical analyzer
			--remove-makefiles      remove all MAKEFILEs
			--gather-makefiles      archive all MAKEFILEs
			--show-all-options      print all options (show this message)
		EOF
		exit 1
	elif [ "\"$i\"" = "\"-h\"" -o "\"$i\"" = "\"-help\"" -o "\"$i\"" = "\"--help\"" ]; then
		cat 1>&2 <<- EOF
			Usage: setup_hecmw.sh [-options]
			-g, --debug             debug mode
			-p, --parallel          for parallel environment with MPI
			--with-tools            compile tools/ (build library only)
			--with-metis            compile with METIS
			--with-parmetis         compile with ParMETIS
			--with-refiner          compile with REVOCAP_Refiner
			--with-paracon          for parallel contact
			--with-mumps            compile with MUMPS
			--with-mkl              compile with MKL PARDISO
			--with-ml               compile with ML
			--with-lapack           compile with LAPACK
			-h, --help              show help (this message)
		EOF
		exit 1
	#else
	#	echo "Unknown parameter: " $i " (ignored, -h:help)"
	fi
done

#------------------------------------------------------------------------------#
###
### User Defined Environment Variables
###
sed -e "s!\([[:alnum:]_]\)[[:blank:]]*=[[:blank:]]*\(.*\)!\1='\2'!g" \
	${HECMW_CONFIGFILE} > ${INTERMED_CONFIGFILE}
sed -e "s!\([[:alnum:]_]\)[[:blank:]]*=[[:blank:]]*\(.*\)!\1='\2'!g" \
	${USER_CONFIGFILE} >> ${INTERMED_CONFIGFILE}

. ./${INTERMED_CONFIGFILE}

#------------------------------------------------------------------------------#
###
### Maintenance Options
###
#
# Gather MAKEFILEs and archive them
#
if [ ${GATHERMAKEFILES} -eq 1 ]; then
	tar cvf ${MAKEFILEARCHIVE} ./${USER_CONFIGFILE} ./${HECMW_CONFIGFILE} ./${SETUPFILE}
	for i in ${BUILDDIRS}
	do
		tar rvf ${MAKEFILE_ARCHIVE} $i/${MAKEFILE_SETUPFILE}
	done
	${RM} ${INTERMED_CONFIGFILE}
	exit 0
fi

#
# Remove MAKEFILEs
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
# Create message file
#
if [ ${MESSAGEONLY} -eq 1 ]; then
	ALLBUILDTARGET="${MESSAGETARGET}"
fi

#
# Perform lexical analyzer
#
if [ ${LEXONLY} -eq 1 ]; then
	ALLBUILDTARGET="${LEXTARGET}"
fi

#------------------------------------------------------------------------------#
###
### General Options
###
if [ ${MESSAGEONLY} -eq 0 -a ${LEXONLY} -eq 0 ]; then

	ALLBUILDTARGET="${CREATEINCLUDEDIRTARGET} ${CREATELIBDIRTARGET}"

	#
	# debug mode
	#
	if [ ${DEBUGMODE} -eq 0 ]; then
		OPTFLAGS="${OPTFLAGS} ${NODEBUG_OPTFLAGS}"
		F90OPTFLAGS="${F90OPTFLAGS}"
	else
		OPTFLAGS="${OPTFLAGS} ${DEBUGFLAGS} ${DEBUG_OPTFLAGS}"
		F90OPTFLAGS="${F90OPTFLAGS} ${F90DEBUGFLAGS}"
	fi

	#
	# with METIS / with ParMETIS
	#
	if [ ${WITHMETIS} -eq 1 ]; then
		F90FLAGS="${F90FLAGS} -DHECMW_WITH_METIS -DHECMW_METIS_VER=${HECMW_METIS_VER}"
	else
		METISDIR=""
		METISLIBDIR=""
		METISINCDIR=""
		METISLIBS=""
		METIS_CFLAGS=""
		METIS_LDFLAGS=""
		PARTITIONER_OPTFLAGS=""
	fi

	if [ ${WITHPARMETIS} -eq 0 ]; then
		PARMETISDIR=""
		PARMETISLIBDIR=""
		PARMETISINCDIR=""
		PARMETISLIBS=""
		PARMETIS_CFLAGS=""
		PARMETIS_LDFLAGS=""
	fi

	#
	# create message files
	#
	if [ ${WITHMESSAGE} -eq 1 ]; then
		ALLBUILDTARGET="${ALLBUILDTARGET} ${MESSAGETARGET}"
	fi

	#
	# perform flex
	#
	if [ ${WITHLEX} -eq 1 ]; then
		ALLBUILDTARGET="${ALLBUILDTARGET} ${LEXTARGET}"
	fi

	#
	# with TOOLS
	#
	if [ ${WITHTOOLS} -eq 1 ]; then
		ALLBUILDTARGET="${ALLBUILDTARGET} ${CREATEBINDIRTARGET} ${LIBSTARGET} ${TOOLSTARGET}"
	else
		ALLBUILDTARGET="${ALLBUILDTARGET} ${LIBSTARGET}"
	fi

	#
	# with FORTRAN and SERIAL setting
	#
	if [ ${WITHFORTRAN} -eq 1 ]; then
		HECMWLIBS="-lfhecmw -lhecmw"
		BUILDTARGET="build-default"
		if [ ${SERIAL} -eq 1 ]; then
			OPTFLAGS="${OPTFLAGS} ${SERIAL_OPTFLAGS}"
			F90OPTFLAGS="${F90OPTFLAGS} ${SERIAL_OPTFLAGS}"
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
		ALLBUILDTARGET="${ALLBUILDTARGET} ${BUILDTARGET}"
	else
		HECMWLIBS="-lhecmw"
		BUILDTARGET="build-without-f"
		ALLBUILDTARGET="${ALLBUILDTARGET} ${BUILDTARGET}"
		if [ ${SERIAL} -eq 1 ]; then
			OPTFLAGS="${OPTFLAGS} ${SERIAL_OPTFLAGS}"
			F90OPTFLAGS="${F90OPTFLAGS} ${SERIAL_OPTFLAGS}"
			MPI_CFLAGS=""
			MPI_LDFLAGS=""
			MPI_F90FLAGS=""
			MPI_F90LDFLAGS=""
		fi
	fi

	#
	# with Refiner
	#
	if [ ${WITHREFINER} -eq 1 ]; then
		HECMW_CFLAGS="${HECMW_CFLAGS} ${REFINER_CFLAGS}"
		HECMW_F90FLAGS="${HECMW_F90FLAGS} ${REFINER_CFLAGS}"
		HECMW_LDFLAGS="${HECMW_LDFLAGS} ${REFINER_LDFLAGS}"
		HECMW_F90LDFLAGS="${HECMW_F90LDFLAGS} ${REFINER_LDFLAGS}"
	fi

	#
	# with paracon
	#
	if [ ${WITHPARACON} -eq 1 ]; then
		HECMW_CFLAGS="${HECMW_CFLAGS} -DPARA_CONTACT"
	fi

	#
	# with MUMPS
	#
	if [ ${WITHMUMPS} -eq 1 ]; then
		F90FLAGS="${F90FLAGS} -DHECMW_WITH_MUMPS"
	else
		MUMPS_CFLAGS=""
		MUMPS_LDFLAGS=""
		MUMPS_F90FLAGS=""
		MUMPS_F90LDFLAGS=""
	fi

	#
	# with MKL PARDISO
	#
	if [ ${WITHMKL} -eq 1 ]; then
		F90FLAGS="${F90FLAGS} -DHECMW_WITH_MKL"
	else
		MKL_CFLAGS=""
		MKL_LDFLAGS=""
		MKL_F90FLAGS=""
		MKL_F90LDFLAGS=""
	fi

	#
	# with ML
	#
	if [ ${WITHML} -ne 1 ]; then
		ML_CFLAGS=""
		ML_LDFLAGS=""
		ML_F90FLAGS=""
		ML_F90LDFLAGS=""
	fi

	#
	# with LAPACK
	#
	if [ ${WITHLAPACK} -eq 1 ]; then
		F90FLAGS="${F90FLAGS} -DHECMW_WITH_LAPACK"
	fi

	#
	# old res format
	#
	if [ ${OLDRESFORMAT} -eq 1 ]; then
		HECMW_CFLAGS="${HECMW_CFLAGS} -DOLD_RES_FORMAT"
	fi

	#
	# C linker
	#
	if [ "${CLINKER}" = "" ]; then
		CLINKER="${CC}"
	fi
fi

#------------------------------------------------------------------------------#
###
### Create MAKEFILES
###
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
		-e "s!@cc@!${CC}!" \
		-e "s!@cflags@!${CFLAGS}!" \
		-e "s!@ldflags@!${LDFLAGS}!" \
		-e "s!@optflags@!${OPTFLAGS}!" \
		-e "s!@clinker@!${CLINKER}!" \
		-e "s!@cpp@!${CPP}!" \
		-e "s!@cppflags@!${CPPFLAGS}!" \
		-e "s!@cppldflags@!${CPPLDFLAGS}!" \
		-e "s!@cppoptflags@!${CPPOPTFLAGS}!" \
		-e "s!@f90@!${F90}!" \
		-e "s!@f90flags@!${F90FLAGS}!" \
		-e "s!@f90ldflags@!${F90LDFLAGS}!" \
		-e "s!@f90fpp@!${F90FPP}!" \
		-e "s!@f90optflags@!${F90OPTFLAGS}!" \
		-e "s!@make@!${MAKE}!" \
		-e "s!@ar@!${AR}!" \
		-e "s!@cp@!${CP}!" \
		-e "s!@rm@!${RM}!" \
		-e "s!@mv@!${MV}!" \
		-e "s!@ranlib@!${RANLIB}!" \
		-e "s!@mkdir@!${MKDIR}!" \
		-e "s!@lex@!${LEX}!" \
		-e "s!@hecmwlib_targetfile@!${HECMWLIB_TARGETFILE}!" \
		-e "s!@hecmwlib_f90targetfile@!${HECMWLIB_F90TARGETFILE}!" \
		-e "s!@partitioner_targetfile@!${PARTITIONER_TARGETFILE}!" \
		-e "s!@visualizer_targetfile@!${VISUALIZER_TARGETFILE}!" \
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
		-e "s!@metisdir@!${METISDIR}!" \
		-e "s!@metislibdir@!${METISLIBDIR}!" \
		-e "s!@metisincdir@!${METISINCDIR}!" \
		-e "s!@metislibs@!${METISLIBS}!" \
		-e "s!@hecmw_metis_ver@!${HECMW_METIS_VER}!" \
		-e "s!@metis_cflags@!${METIS_CFLAGS}!" \
		-e "s!@metis_ldflags@!${METIS_LDFLAGS}!" \
		-e "s!@metis_f90flags@!${METIS_F90FLAGS}!" \
		-e "s!@metis_f90ldflags@!${METIS_F90LDFLAGS}!" \
		-e "s!@parmetisdir@!${PARMETISDIR}!" \
		-e "s!@parmetislibdir@!${PARMETISLIBDIR}!" \
		-e "s!@parmetisincdir@!${PARMETISINCDIR}!" \
		-e "s!@parmetislibs@!${PARMETISLIBS}!" \
		-e "s!@parmetis_cflags@!${PARMETIS_CFLAGS}!" \
		-e "s!@parmetis_ldflags@!${PARMETIS_LDFLAGS}!" \
		-e "s!@parmetis_f90flags@!${PARMETIS_F90FLAGS}!" \
		-e "s!@parmetis_f90ldflags@!${PARMETIS_F90LDFLAGS}!" \
		-e "s!@partitioner_cflags@!${PARTITIONER_CFLAGS}!" \
		-e "s!@partitioner_ldflags@!${PARTITIONER_LDFLAGS}!" \
		-e "s!@partitioner_optflags@!${PARTITIONER_OPTFLAGS}!" \
		-e "s!@partitioner_f90flags@!${PARTITIONER_F90FLAGS}!" \
		-e "s!@partitioner_f90ldflags@!${PARTITIONER_F90LDFLAGS}!" \
		-e "s!@partitioner_f90optflags@!${PARTITIONER_F90OPTFLAGS}!" \
		-e "s!@visualizer_cflags@!${VISUALIZER_CFLAGS}!" \
		-e "s!@visualizer_ldflags@!${VISUALIZER_LDFLAGS}!" \
		-e "s!@visualizer_f90flags@!${VISUALIZER_F90FLAGS}!" \
		-e "s!@visualizer_f90ldflags@!${VISUALIZER_F90LDFLAGS}!" \
		-e "s!@refinerdir@!${REFINERDIR}!" \
		-e "s!@refinerincdir@!${REFINERINCDIR}!" \
		-e "s!@refinerlibdir@!${REFINERLIBDIR}!" \
		-e "s!@refinerlibs@!${REFINERLIBS}!" \
		-e "s!@mumpsdir@!${MUMPSDIR}!" \
		-e "s!@mumpslibdir@!${MUMPSLIBDIR}!" \
		-e "s!@mumpsincdir@!${MUMPSINCDIR}!" \
		-e "s!@mumpslibs@!${MUMPSLIBS}!" \
		-e "s!@mumps_cflags@!${MUMPS_CFLAGS}!" \
		-e "s!@mumps_ldflags@!${MUMPS_LDFLAGS}!" \
		-e "s!@mumps_f90flags@!${MUMPS_F90FLAGS}!" \
		-e "s!@mumps_f90ldflags@!${MUMPS_F90LDFLAGS}!" \
		-e "s!@mkldir@!${MKLDIR}!" \
		-e "s!@mkllibdir@!${MKLLIBDIR}!" \
		-e "s!@mklincdir@!${MKLINCDIR}!" \
		-e "s!@mkllibs@!${MKLLIBS}!" \
		-e "s!@mkl_cflags@!${MKL_CFLAGS}!" \
		-e "s!@mkl_ldflags@!${MKL_LDFLAGS}!" \
		-e "s!@mkl_f90flags@!${MKL_F90FLAGS}!" \
		-e "s!@mkl_f90ldflags@!${MKL_F90LDFLAGS}!" \
		-e "s!@mldir@!${MLDIR}!" \
		-e "s!@mllibdir@!${MLLIBDIR}!" \
		-e "s!@mlincdir@!${MLINCDIR}!" \
		-e "s!@mllibs@!${MLLIBS}!" \
		-e "s!@ml_cflags@!${ML_CFLAGS}!" \
		-e "s!@ml_ldflags@!${ML_LDFLAGS}!" \
		-e "s!@ml_f90flags@!${ML_F90FLAGS}!" \
		-e "s!@ml_f90ldflags@!${ML_F90LDFLAGS}!" \
		-e "s!@all_build_target@!${ALLBUILDTARGET}!" \
		-e "s!@build_target@!${BUILDTARGET}!" \
		$i/${MAKEFILE_SETUPFILE} > $i/${MAKEFILE_NAME}
done

${RM} ${INTERMED_CONFIGFILE}
