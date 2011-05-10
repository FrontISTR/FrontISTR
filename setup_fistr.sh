#!/bin/sh

REMOVEMAKEFILES=0

BUILDTARGET="build-default"

SETUPFILE="setup_fistr.sh"
USER_CONFIGFILE="Makefile.conf"
HECMW_CONFIGFILE="Makefile.dev"
INTERMED_CONFIGFILE="Makefile.mid"
MAKEFILE_NAME="Makefile"
MAKEFILE_SETUPFILE="Makefile.am"

LIBSRCDIRS="\
    fistr"

BUILDDIRS="${LIBSRCDIRS} ."

for i in $*
do
	if [ "\"$i\"" = "\"-remove-makefiles\"" -o "\"$i\"" = "\"--remove-makefiles\"" ]; then
		REMOVEMAKEFILES=1
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

#------------------------------------------------------------------------------#
#
# create Makefile
#
for i in ${BUILDDIRS}
do
	sed	-e "s!@prefix@!${PREFIX}!" \
		-e "s!@bindir@!${BINDIR}!" \
		-e "s!@libdir@!${LIBDIR}!" \
		-e "s!@includedir@!${INCLUDEDIR}!" \
		-e "s!@cc@!${CC}!" \
		-e "s!@cflags@!${CFLAGS}!" \
		-e "s!@base_cflags@!${BASE_CFLAGS}!" \
		-e "s!@ldflags@!${LDFLAGS}!" \
		-e "s!@optflags@!${OPTFLAGS}!" \
		-e "s!@make@!${MAKE}!" \
		-e "s!@ar@!${AR}!" \
		-e "s!@cp@!${CP}!" \
		-e "s!@rm@!${RM}!" \
		-e "s!@ranlib@!${RANLIB}!" \
		-e "s!@mkdir@!${MKDIR}!" \
		-e "s!@fstrexec_targetfile@!${FSTREXEC_TARGETFILE}!" \
		-e "s!@cobjfilepostfix@!${COBJFILEPOSTFIX}!" \
		-e "s!@build_target@!${BUILDTARGET}!" \
		$i/${MAKEFILE_SETUPFILE} > $i/${MAKEFILE_NAME}
done

${RM} ${INTERMED_CONFIGFILE}
