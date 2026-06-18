#!/bin/bash
#
# Reference result creator using built fistr1
#

echo_err () {
  ESC=$(printf '\033')
  echo "${ESC}[31m$1${ESC}[m" >&2
}

check_executable () {
  if [[ ! -x $1 ]]; then
    echo_err "$1 is not executable"
    exit 1
  fi
}

usage () {
  cat <<EOM
Usage: $(basename "$0") [OPTION]...
  -h          Display this help
  -f VALUE    fistr1 binary path       (Default: ../build/fistr1/fistr1)
EOM
}

fistr1=`pwd`/../build/fistr1/fistr1
errors=0
while getopts ":d:p:t:f:e:r:h" optKey; do
  case "$optKey" in
    f)
      fistr1=${OPTARG};;
    h)
      usage; exit 0;;
    *)
      usage; exit 1;;
  esac
done

check_executable $fistr1

for path in $(find . -not -path "./_archive/*" -type f -name "*.msh" ); do
  dir=$(dirname $path)
  mesh=$(basename $path)
  cnt=${mesh%.msh}.cnt
  res=${mesh%.msh}.res

  if ls $dir/${res}.* > /dev/null 2>&1; then
    echo_err ${res}" exists. Skip this mesh file."
    continue
  fi
  # Remove previous results if exists
  rm -fr $dir/${mesh}_psf*
  rm -fr $dir/${res}*

  if [ ! -e $dir/$cnt ]; then
    echo_err "*.cnt file for $path is not found. Skip this mesh file."
    continue
  fi

  pushd $dir
  # frequency-response: regenerate the eigenmode fixtures first (mode-superposition input)
  if [ -e ${mesh%.msh}_eigen.cnt ]; then
    rm -f ${mesh%.msh}_eig.res.* ${mesh%.msh}_eig.log
cat <<EOLE > hecmw_ctrl.dat
!MESH, NAME=fstrMSH,TYPE=HECMW-ENTIRE
${mesh}
!CONTROL,NAME=fstrCNT
${mesh%.msh}_eigen.cnt
!RESULT,NAME=fstrRES,IO=OUT
${mesh%.msh}_eig.res
EOLE
    $fistr1 -t 1
    mv 0.log ${mesh%.msh}_eig.log
  fi
cat <<EOL > hecmw_ctrl.dat
!MESH, NAME=fstrMSH,TYPE=HECMW-ENTIRE
${mesh}
!CONTROL,NAME=fstrCNT
${cnt}
!RESULT,NAME=fstrRES,IO=OUT
${res}
!RESULT,NAME=vis_out,IO=OUT
${mesh}

EOL
  # frequency-response: map eigenmode fixtures (result-in) and dynamic result (fstrDYNA)
  if [ -e ${mesh%.msh}_eigen.cnt ]; then
    printf '!RESULT,NAME=result-in,IO=IN\n%s\n' "${mesh%.msh}_eig.res" >> hecmw_ctrl.dat
    printf '!RESULT,NAME=fstrDYNA,IO=OUT\n%s\n' "${mesh%.msh}_dyna.res" >> hecmw_ctrl.dat
  fi
  $fistr1 -t 1
  rm -f hecmw_ctrl.dat FSTR.msg 0.log FSTR.sta FSTR.dbg.0 hecmw_vis.ini dyna*.txt
  rm -fr ${mesh%.msh}_dyna.res* ${mesh%.msh}.vis*
  popd
done
