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

for path in $(find . -type f -name "*.msh"); do
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
  $fistr1 -t 1
  rm -f hecmw_ctrl.dat FSTR.msg 0.log FSTR.sta FSTR.dbg.0 hecmw_vis.ini dyna*.txt
  popd
done
