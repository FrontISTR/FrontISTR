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

write_hecmw_ctrl () {
  ctrl_cnt=$1
  ctrl_res=$2
cat <<EOL > hecmw_ctrl.dat
!MESH, NAME=fstrMSH,TYPE=HECMW-ENTIRE
${mesh}
!CONTROL,NAME=fstrCNT
${ctrl_cnt}
!RESULT,NAME=fstrRES,IO=OUT
${ctrl_res}
!RESULT,NAME=vis_out,IO=OUT
${mesh}

EOL
  if [ -e $cnt_eig ]; then
cat <<EOL >> hecmw_ctrl.dat
!RESULT,NAME=fstrDYNA,IO=OUT
${dres}
!RESULT,NAME=result-in,IO=IN
${eres}
EOL
  fi
}

for path in $(find . -not -path "./_archive/*" -type f -name "*.msh" ); do
  dir=$(dirname $path)
  mesh=$(basename $path)
  cnt=${mesh%.msh}.cnt
  cnt_eig=${mesh%.msh}_eigen.cnt
  res=${mesh%.msh}.res
  eres=${mesh%.msh}_eigen.res
  dres=${mesh%.msh}_dyna.res

  if ls $dir/${res}.* > /dev/null 2>&1; then
    echo_err ${res}" exists. Skip this mesh file."
    continue
  fi
  # Remove previous results if exists
  rm -fr $dir/${mesh}_psf*
  rm -fr $dir/${res}*
  rm -fr $dir/${dres}*

  if [ ! -e $dir/$cnt ]; then
    echo_err "*.cnt file for $path is not found. Skip this mesh file."
    continue
  fi

  pushd $dir
  #
  # Frequency response analysis reads the modes of a preceding eigen analysis,
  # so run it first and hand over its log file and result files
  #
  if [ -e $cnt_eig ]; then
    write_hecmw_ctrl $cnt_eig $eres
    $fistr1 -t 1
    mv 0.log eigen_0.log
  fi
  write_hecmw_ctrl $cnt $res
  $fistr1 -t 1
  rm -f hecmw_ctrl.dat FSTR.msg 0.log eigen_0.log FSTR.sta FSTR.dbg.0 hecmw_vis.ini dyna*.txt
  rm -f ${eres}.*
  popd
done
