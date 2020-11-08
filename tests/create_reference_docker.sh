#!/bin/bash
#
# Reference result creator using Docker
#
# Requirements
# -------------
# - `docker` command works well, e.g. `docker ps` does not returns error
#
# Input environment variable
# ---------------------------
# - REFERENCE_IMAGE : Docker image tag which will be used for generating reference result
#

: ${REFERENCE_IMAGE:=registry.gitlab.com/frontistr-commons/frontistr/fistr1:master}

echo_err () {
  echo -e "\e[31m$1\e[m" >&2
}

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
  docker run --rm        \
    -u $(id -u):$(id -g) \
    -v $PWD:$PWD         \
    -w $PWD              \
    ${REFERENCE_IMAGE}   \
    fistr1 -t 1
  rm -f hecmw_ctrl.dat FSTR.msg 0.log FSTR.sta FSTR.dbg.0 hecmw_vis.ini dyna*.txt
  popd
done

