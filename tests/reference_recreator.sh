#!/bin/bash
# REFERENCE DATA CREATOR
DIR=${1:-.}
for path in $(find $DIR -type f -name "*.msh"); do
  dir=$(dirname $path)
  mesh=$(basename $path)
  cnt=${mesh%.msh}.cnt
  res=${mesh%.msh}.res
  rm -fr $dir/${mesh}_psf*
  rm -fr $dir/${res}*
  [ ! -e $dir/$cnt ] && continue
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
  docker run -it --sig-proxy=false --rm -u $UID -v $PWD:$PWD -w $PWD registry.gitlab.com/frontistr-commons/frontistr/fistr1:master fistr1 -t 1
  rm -f hecmw_ctrl.dat FSTR.msg 0.log FSTR.sta FSTR.dbg.0 hecmw_vis.ini dyna*.txt
  popd
done


