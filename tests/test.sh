#!/bin/bash

function usage {
  cat <<EOM
Usage: $(basename "$0") [OPTION]...
  -h          Display help
  -d VALUE    Target test case dir     (Default: .)
  -p VALUE    MPI processes            (Default: 1)
  -t VALUE    OMP threads              (Default: 1)
  -f VALUE    fistr1 binary path       (Default: fistr1)
  -e VALUE    hecmw_part1 binary path  (Default: hecmw_part1)
  -r VALUE    rmerge binary path       (Default: rmerge)
  -v          Verbose test message
  -c          Go ahead if test failed
EOM
  exit 127
}

function die {
  exit $ERRORS
}


target=.
fistr1=fistr1
hecmw_part1=hecmw_part1
rmerge=rmerge
np=1
nt=1
VERBOSE=0
goahead=0
ERRORS=0
while getopts ":d:p:t:f:e:r:vch" optKey; do
  case "$optKey" in
    d)
      target=${OPTARG};;
    p)
      np=${OPTARG};;
    t)
      nt=${OPTARG};;
    f)
      fistr1=${OPTARG};;
    e)
      hecmw_part1=${OPTARG};;
    r)
      rmerge=${OPTARG};;
    v)
      VERBOSE=$(($VERBOSE+1));;
    c)
      goahead=1;;
    '-h'|'--help'|* )
      usage;;
  esac
done

test_dir=run_test/$$
output=""
compare_res=$(pwd)/compare_res.pl
rm -fr $test_dir
mkdir -p $test_dir
tdir=$(pwd)/$test_dir

for path in $(find $target -type f -name "*.msh"); do
  SECONDS=0
  cdir=$(cd $(dirname $path) && pwd)
  mesh=$(basename $path)
  cnt=${mesh%.msh}.cnt
  res=${mesh%.msh}.res
  [ ! -e $cdir/$cnt ] && continue
  echo $path:
  cp -r $cdir/$mesh $tdir
  cp -r $cdir/$cnt $tdir
  pushd $tdir > /dev/null
  [ $np -gt 1 ] && MESHTYPE=HECMW-DIST || MESHTYPE=HECMW-ENTIRE
cat <<EOL > hecmw_ctrl.dat
!MESH, NAME=fstrMSH,TYPE=$MESHTYPE
${mesh}
!CONTROL,NAME=fstrCNT
${cnt}
!RESULT,NAME=fstrRES,IO=OUT
${res}
!RESULT,NAME=vis_out,IO=OUT
${mesh}
!MESH, NAME=part_in,TYPE=HECMW-ENTIRE
${mesh}
!MESH, NAME=part_out,TYPE=HECMW-DIST
${mesh}
EOL
  echo "!PARTITION,TYPE=NODE-BASED,METHOD=KMETIS,DOMAIN=$np" > hecmw_part_ctrl.dat
  [ $VERBOSE -le 2 ] && OUT1=" >fistr1.log 2>&1" && OUT2=" >hecmw_part1.log 2>&1"
  if [ $np -gt 1 ]; then
    sh -c "$hecmw_part1 $OUT2"
    sh -c "mpirun --oversubscribe --allow-run-as-root -n $np $fistr1 -t $nt $OUT1"
  else
    sh -c "$fistr1 -t $nt $output $OUT1"
  fi
  COMPARE=0
  find . -name "*.res.*"|awk -F. '{print $NF}'|sort|uniq|xargs -I{} $rmerge -n $np -s {} -e {} ${res} >/dev/null 2>&1
  for t in $(find . -name "*.res.*"|awk -F. '{print $NF}'|sort|uniq); do
    [ $VERBOSE -ge 1 ] && echo "  "$res.$t:
    TAGT=$PWD/$res.$t
    REF=$cdir/$res.0.$t
    [ $VERBOSE -le 1 ] && OUT3=" >/dev/null 2>&1"
    sh -c "perl $compare_res $REF $TAGT $OUT3"
    COMPARE=$(($COMPARE+$?))
    [ $VERBOSE -ge 1 ] && if [ $COMPARE -eq 0 ]; then echo "    result: sucess";else echo "    result: failure"; ((ERRORS++)); fi
  done
  if [ $COMPARE -eq 0 ]; then echo "  result: sucess";else echo "  result: failure"; ((ERRORS++)); [ $goahead -eq 0 ] && die; fi
  [ $VERBOSE -ge 1 ] && echo "  time: ${SECONDS}s";
  rm -rf $tdir/*
  popd > /dev/null
done

rm -fr $test_dir
[ $ERRORS -ge 1 ] && die
exit 0

