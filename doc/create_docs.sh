#!/bin/bash

set -eux

export LC_ALL=C.UTF-8
export LANG=C.UTF-8

script_dir=$(cd $(dirname $0); pwd)
manual_dir=${script_dir}/manuals

rm -rf $manual_dir
mkdir -p $manual_dir

cd $script_dir/FrontISTR_manual/markdown_files/ja/
  mkdocs build --clean
  mv site $manual_dir/manual_ja
cd -

cd $script_dir/FrontISTR_manual/markdown_files/en/
  mkdocs build --clean
  mv site $manual_dir/manual_en
cd -
