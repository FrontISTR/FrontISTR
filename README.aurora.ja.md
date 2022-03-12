# SX-Aurora TSUBASA 対応について
SXAT2ブランチは、NEC SX-Aurora TSUBASAでの実行用に最適化を行なったバージョンの
ソースコードが格納されています。

本バージョンは、vhcallと呼ばれる仕組みを用いて、剛性行列の作成部分をCPUで
それ以降の計算をベクトルエンジン上で行なう実行方法を採用しています。

vhcallではベクトルエンジンでの実行を担当する通常のプログラムから、CPUで実行される処理をまとめた共有ライブラリを呼び出すことでこのようなヘテロジニアスな実行形態を実現しています。
このため、CPUでの実行部分を担当する共有ライブラリと、ベクトルエンジンで実行されこの共有ライブラリを呼び出す、実行ファイルの2つをそれぞれビルドする必要があります。

また、並列実行を行なう際に使用する領域分割ツール `hecmw_part1` はベクトルエンジン上での実行に不向きなプログラムなので、これらのユーティリティツールも別途CPU用にビルドして使用することを強く推奨します。


## ビルド方法
ソースコードのトップレベルにビルド用のディレクトリを作成し、そのディレクトリ内からcmakeを実行するものとします。

1. CPU上で動作するユーティリティツールのビルド(intelコンパイラを使用する場合)

```
CC=icc CXX=icpc FC=ifort cmake \
  -DTOOLS_ONLY=YES\
  -DCMAKE_Fortran_FLAGS="-fpp"\
  -DWITH_METIS=YES\
  -DWITH_MKL=OFF\
  -DWITH_LAPACK=OFF\
  -DCMAKE_EXE_LINKER_FLAGS=-lrt\
  ..
make
make install
```

2. CPU上での動作部分を担当する共有ライブラリのビルド(intelコンパイラを使用する場合)

```
CC=icc CXX=icpc FC=ifort cmake \
  -DCMAKE_Fortran_FLAGS="-fpp"\
  -DBUILD_SHARED=YES\
  -DWITH_DOC=OFF\
  -DWITH_METIS=OFF\
  -DWITH_MKL=OFF\
  -DWITH_LAPACK=OFF\
  -DWITH_TOOLS=OFF\
  -DCMAKE_INSTALL_PREFIX=${REMOTE_INSTALLDIR}\
  -DCMAKE_EXE_LINKER_FLAGS=-lrt\
  ..
make
make install
```

3. ベクトルエンジン上で動作する実行ファイルのビルド

```
cmake \
  -DUSE_VHCALL=YES\
  -DCMAKE_Fortran_FLAGS="-fpp"\
  -DUSE_HETERO_SOLVER=ON\
  -DWITH_DOC=OFF\
  -DWITH_MKL=OFF\
  -DWITH_METIS=OFF\
  -DWITH_LAPACK=OFF\
  -DWITH_TOOLS=OFF\
  -DCMAKE_INSTALL_PREFIX=${REMOTE_INSTALLDIR}\
  -DCMAKE_TOOLCHAIN_FILE=../../cmake/AuroraToolchain.cmake\
  -DCMAKE_EXE_LINKER_FLAGS="-lrt -lvhcall_fortran -lheterosolver_mpi_openmp -lsblas_sequential"\
  ..
make
make install
```

ビルドに用いるディレクトリは全て異なるディレクトリを用いるか、1段階ビルドが終わる毎に削除してから次の段階のビルドを行なってください。

全てのビルドが正常に終了すると、インストール先のディレクトリには次のファイルが格納されています。

- fistr1        ---- ベクトルエンジン用の実行ファイル
- hec2rcap      ---- CPU用の実行ファイル
- hecmw\_part1  ---- CPU用の実行ファイル
- hecmw\_vis1   ---- CPU用の実行ファイル
- libfistr.so   ---- fistr1から呼び出す共有ライブラリファイル
- -rconv        ---- CPU用の実行ファイル
- rmerge        ---- CPU用の実行ファイル


## 実行方法
実行ディレクトリに、前述のlibfistr.soをコピーするか、このファイルへのシンボリックリンクを作成してください。
なお、のファイルは、通常の共有ライブラリの探索手順とは異なり環境変数 `LD_LIBRARY_PATH` に指定されたディレクトリの探索を行ないません。必ず共有ライブラリファイルのコピーまたはシンボリックリンクをカレントディレクトリに置いた状態で実行してください。


SX-Aurora TSUBASA版 FrontISTRでは線形ソルバとして、NEC製直接法ソルバライブラリ `HeteroSolver` を使用することができます。
本ソルバを使用する時は、解析制御データの `METHOD`  行に `HeteroSolver` を指定してください。

その他のプログラムの実行方法は通常のFrontISTRに準じます。

