# FrontISTR のコンパイル

FrontISTRをソースからコンパイルする手順を説明します。

## ライセンス

FrontISTRはMITライセンスで提供されています。詳しくは`License.txt`を参照してください。

## 準備

FrontISTRは、コンパイルに`cmake`を使います。`cmake`はバージョン 2.8.11 以上のものが必要です。

~~~txt
% cmake --version
cmake version 2.8.12.2
~~~

## コンパイル

最小限のFrontISTRは、以下の手順でコンパイルすることができます。

~~~txt
% tar xvf FrontISTR.tar.gz
% cd FrontISTR
% mkdir build
% cd build
% cmake ..
% make
% make test
% sudo make install
~~~

`/usr/local/bin`以下にコンパイルされた`fistr1`がインストールされます。

## `cmake`のオプション

他の機能を有効にしたい場合、`cmake`に以下のオプションを付けてください。ライブラリが既にインストールされていれば、自動的にその場所を探します。

| オプション        | 説明                                          | 備考                      |
|-------------------|-----------------------------------------------|---------------------------|
| -DWITH_TOOLS=1    | パーティショナなどのツールもコンパイルします  | hecmw_part1など           |
| -DWITH_MPI=1      | MPIを有効にします                             | ライブラリが必要          |
| -DWITH_OPENMP=1   | OpenMPを有効にします                          | コンパイラの対応が必要    |
| -DWITH_REFINER=1  | REVOCAP_Refinerの機能を有効にします           | ライブラリが必要          |
| -DWITH_REVOCAP=1  | REVOCAP_Couplerの機能を有効にします           | ライブラリが必要          |
| -DWITH_PARAC1=1   | 並列接触解析の機能を有効にします              | (未実装)                  |
| -DWITH_METIS=1    | METISの機能を有効にします                     | 4.0.3と5.1.0に対応        |
| -DMETIS_VER_4=1   | metis-4.0.3を使う場合に設定                   | metis-5.1.0の場合指定不要 |
| -DWITH_PARMETIS=1 | ParMETISの機能を有効にします                  | (未実装)                  |
| -DWITH_MKL=1      | MKL PARDISOの機能を有効にします               | (未実装)                  |
| -DWITH_MUMPS=1    | MUMPSの機能を有効にします                     | ライブラリが必要          |
| -DWITH_LAPACK=1   | LAPACKの機能を有効にします                    | ライブラリが必要          |
| -DWITH_ML=1       | Trilinos MLの機能を有効にします               | ライブラリが必要          |
| -DBUILD_DOC=1     | FrontISTRのソースコードをドキュメント化します | doxygenとgraphvizが必要   |

## LinuxなどのUnix系プラットフォームでの構築例

### 外部ライブラリの場所の指定

例えば`$HOME/tools`以下に、`include/`,`lib/`のようなディレクトリ構造で外部ライブラリをインストール
している場合

~~~txt
% CMAKE_INSTALL_PATH=$HOME/tools cmake ..
~~~

`$HOME/local`や`$HOME/.local`にインストールしている場合、指定は不要

### パーティショナ、MPI、OpenMPを有効にする

~~~txt
% cmake -DWITH_TOOLS=1 -DWITH_METIS=1 -DWITH_MPI=1 -DWITH_OPENMP=1 ..
~~~

### MUMPSを有効にする

~~~txt
% cmake -DWITH_MUMPS=1 -DWITH_MPI=1 ..
~~~

### MLを有効にする

~~~txt
% cmake -DWITH_ML=1 -DWITH_MPI=1 ..
~~~

### LaPACK, BLASを指定する

~~~txt
% cmake -DWITH_LAPACK=1 -DBLAS_LIBRARIES=$HOME/tools/lib/libopenblas.a -DLAPACK_LIBRARIES=$HOME/tools/lib/libopenblas.a ..
~~~txt

### Intel コンパイラを使う

~~~txt
% cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort ..
~~~

### MKLをLaPACK, BLASとして利用する

~~~txt
% source /opt/intel/mkl/bin/mklvar.sh intel64
% source /opt/intel/bin/compilervars.sh intel64
% echo $LD_LIBRARY_PATH (ライブラリが見えていることを確認)
% export BLA_VENDOR="Intel"
% cmake -DWITH_LAPACK=1 ..
~~~

### ソースコードのドキュメントを構築し参照する

FrontISTRのソースコードのドキュメンテーションを参照することができます。ただし、予め `doxygen`と`graphviz` をインストールしておく必要があります。

~~~txt
% cmake -DBUILD_DOC=1 ..
% make doc
% firefox doc/html/index.html
~~~

### 全てのオプションを有効にする

~~~txt
% cmake -DWITH_TOOLS=1 -DWITH_MPI=1 -DWITH_OPENMP=1 \
        -DWITH_REFINER=1 -DWITH_REVOCAP=1 -DWITH_METIS=1 \
        -DWITH_MUMPS=1 -DWITH_LAPACK=1 -DWITH_ML=1 -DBUILD_DOC=1 ..
~~~

### インストールディレクトリを変更する

`make install`でインストールする場所は`-DCMAKE_INSTALL_PREFIX`で指定することができます。デフォルトは`/usr/local`です。

~~~txt
% cmake -DCMAKE_INSTALL_PREFIX=$HOME/local ..
~~~

### デバッグを有効にする

デフォルトは`RELEASE`で構築します。デバッグを有効にしたい場合、指定してください。

~~~txt
% cmake -DCMAKE_BUILD_TYPE="DEBUG" ..
~~~

## Windowsでの構築例

FrontISTRは、MSYS2によるMinGW-w64でコンパイルが出来ることを確認しています。

MPIにMicrosoft MPIを使った場合の例を示します。

### Microsoft MPIのインストール

msmpisdk.msiとMSMpiSetup.exeを予めインストールしておきます。

~~~txt
% mkdir -p $HOME/local/lib
% mkdir $HOME/local/include
% cd $HOME/local/lib
% gendef /c/Windows/System32/msmpi.dll
% dlltool -d msmpi.def -l libmsmpi.a -D /c/Windows/System32/msmpi.dll
% cd $HOME/local/include
% cp "/c/Program Files (x86)/Microsoft SDKs/MPI/Include/"*.h .
% cp "/c/Program Files (x86)/Microsoft SDKs/MPI/Include/x64/"*.h .
~~~

コピーした `mpi.h` と `mpif.h` に変更を加えます。

~~~txt
% vi mpi.h
#ifndef MPI_INCLUDED
#define MPI_INCLUDED
のすぐ下に
#include <stdint.h>
を追加

% vi mpif.h
PARAMETER (MPI_ADDRESS_KIND=INT_PTR_KIND())
を
PARAMETER (MPI_ADDRESS_KIND=8)
に修正
~~~

変更が済んだらコンパイルしてください。

~~~txt
% cmake -G "MSYS Makefiles" -DWITH_MPI=1 \
        -DMPI_C_INCLUDE_PATH=$HOME/local/include \
        -DMPI_C_LIBRARIES=$HOME/local/lib/libmsmpi.a \
        -DMPI_CXX_INCLUDE_PATH=$HOME/local/include \
        -DMPI_CXX_LIBRARIES=$HOME/local/lib/libmsmpi.a \
        ..
% make
~~~

または

~~~txt
% cmake -G "MinGW Makefiles" -DWITH_MPI=1 \
        -DMPI_C_INCLUDE_PATH=$HOME/local/include \
        -DMPI_C_LIBRARIES=$HOME/local/lib/libmsmpi.a \
        -DMPI_CXX_INCLUDE_PATH=$HOME/local/include \
        -DMPI_CXX_LIBRARIES=$HOME/local/lib/libmsmpi.a \
        ..
% mingw32-make
~~~

でコンパイルしてください。

## OS X上でコンパイルする場合

(未検証)

## ライブラリをインストールする場合のTips

上記に示した機能を有効にする場合、予めライブラリをインストールしておく必要があります。

ライブラリやヘッダファイルは

1. システム上のパス(`/usr/lib`, `/usr/include`, `/usr/local/lib`, `/usr/local/include` など)
2. ホームディレクトリ上のパス(`$HOME/local/lib`, `$HOME/.local/lib` など)
3. 環境変数で指定したディレクトリを起点とするパス(`$METIS_ROOT`, `$MUMPS_ROOT` など)
4. FrontISTRのソースを展開した場所と同レベルのパス

などを探索します。見つからない場合や、意図しないパスの場合、手動で指定する必要があります。

例えば、LAPACKにOpenBLASを利用して、`$HOME/local`にインストールした場合、以下のような指定が必要になります。

~~~txt
% cmake -DWITH_LAPACK=1 \
        -DBLAS_LIBRARIES=$HOME/local/lib/libopenblas.a \
        -DLAPACK_LIBRARIES=$HOME/local/lib/libopenblas.a ..
~~~

`ccmake`や`cmake-gui`を使うと、これらを簡単に設定することができます。


## 外部ライブラリ (参考)

オプションを指定する際、予めライブラリをインストールしておく必要があるものがあります。

外部ライブラリは、`$HOME/local`へインストールするか、`FrontISTR`と同じディレクトリに置いておくと自動的に検索することができます。

### REVOCAP_Refiner

メッシュの細分化機能を有効にするには、REVOCAP_Refinerを構築する必要があります。

~~~txt
% tar xvf REVOCAP_Refiner-1.1.04.tar.gz
% cd REVOCAP_Refiner-1.1.04
% make
~~~

### REVOCAP_Coupler

FrontFlowなどと連成を行うためのカップラーです。

~~~txt
% tar xvf REVOCAP_Coupler-2.1.tar.gz
% cd REVOCAP_Coupler-2.1
% env REFINER_INCLUDE="-I$HOME/work/REVOCAP_Refiner-1.1.04/Refiner" \
     REFINER_LIBS="-L$HOME/work/REVOCAP_Refiner-1.1.04/lib/x86_64-linux -lRcapRefiner -lstdc++" \
     ./configure
~~~

### OpenBLAS

OpenBLASにはLAPACKも含まれているので、LAPACKを指定する部分をこれに置き換えることもできます。

~~~txt
% tar xvf OpenBLAS-0.2.18.tar.gz
% cd OpenBLAS-0.2.18
% make BINARY=64 NO_SHARED=1 USE_OPENMP=1
% make PREFIX=$HOME/local install
~~~

### Metis-5.1.0

Metis-5.1.0の`CMakeList.txt`には不具合があるため、トップディレクトリにある`CMakeLists.txt`を修正してからライブラリを構築してください。

~~~txt
set(GKLIB_PATH "GKlib" CACHE PATH "path to GKlib")
~~~

の部分を

~~~txt
set(GKLIB_PATH "${CMAKE_SOURCE_DIR}/GKlib" CACHE PATH "path to GKlib")
~~~

と修正してください。構築は一般的な `cmake` を使ったソフトウェアと同じです。

~~~txt
% cd build
% cmake -DCMAKE_INSTALL_PREFIX=$HOME/local -DOPENMP=1 ..
% make
% make install
~~~

### MUMPS-5.0.1

直接法ソルバーMUMPSの構築方法です。

`Make.inc`ディレクトリにある、`Makefile.inc.generic`を元に`Makefile.inc`を環境に合わせた内容へ書き換えます。

~~~txt
% cp Make.inc/Makefile.inc.generic
% vi Makefile.inc
LMETISDIR = $(HOME)/local
IMETIS    = -I$(LMETISDIR)/include
LMETIS    = -L$(LMETISDIR)/lib -lmetis

ORDERINGSF  = -Dmetis -Dpord

CC      = gcc -fopenmp
FC      = gfortran -fopenmp
FL      = gfortran -fopenmp

SCALAP  = -L$(HOME)/local/lib -lscalapack

INCPAR = -I/usr/lib/openmpi/include

LIBPAR  = $(SCALAP) -L/usr/lib/openmpi -lmpi -lmpi_mpifh

OPTF    = -O -DMUMPS_OPENMP
OPTC    = -O -I. -DMUMPS_OPENMP
OPTL    = -O
~~~

ファイルの修正が済んだらコンパイルをします。

~~~txt
% make
~~~

### scalapack-2.0.2

MPIを有効にしたMUMPSを構築するには、scalapack を事前に構築しておく必要があります。
また、現在の blacs は scalapack に統合されているため、別途 blacs をインストールする必要はありません。

~~~txt
% cd scalapack-2.0.2
% mkdir build
% cmake ..
% make
% cp libscalapack.a $HOME/local/lib
~~~

MinGWの場合、cmakeでビルド出来ないので`SLmake.inc`を編集する旧来の方法を用いる方が良いようです。

### Trilinos ML

FrontISTR で ML による疎行列のプリコンディショナーを利用する場合、以下の方法で Trilinos ML を構築してください。

~~~txt
% cmake -DTrilinos_ENABLE_ALL_OPTI1AL_PACKAGES=OFF -DTrilinos_ENABLE_ML=1 \
        -DTrilinos_ENABLE_Zoltan=1 -DTrilinos_ENABLE_OpenMP=1 -DTPL_ENABLE_MPI=1 ..
~~~

ただし、シリアル版を作成する場合は `-DTPL_ENABLE_MPI=1` の指定は不要です。

### パッケージの作成 (初期サポート)

構築したバイナリを、各プラットフォームでサポートするパッケージにすることできます。

TBZ2, DEB, RPM,

~~~txt
% cd build
% cpack
~~~
