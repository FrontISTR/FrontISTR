# FrontISTR のコンパイル

FrontISTRをソースからコンパイルする手順を説明します。

## ライセンス

FrontISTRはMITライセンスで提供されています。詳しくは`License.txt`を参照してください。

## 準備

FrontISTRは、コンパイルに`cmake`を使います。`cmake`はバージョン 2.8.11よりも新しいものが必要です。

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

| オプション         | 説明                                            | 備考                      |
|--------------------|-------------------------------------------------|---------------------------|
| -DWITH_TOOLS=ON    | パーティショナなどのツールもコンパイルします。  | hecmw_part1など           |
| -DWITH_MPI=ON      | MPIを有効にします。                             | ライブラリが必要          |
| -DWITH_OPENMP=ON   | OpenMPを有効にします。                          | コンパイラの対応が必要    |
| -DWITH_REFINER=ON  | REVOCAP_Refinerの機能を有効にします。           | ライブラリが必要          |
| -DWITH_REVOCAP=ON  | REVOCAP_Couplerの機能を有効にします。           | ライブラリが必要          |
| -DWITH_PARACON=ON  | 並列接触解析の機能を有効にします。              | (未実装)                  |
| -DWITH_METIS=ON    | METISの機能を有効にします。                     | 4.0.3と5.1.0に対応        |
| -DMETIS_VER_4=ON   | metis-4.0.3を使う場合に設定                     | metis-5.1.0の場合指定不要 |
| -DWITH_PARMETIS=ON | ParMETISの機能を有効にします。                  | (未実装)                  |
| -DWITH_MKL=ON      | MKL PARDISOの機能を有効にします。               | (未実装)                  |
| -DWITH_MUMPS=ON    | MUMPSの機能を有効にします。                     | ライブラリが必要          |
| -DWITH_LAPACK=ON   | LAPACKの機能を有効にします。                    | ライブラリが必要          |
| -DWITH_ML=ON       | Trilinos MLの機能を有効にします。               | ライブラリが必要          |
| -DBUILD_DOC=ON     | FrontISTRのソースコードをドキュメント化します。 | doxygenとgraphvizが必要   |

## LinuxなどのUnix系プラットフォームでの構築例

### パーティショナ、MPI、OpenMPを有効にする

~~~txt
% cmake -DWITH_TOOLS=ON -DWITH_MPI=ON -DWITH_OPENMP=ON ..
~~~

### MUMPSを有効にする

~~~txt
% cmake -DWITH_MUMPS=ON ..
~~~

### ソースコードのドキュメントを構築し参照する

FrontISTRのソースコードのドキュメンテーションを参照することができます。ただし、予め `doxygen`と`graphviz` をインストールしておく必要があります。

~~~txt
% cmake -DBUILD_DOC=ON ..
% make doc
% firefox doc/html/index.html
~~~

### 全てのオプションを有効にする

~~~txt
% cmake -DWITH_TOOLS=ON -DWITH_MPI=ON -DWITH_OPENMP=ON \
        -DWITH_REFINER=ON -DWITH_REVOCAP=ON -DWITH_METIS=ON \
        -DWITH_MUMPS=ON -DWITH_LAPACK=ON -DWITH_ML=ON -DBUILD_DOC=ON ..
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

## Windows上でコンパイルする場合

MSYS2による、MinGW-w64でコンパイルが出来ることを確認しています。

~~~txt
% cmake -G "MSYS Makefiles" ..
% make
~~~

または

~~~txt
% cmake -G "MinGW Makefiles" ..
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
3. 環境変数で指定したディレクトリを起点とするパス (`$METIS_ROOT`, `$MUMPS_ROOT` など)
4. FrontISTRのソースを展開した場所と同レベルのパス

などを探索します。

見つからない場合や、意図しないパスの場合、手動で指定する必要があります。

例えば、LAPACKにOpenBLASを利用して、`$HOME/local`にインストールした場合、以下のような指定が必要になります。

~~~txt
% cmake -DWITH_LAPACK=ON \
        -DBLAS_LIBRARIES=$HOME/local/lib/libopenblas.a \
        -DLAPACK_LIBRARIES=$HOME/local/lib/libopenblas.a ..
~~~

`ccmake`や`cmake-gui`を使うと、これらを簡単に設定することができます。


## 外部ライブラリ (参考)

オプションを指定する際、予めライブラリをインストールしておく必要があるものがあります。

外部ライブラリは、`$HOME/local`へインストールするか、`FrontISTR`と同じディレクトリに置いておくと自動的に検索することができます。

### OpenBLAS

OpenBLASにはLaPACKも含まれているので、LAPACKを指定する部分をこれに置き換えることもできます。

~~~txt
% tar xvf OpenBLAS-0.2.18.tar.gz
% cd OpenBLAS-0.2.18
% make BINARY=64 NO_SHARED=1 USE_OPENMP=1
% make PREFIX=/path/to/install install
~~~

### Trilinos ML

FrontISTR で ML による疎行列のプリコンディショナーを利用する場合、以下の方法で Trilinos ML を構築してください。

~~~txt
% cmake -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF -DTrilinos_ENABLE_ML=ON \
        -DTrilinos_ENABLE_Zoltan=ON -DTrilinos_ENABLE_OpenMP=ON -DTPL_ENABLE_MPI=ON ..
~~~

ただし、シリアル版を作成する場合は `-DTPL_ENABLE_MPI=ON` の指定は不要。

### Metis-5.1.0

Metis-5.1.0の CMakeList.txt には不具合があるため、CMakeLists.txtを修正してからライブラリを構築する必要がある。 トップディレクトリにある CMakeLists.txt を修正する。

~~~txt
set(GKLIB_PATH "GKlib" CACHE PATH "path to GKlib")
~~~

となっているところを

~~~txt
set(GKLIB_PATH "${CMAKE_SOURCE_DIR}/GKlib" CACHE PATH "path to GKlib")
~~~

と修正してください。構築は一般的な `cmake` を使ったソフトウェアと同じです。

~~~txt
% cd build
% cmake -DOPENMP=ON ..
% make
% make install
~~~

### scalapack-2.0.2

MPIを有効にしたMUMPSを構築するには、scalapackを事前に構築しておく必要があります。
また、現在blacsとscalapackは統合されているため、blacsを別途インストールする必要はありません。

~~~txt
% cd scalapack-2.0.2
% mkdir build
% cmake ..
% make
% make install
~~~

MinGWの場合、cmakeでビルド出来ないので`SLmake.inc`を編集する旧来の方法を用いる方が良いようです。

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
% env REFINER_INCLUDE="-I../REVOCAP_Refiner-1.1.04/Refiner" \
     REFINER_LIBS="-L../REVOCAP_Refiner-1.1.04/lib/x86_64-linux -lRcapRefiner -lstdc++" \
     ./configure
~~~

### パッケージの作成 (初期サポート)

構築したバイナリを、各プラットフォームでサポートするパッケージにすることできます。

TBZ2, DEB, RPM, 

~~~txt
% cd build
% cpack
~~~


