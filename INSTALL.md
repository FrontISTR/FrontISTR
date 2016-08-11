# FrontISTRのインストール {#mainpage}

FrontISTRの一般的なインストール方法を説明します。<br>
オプションを指定する事で、様々な機能を有効にする事ができます。 詳細はマニュアルを参照して下さい。

## Quick Start

FrontISTRは、以下の手順でインストールすることが出来ます。

~~~txt
$ mkdir build
$ cd build
$ cmake ..
$ make
$ make test
$ sudo make install
~~~

## cmake のオプション

オプションを指定してビルドすることで、様々な機能が有効になります。

- ツールをビルドする `-DWITH_TOOLS=ON`
- MPI を有効にする `-DWITH_MPI=ON`
- Metis を有効にする `-DWITH_METIS=ON`

  - Metis のライブラリを指定する `-DMETIS_LIBRARY=/path/libmetis.a`
  - Metis のインクルードパスを指定する `-DMETIS_INCLUDE_DIR=/path/include`
  - Metis-4系列を利用する `-DMETIS_VER_4=ON` (Metis-5系列の場合、指定する必要なし)

- OpenMP を有効にする `-DWITH_OPENMP=ON`

- REVOCAP_Refiner を有効にする `-DWITH_REFINER=ON`

- REVOCAP_Coupler を有効にする `-DWITH_COUPLER=ON`

- LaPACK を有効にする `-DWITH_LAPACK=ON`

- ML を有効にする `-DWITH_ML=ON`

- Doxygen によるドキュメントの生成 `-DBUILD_DOCUMENTATION=ON`

  - ドキュメントを生成するには `$ make doc`
  - 生成されたドキュメントは `firefox doc/html/index.html` としてブラウザで開くことが出来る。

- インストール先の指定 `-DCMAKE_INSTALL_PREFIX=$HOME/local`

  - METISやMLなどの外部ライブラリのインストール先を `${CMAKE_INSTALL_PREFIX}` にすると、自動でライブラリとインクルードファイルの場所を認識する( `${CMAKE_INSTALL_PREFIX}/lib` や `${CMAKE_INSTALL_PREFIX}/include` にファイルが有る場合)。

## 構築例

代表的な構築例を示します。

### Linux (Parallel)

MPI実行可能なソルバを構築します。

~~~txt
$ cmake -DWITH_MPI=ON ..
~~~

### Linux (Parallel with Metis)

ソルバに加えて、メッシュファイルを分割するためのパーティショナも構築します。

~~~txt
$ cmake -DWITH_TOOLS=ON -DWITH_MPI=ON -DWITH_METIS -DMETIS_INCLUDE_DIR=$HOME/metis-4.0.3/Lib -DMETIS_LIBRARY=$HOME/metis-4.0.3/libmetis.a ..
~~~

### Linux (Parallel with OpenMP Metis LaPACK and ML)

高性能な反復法ソルバのプリコンディショナーMLを有効にします。

~~~txt
$ cmake -DWITH_TOOLS=ON -DWITH_MPI=ON -DWITH_OPENMP=ON -DWITH_METIS=ON -DWITH_LAPACK=ON -DWITH_ML=ON ..
~~~

### MSYS

Windowsで実行可能なバイナリを作成します。

~~~txt
$ cmake -G "MinGW Makefiles" ..
$ mingw32-make
~~~

### Mac

## 外部ライブラリ (参考)

オプションを指定する際、予めライブラリをインストールしておく必要があるものがある。

### OpenBLAS

OpenBLASにはLaPACKも含まれているので、LaPACKを指定する部分をこれに置き換えることもできます。

~~~txt
$ tar xvf OpenBLAS-0.2.18.tar.gz
$ cd OpenBLAS-0.2.18
$ make BINARY=64 NO_SHARED=1 USE_OPENMP=1
$ make PREFIX=/path/to/install install
~~~

### Trilinos ML

FrontISTR で ML による疎行列のプリコンディショナーを利用する場合、以下の方法で Trilinos ML を構築する必要がある。

~~~txt
$ cmake -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF -DTrilinos_ENABLE_ML=ON -DTrilinos_ENABLE_Zoltan=ON -DTrilinos_ENABLE_OpenMP=ON -DTPL_ENABLE_MPI=ON ..
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

と修正する。構築は一般的な `cmake` を使ったソフトウェアと同じ。

~~~txt
$ cd build
$ cmake -DOPENMP=ON ..
$ make
$ make install
~~~

### scalapack-2.0.2

MPIを有効にしたMUMPSを構築するには、scalapackを事前に構築しておく必要があります。<br>
また、現在blacsとscalapackは統合されているため、blacsを別途インストールする必要はありません。

~~~txt
$ cd scalapack-2.0.2
$ mkdir build
$ cmake ..
$ make
$ make install
~~~

MinGWの場合、cmakeでビルド出来ないので`SLmake.inc`を編集する旧来の方法を用いる方が良い。

### MUMPS-5.0.1

### REVOCAP_Refiner

メッシュの細分化機能を有効にするには、REVOCAP_Refinerを構築する必要があります。

~~~txt
$ tar xvf REVOCAP_Refiner-1.1.04.tar.gz
$ cd REVOCAP_Refiner-1.1.04
$ make
~~~

### REVOCAP_Coupler

FrontFlowなどと連成を行うためのカップラー。MinGWでは少々修正しないといけない。

~~~txt
$ tar xvf REVOCAP_Coupler-2.1.tar.gz
$ cd REVOCAP_Coupler-2.1
$ env REFINER_INCLUDE="-I../REVOCAP_Refiner-1.1.04/Refiner" \
     REFINER_LIBS="-L../REVOCAP_Refiner-1.1.04/lib/x86_64-linux -lRcapRefiner -lstdc++" \
     ./configure
~~~
