# SX-Aurora TSUBASA 対応について
SXAT2ブランチは、NEC SX-Aurora TSUBASAでの実行用に最適化を行なったバージョンの
ソースコードが格納されています。

本バージョンは、vhcallと呼ばれる仕組みを用いて、剛性行列の作成部分をCPUで
それ以降の計算をベクトルエンジン上で行なう実行方法を採用しています。

vhcallではベクトルエンジンでの実行を担当する通常のプログラムから、CPUで実行される処理をまとめた共有ライブラリを呼び出すことでこのようなヘテロジニアスな実行形態を実現しています。
このため、CPUでの実行部分を担当する共有ライブラリと、ベクトルエンジンで実行されこの共有ライブラリを呼び出す、実行ファイルの2つをそれぞれビルドする必要があります。

また、並列実行時に使用する領域分割ツール `hecmw_part1` などベクトルエンジン上での実行に不向きなツールもありますので、これらのツール類はCPU用にビルドしたものを用いることを推奨します。


## ビルド方法
ソースコードのトップレベルにビルド用のディレクトリを作成し、そのディレクトリ内からcmakeを実行するものとします。

1. CPU上で動作するユーティリティツール、共有ライブラリのビルド(intelコンパイラを使用する場合)

```
CC=icc CXX=icpc FC=ifort cmake \
  -DBUILD_SHARED=YES\
  -DWITH_DOC=OFF\
  -DWITH_METIS=YES\
  -DWITH_MKL=OFF\
  -DWITH_LAPACK=OFF\
  ..
make
make install
```

2. ベクトルエンジン上で動作する実行ファイルのビルド

```
cmake \
  -DUSE_VHCALL=YES\
  -DUSE_HETERO_SOLVER=ON\
  -DWITH_DOC=OFF\
  -DWITH_MKL=OFF\
  -DWITH_METIS=OFF\
  -DWITH_LAPACK=OFF\
  -DWITH_TOOLS=OFF\
  -DCMAKE_TOOLCHAIN_FILE=../../cmake/AuroraToolchain.cmake\
  ..
make
make install
```

ビルドに用いるディレクトリはCPU用とベクトルエンジン用でそれぞれ別のディレクトリを作成するか
一方のビルドが終わった後でディレクトリ内を全て削除してから次のビルドに進んでください。

全てのビルドが正常に終了すると、インストール先のディレクトリには次のファイルが格納されています。

- fistr1        ---- ベクトルエンジン用の実行ファイル
- hec2rcap      ---- CPU用の実行ファイル
- hecmw\_part1  ---- CPU用の実行ファイル
- hecmw\_vis1   ---- CPU用の実行ファイル
- libfistr.so   ---- fistr1から呼び出す共有ライブラリファイル
- rconv         ---- CPU用の実行ファイル
- rmerge        ---- CPU用の実行ファイル


## 実行方法
実行時に環境変数 `LD_LIBRARY_PATH` に、libfistr.soを保存したパスを追加してください。

プログラムの実行は、通常のSX-Aurora TSUBASA上で動かすプログラムと同様
mpiexecを使用してください。

SX-Aurora TSUBASA版 FrontISTRでは線形ソルバとして、NEC製直接法ソルバライブラリ `HeteroSolver` を使用することができます。
本ソルバを使用する時は、解析制御データの `METHOD`  行に `HeteroSolver` を指定してください。

その他のプログラムの実行方法は通常のFrontISTRに準じます。
