# valgrind によるメモリ検査

## 概要

`fistr1` を valgrind の memcheck で実行すると、FrontISTR 自身のメモリ確保に加えて、
リンクしているライブラリの確保も報告される。特に MPI 実行では、Open MPI が
MPI_Init で登録した MCA コンポーネントや、hwloc が構築したトポロジ情報が
プロセス終了まで保持されるため、小さなモデルでも1プロセスあたり数百件の記録が出る。

このディレクトリには、そうした FrontISTR の管理外のライブラリに由来する記録を抑止する
valgrind の suppression ファイルを、対象ごとに置いてある。指定すると、
残る記録は FrontISTR 自身が確保したブロックだけになる。

| ファイル | 抑止する対象 |
|:---------|:-------------|
| `openmpi.supp` | Open MPI 本体（libmpi、libmpi_mpifh、libopen-pal、libopen-rte、libpmix）と、それが引き込む libevent・hwloc とそのプラグイン・libpciaccess・libudev・libOpenCL、およびプラグイン読み込み時の `dlopen` |
| `libgomp.supp` | GCC の OpenMP ランタイム |

ファイルを分けてあるのは、対象が入れ替わりうるためである。
Open MPI 以外の MPI 実装ではライブラリ名が異なるので `openmpi.supp` は効かず、
`libgomp.supp` だけを使うことになる。同様に、Intel コンパイラでビルドした場合の
OpenMP ランタイムは libiomp5 なので `libgomp.supp` は効かない。

各エントリは、確保時のスタックがそのライブラリの中を通るブロックだけに一致する。
FrontISTR 自身が確保したブロックは、確保時のスタックが `fistr1` の中にとどまるため、
どのエントリにも一致しない。

## 実行例

以下では、FrontISTR のソースツリーを `$FRONTISTR_HOME`、ビルドディレクトリを
`$FRONTISTR_HOME/build` とし、解析の入力ファイル一式がカレントディレクトリに
置かれているものとする。`--suppressions` は必要な数だけ並べられる。

### 逐次実行

MPI を有効にしてビルドした `fistr1` は、逐次実行でも起動時に MPI_Init を呼ぶため、
Open MPI の分も抑止する必要がある。

```sh
valgrind --leak-check=full --show-leak-kinds=all \
  --suppressions=$FRONTISTR_HOME/doc/valgrind/openmpi.supp \
  --suppressions=$FRONTISTR_HOME/doc/valgrind/libgomp.supp \
  $FRONTISTR_HOME/build/fistr1/fistr1
```

### MPI 実行

プロセスごとに別のログファイルへ出力する。`%q{OMPI_COMM_WORLD_RANK}` は
Open MPI が各プロセスに与える環境変数を展開するもので、ランク番号が
ファイル名に入る。他の MPI 実装ではプロセスIDを使う `%p` を指定する。

```sh
mpirun -n 2 valgrind --leak-check=full --show-leak-kinds=all \
  --suppressions=$FRONTISTR_HOME/doc/valgrind/openmpi.supp \
  --suppressions=$FRONTISTR_HOME/doc/valgrind/libgomp.supp \
  --log-file=valgrind.%q{OMPI_COMM_WORLD_RANK}.log \
  $FRONTISTR_HOME/build/fistr1/fistr1
```

valgrind は実行を数十倍遅くするため、検査には小さなモデルを使い、
`OMP_NUM_THREADS` も小さくしておくとよい。

残った記録の呼び出し元をさらに遡りたい場合は `--num-callers` を指定する
（既定は12フレーム）。

## 出力の読み方

各ログの末尾に集計が出る。`suppressed` が抑止された分である。

```
==10== LEAK SUMMARY:
==10==    definitely lost: 0 bytes in 0 blocks
==10==    indirectly lost: 0 bytes in 0 blocks
==10==      possibly lost: 0 bytes in 0 blocks
==10==    still reachable: 26,922 bytes in 36 blocks
==10==         suppressed: 89,678 bytes in 837 blocks
```

- `definitely lost` と `indirectly lost` は、そのブロックを指すポインタが
  どこにも残っていないもので、解放漏れとして直す対象である
- `still reachable` は、終了時点でポインタが残っているものである。
  プロセス終了で回収されるため実害は無いが、解析の途中で確保した作業領域が
  ここに現れる場合は、解放忘れの可能性がある
- `possibly lost` は、ブロックの先頭以外を指すポインタしか残っていないもので、
  Fortran の配列記述子などで生じることがある

## 抑止されない報告が出た場合

ライブラリの構成が異なる環境では、用意した suppression が想定していない経路から
報告が出ることがある。その場合は `--gen-suppressions=all` を付けて実行すると、
各報告の直後に、それを抑止するための suppression が出力される。

```sh
valgrind --leak-check=full --show-leak-kinds=all --gen-suppressions=all \
  $FRONTISTR_HOME/build/fistr1/fistr1
```

出力された suppression のうち、FrontISTR の管理外のライブラリに由来するものだけを、
対象に応じたファイルへ既存の書き方に合わせて追加する。すなわち、確保元のライブラリを
`obj:*/<ライブラリ名>.so*` の1フレームで指定し、その前後を `...` で挟む。
関数名やアドレスをそのまま並べると、ライブラリのバージョンが変わるたびに
一致しなくなるため避ける。
別の MPI 実装や別の OpenMP ランタイムを対象にする場合は、
`openmpi.supp` に倣って新しいファイルを追加する。
