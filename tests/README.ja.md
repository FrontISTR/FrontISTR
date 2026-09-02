[English Version](./README.md)

テストの実行方法
-----------------

[ctest][ctest]を使います。
FrontISTRプロジェクト自体を[cmake][cmake]で次の様にビルドしたと想定します

```
git clone https://gitlab.com/FrontISTR-Commons/FrontISTR
cd FrontISTR   # このディレクトリの事を以降 FRONTISTR_HOME と呼びます
cmake . -Bbuild
cmake --build build/ -j $(nproc)
```

これで `${FRONTISTR_HOME}/build` 以下に `fistr1` 等の実行バイナリが生成されています。
この状態で

```
cd build/
ctest
```

これで全てのテストが実行されます。

### テストのラベル

テストはラベルによって管理されていて、並列化について以下のラベルが存在します

| label | OpenMP | MPI |
|:------|:------:|:---:|
|serial | OFF    | OFF |
|openmp | ON     | OFF |
|mpi    | OFF    | ON  |
|hybrid | ON     | ON  |

これらのラベルが付いたテストだけを実行するには例えば

```
ctest -L mpi
```

の様に `-L` (`--label-regex`) を付けて実行します。

テストには並列化に関するものの他に以下に説明するように `target` によるものがあります。
テストは下記の節で説明されているように [cmake][cmake] によって自動的に追加されますが、
例えば `${FRONTISTR_HOME}/tests/analysis/eigen/exK/` ディレクトリにあるテストには
`${FRONTISTR_HOME}/tests/` からの相対パスをそのまま文字列として用いて `analysis/eigen/exK` というラベルを付けます
このディレクトリにあるテストだけを実行するには次の様にします

```
ctest -L analysis/eigen/exK
```

また `-L` は部分マッチで選択出きるので、

```
ctest -L analysis
```

とすれば `${FRONTISTR_HOME}/tests/analysis` 以下の全てのテストを実行します。
　
### 出力の調整

`ctest` の一般的な使い方として、

```
ctest -V
```

は全ての出力を表示、

```
ctest --output-on-failure
```

は失敗したテストの結果のみ出力を表示します。詳しくは `ctest -h` を確認してください。
　
テストの追加方法
-----------------

テストは [cmake][cmake] が `${FRONTISTR_HOME}/tests/` ディレクトリ以下にあるメッシュデータ(`*.msh`)を自動的に走査して、
発見したものを自動的に登録します。
このテストは信頼できるリファレンスと現在のソースコードに基づく `fistr1` での計算結果を比較して、
それが十分に小さいかどうかを判定します。

したがって、テストを新たに追加するには

1. `${FRONTISTR_HOME}/tests` 以下に新たにディレクトリを作成する
    - analysis, lib, solver配下に作成されたディレクトリは常にテスト実行対象となります
    - with_[mkl|mumps|ml]配下に作成されたディレクトリはcmakeオプション-DWITH_[MKL|MUMPS|ML]がONのときにテスト実行対象となります
    - _archive配下に作成されたディレクトリはテスト対象に含まれません
2. `*.msh` ファイル、`*.cnt` ファイルを追加する
3. `${FRONTISTR_HOME}/tests/create_reference.sh` でリファレンスデータを生成する
4. 計算結果が妥当かどうかを手動で確認する

というプロセスを経ます。

`create_reference.sh` の実行には`${FRONTISTR_HOME}/build/fistr/fistr1`をデフォルトで利用しますので事前にビルドが必要です。
`create_reference_docker.sh` の実行にはFrontISTRの公式リリースのコンテナイメージを用いるので
[Docker][docker] の実行権限が必要です。

### 周波数応答解析のテスト

周波数応答解析（`!SOLUTION,TYPE=DYNAMIC` で `!DYNAMIC` の1行目第2項が `2`）は、
先に実行した固有値解析の固有モードを重ね合わせて応答を求めます。
`<メッシュ名>.msh` と同じディレクトリに `<メッシュ名>_eigen.cnt` を置くと、
`test.sh` と `create_reference.sh` は周波数応答解析の前に固有値解析を実行し、
その出力を次のように受け渡します。

| 受け渡すもの | ファイル名 |
|:-------------|:-----------|
| 固有値（ログ）       | `eigen_0.log` |
| 固有ベクトル（結果） | `<メッシュ名>_eigen.res` |
| 時刻歴の出力先       | `<メッシュ名>_dyna.res` |

したがって `<メッシュ名>.cnt` の `!EIGENREAD` にはログファイルとして `eigen_0.log` を指定します。
リファレンスとの比較は、周波数掃引の結果 `<メッシュ名>.res.0.*` と
時刻歴の結果 `<メッシュ名>_dyna.res.0.*` の両方について行われます。

`${FRONTISTR_HOME}/tests/analysis/freq` 配下のテストは、変化させる条件ごとにディレクトリを分けています。
テストケースを追加するときは、確認したい条件のディレクトリに追加し、
それ以外の条件は下記の基準構成のままにすることで、変化させた条件の影響だけが現れるようにします。

| ディレクトリ | 変化させる条件 | ケース |
|:-------------|:---------------|:-------|
| `element`  | 要素タイプ | `FQ341` `FQ342` `FQ351` `FQ352` `FQ361` `FQ362` |
| `load`     | `!FLOAD`   | `FQL01` 実部、`FQL02` 虚部、`FQL03` 複数群・複数自由度の複素荷重、`FQL04` 面群 |
| `boundary` | `!BOUNDARY`| `FQB01` 片持ち、`FQB02` すべり支持を追加、`FQB03` 両端拘束 |
| `modal`    | 減衰とモード範囲 | `FQM01` 質量比例減衰、`FQM02` 質量+剛性比例減衰、`FQM03` `!EIGENREAD` を2次モードから |

基準構成は、361要素の 1.0 x 0.2 x 0.1 の片持ちはりを `FIX` で固定し、
`LOADP` にz方向の荷重を与え、5モード、レイリー減衰 alpha=0・beta=1.0E-4、
250Hzまでを5点で掃引するものです。


[cmake]: https://cmake.org/cmake/help/latest/manual/cmake.1.html
[ctest]: https://cmake.org/cmake/help/latest/manual/ctest.1.html
[docker]: https://www.docker.com/
