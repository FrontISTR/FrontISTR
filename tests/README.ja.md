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


[cmake]: https://cmake.org/cmake/help/latest/manual/cmake.1.html
[ctest]: https://cmake.org/cmake/help/latest/manual/ctest.1.html
[docker]: https://www.docker.com/
