# Contributing to FrontISTR

Thank you for your interest in contributing to FrontISTR.
Development takes place on GitLab: https://gitlab.com/FrontISTR-Commons/FrontISTR

## Reporting a bug / requesting a feature

- Open an Issue on GitLab.
- Attaching a minimal input set (mesh `*.msh` and control `*.cnt`) that reproduces the issue helps us a lot.

## Contributing a change (bug fix / new feature)

1. Fork the repository on GitLab and create a branch.
2. Make your change. If it adds or changes behavior, please update the documentation as well: https://gitlab.com/FrontISTR-Commons/frontistr_documents
3. Build and run the tests:
   ```sh
   mkdir build && cd build
   cmake .. && make -j$(nproc)
   ctest
   ```
   The same tests run in CI, so please make sure they pass.
4. Adding a small model (mesh / cnt) under `tests/` helps keep your change covered: a case that exercises a new feature, or the input that reproduced a fixed bug kept as a regression test. Encouraged, not required.
5. Open a Merge Request against `master`.

## Where new code goes

Each file holds roughly one module, and each directory has a clear role: analysis-type-specific drivers live under `fistr1/src/analysis/`, reusable analysis parts (elements, material models, contact, ...) live in `lib/`, and control-file I/O and setup live in `common/`. Following the existing layout keeps a change easy to review and to build on.

If you are not sure where something belongs, that is fine — open the Merge Request anyway and we will help you place it. Finding the right boundary is part of review, not something you need to settle first.

Contributors are listed in `CONTRIBUTORS.txt` — feel free to add your name.

---

# FrontISTR への貢献

FrontISTR への貢献に興味を持っていただきありがとうございます。
開発は GitLab で行っています: https://gitlab.com/FrontISTR-Commons/FrontISTR

## バグ報告・機能要望

- GitLab で Issue を立ててください。
- 問題を再現する最小の入力（メッシュ `*.msh`・制御 `*.cnt`）を添えていただけると助かります。

## 変更の提案（バグ修正・新機能）

1. GitLab でリポジトリを fork し、ブランチを作成する。
2. 変更を行う。挙動の追加・変更を伴う場合はドキュメントも更新する: https://gitlab.com/FrontISTR-Commons/frontistr_documents
3. ビルドしてテストを実行する:
   ```sh
   mkdir build && cd build
   cmake .. && make -j$(nproc)
   ctest
   ```
   同じテストが CI でも実行されるので、パスすることを確認してください。
4. `tests/` に小さなモデル（mesh / cnt）を追加すると、変更がカバーされ続けます。新機能を動かすケースや、バグを再現した入力を回帰テストとして残すかたちです。必須ではありませんが、ぜひ。
5. `master` 向けに Merge Request を作成する。

## 新しいコードの置き場所

1ファイルはおおよそ1モジュール、ディレクトリごとに役割が分かれています。解析種別ごとのドライバは `fistr1/src/analysis/` に、解析のパーツ（要素・構成則・接触など）は `lib/` に、制御ファイルの入出力やセットアップは `common/` に置きます。既存の構成にならうと、レビューも以後の拡張もしやすくなります。

置き場所に迷っても大丈夫です。まず Merge Request を出していただければ、一緒に置き場所を考えます。適切な境界を見つけるのはレビューの一部で、最初に決めきる必要はありません。

貢献者は `CONTRIBUTORS.txt` に記載されます。お名前を追加していただいて構いません。
