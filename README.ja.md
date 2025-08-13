# FrontISTR

オープンソースの**大規模並列FEM（有限要素法）構造解析ソフトウェア**です。非線形解析・接触・大変形・動解析・固有値・周波数応答・熱伝導など、産業／学術の実務で必要な機能を備えています。ライセンスは **MIT License** です。

* 公式サイト: [https://www.frontistr.com](https://www.frontistr.com)
* 公式ドキュメント: [https://manual.frontistr.com/en/](https://manual.frontistr.com/en/) （英語） / [https://manual.frontistr.com/ja/](https://manual.frontistr.com/ja/) （日本語）
* ソースコード: **GitLab（開発の本流）** [https://gitlab.com/FrontISTR-Commons/FrontISTR](https://gitlab.com/FrontISTR-Commons/FrontISTR) ／ **GitHub（ミラー）** [https://github.com/FrontISTR/FrontISTR](https://github.com/FrontISTR/FrontISTR)

---

## 主な特徴

* **大規模並列**: MPI（分散メモリ）＋ OpenMP（共有メモリ）
* **解析種別**

  * 線形静解析
  * 非線形静解析（弾塑性・ハイパーエラスト・クリープ・粘弾性 等）
  * 動解析（線形／非線形過渡）
  * 固有値解析
  * 周波数応答解析
  * 熱伝導解析（定常・非定常）
  * 接触解析（動的接触を含む）
* **ソルバ**: CG / BiCGSTAB / GMRES / GPBiCG / 直接法（MUMPS / MKL PARDISO）
* **前処理**: SSOR / ILU(0) / AMG（Trilinos-ML） 等
* **要素ライブラリ**: 1次/2次のソリッド（TET/PRISM/HEX）、平面（三角/四角）、ビーム、シェル、トラス、インタフェース要素 等
* **可視化**: 結果から `hecmw_vis` 等でVTK可視化データを生成し、ParaView等で表示可能

> **注**: メッシュ生成・可視化ツールは同梱されません（REVOCAP\_PrePost、FreeCAD FEM\_FrontISTR、ParaView 等の外部ツールを利用）。

---

## 動作環境・依存関係

* **OS**: Linux（推奨）、Windows（MSYS2 / WSL2 を含む）
* **コンパイラ**: C/C++/Fortran90（GCC / Clang / Intel oneAPI 等）
* **必須/推奨ライブラリ**

  * **MPI**（OpenMPI / MPICH）
  * **METIS**（領域分割）※未導入でもRCBで分割可能
  * **HEC-MW**（同梱、自動ビルド）
  * **BLAS/LAPACK**（一部機能で利用。OpenBLAS/MKL 等）
  * **MUMPS**（並列直接法）＋ **ScaLAPACK**（MUMPSに必要）
  * **Trilinos-ML**（AMG前処理）
  * **Intel MKL**（一部の接触解析機能で必須／MKL PARDISOを利用）

> 依存関係の詳細と推奨バージョンはインストールマニュアルを参照してください。

---

## ビルド（CMake）

最小構成の例：

```bash
# 依存ライブラリを事前に用意
# 例: OpenMPI/METIS/OpenBLAS など

git clone https://gitlab.com/FrontISTR-Commons/FrontISTR.git
cd FrontISTR
mkdir build && cd build
cmake ..
make -j$(nproc)
make install   # CMAKE_INSTALL_PREFIX でインストール先指定可
```

主なCMakeオプション（抜粋）:

* `-DWITH_MPI=ON` / `-DWITH_OPENMP=ON`
* `-DWITH_TOOLS=ON`（`hecmw_part1` などのツールを同時ビルド）
* `-DWITH_METIS=ON`（`-DMETIS_VER_4=ON` でmetis 4系を指定可能）
* `-DWITH_MUMPS=ON` / `-DWITH_LAPACK=ON` / `-DWITH_MKL=ON` / `-DWITH_ML=ON`
* `-DWITH_REFINER=ON`（REVOCAP\_Refiner）/ `-DWITH_REVOCAP=ON`（REVOCAP\_Coupler）
* `-DOLD_RES_FORMAT=OFF`（旧結果フォーマットへ戻す場合はON）

> `./fistr1 -v` でビルド時に有効化された機能を確認できます。

---

## クイックスタート

**単一領域（スレッド並列）**

```bash
# 入力ファイルが置かれたディレクトリで
fistr1 -t 4    # 4スレッドで実行
```

**領域分割（MPI並列）**

```bash
# 例題: tutorial/02_elastic_hinge_parallel を想定
hecmw_part1           # メッシュ分割（hecmw_part_ctrl.dat に従う）
mpirun -n 4 fistr1    # 4並列で実行（環境に応じて mpiexec など）
```

---

## 入出力ファイル

最低限必要な入力：

* `hecmw_ctrl.dat`（全体制御：メッシュ/制御ファイルの指定 等）
* `<case>.cnt`（解析制御：境界条件・荷重・ソルバ設定 等）
* `<case>.msh`（メッシュ・材料/セクション定義）

代表的な出力：

* `*.res.*`（解析結果：MPIランクごとに出力）
* `*_vis_*`（可視化用データ。VTK等に変換して可視化）
* `0.log`（ログ）

---

## 可視化・プリポスト

* **ParaView**（VTKの読み込み）
* **REVOCAP\_PrePost**（メッシュ生成・前処理／実行連携）
* **FreeCAD FEM\_FrontISTR**（FrontISTR用Workbench）: [https://github.com/FrontISTR/FEM\_FrontISTR](https://github.com/FrontISTR/FEM_FrontISTR)

---

## バイナリ／配布

* **Linux**: 主要ディストリ向けにソースからのビルドが基本です。
* **Windows**: MSYS2 で `pacman -S mingw-w64-x86_64-frontistr` によりインストール可能（環境により追加設定が必要）。
* **HPC**: 一部スパコン環境ではモジュールとして提供される場合があります（運用ガイドに従ってください）。

---

## プロジェクト構成・貢献

* イシュー／マージリクエストは **GitLab** 側で受け付けています。
* 変更提案の際は、

  1. Fork → ブランチ作成
  2. コーディング規約・テスト方針に従い修正
  3. MR（Merge Request）を作成（再現手順・ベンチマーク・既知の制約を記載）

`CONTRIBUTING` と各種テンプレート（課題・MR）を参照してください。

---

## ライセンス

* 本体：**MIT License**（リポジトリ内 `LICENSE.txt` を参照）
* 連携ライブラリ（METIS/MUMPS/Trilinos-ML/MKL 等）は各ライセンスに従います。

---

## 謝辞

FrontISTRは、大学・研究機関・企業からなる開発コミュニティ（FrontISTR Commons）によって継続的に開発されています。関連する研究開発プロジェクト（革新的シミュレーションソフトウェア等）の支援に感謝します。

---

## 参考リンク（抜粋）

* インストールガイド（CMake / 依存関係 / テスト）: [https://manual.frontistr.com/en/install/](https://manual.frontistr.com/en/install/)
* 解析フロー（実行方法・入出力）: [https://manual.frontistr.com/en/analysis/analysis\_01.html](https://manual.frontistr.com/en/analysis/analysis_01.html)
* 要素ライブラリ一覧: [https://manual.frontistr.com/en/analysis/analysis\_02.html](https://manual.frontistr.com/en/analysis/analysis_02.html)
* チュートリアル一覧: [https://manual.frontistr.com/en/tutorial/](https://manual.frontistr.com/en/tutorial/)
* Doxygen（API/内部構造）: [https://frontistr-commons.gitlab.io/FrontISTR/doxygen/index.html](https://frontistr-commons.gitlab.io/FrontISTR/doxygen/index.html)
