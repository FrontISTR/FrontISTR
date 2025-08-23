# FrontISTR

![Logo](./doc/FrontISTR_logo.svg)

[![CI Status](https://gitlab.com/frontistr-commons/frontistr/badges/master/pipeline.svg)](https://gitlab.com/frontistr-commons/frontistr/-/pipelines)
[![Documentation](https://img.shields.io/badge/docs-latest-blue)](https://manual.frontistr.com/en/)
[![License](https://img.shields.io/badge/license-MIT-green)](License.txt)

FrontISTRはオープンソースの**大規模並列有限要素法構造解析ソフトウェア**です。材料非線形・大変形・接触を含む非線形の変形応力・固有値・周波数応答・熱伝導など、産業／学術の実務で必要な解析機能を備えています。ライセンスは **MIT License** です。

* 公式サイト: [https://www.frontistr.com](https://www.frontistr.com)
* 公式ドキュメント: [https://manual.frontistr.com/en/](https://manual.frontistr.com/en/) （英語） / [https://manual.frontistr.com/ja/](https://manual.frontistr.com/ja/) （日本語）

---

## 主な特徴

* **大規模並列**: MPI（分散メモリ）＋ OpenMP（共有メモリ）
* **解析種別**
  * 静解析、動解析（陰解法・陽解法）
  * 大変形、接触、非線形材料（弾塑性・超弾性・クリープ・粘弾性 等）
  * 固有値解析
  * 周波数応答解析
  * 熱伝導解析（定常・非定常）
* **ソルバ**: CG / BiCGSTAB / GMRES / GPBiCG / 直接法（MUMPS / MKL PARDISO）
* **前処理**: SSOR / ILU(0) / AMG（Trilinos-ML） 等
* **要素ライブラリ**: 1次/2次のソリッド（TET/PRISM/HEX）、平面（三角/四角）、はり、シェル、トラス要素 等
* **可視化**: ParaView等で表示可能

---

## 動作環境・依存関係

* **コンパイラ**: C/C++/Fortran90（GCC / Clang / Intel oneAPI 等）
* **並列計算時の必須ライブラリ**
  * **MPI**（Open MPI / MPICH 等）
  * **METIS**（領域分割）
* **オプションのライブラリ**
  * **BLAS/LAPACK**（一部機能で利用。OpenBLAS/MKL 等）
  * **MUMPS**（並列直接法）＋ **ScaLAPACK**（MUMPSに必要）
  * **Trilinos-ML**（AMG前処理）
  * **Intel MKL**（MKL PARDISOによる並列直接法）

---

## ビルド（CMake）

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

---

## 実行

**単一領域（スレッド並列）**

```bash
# 入力ファイルが置かれたディレクトリで
fistr1
```

**領域分割（MPI並列）**

```bash
hecmw_part1           # メッシュ分割
mpirun -n 4 fistr1    # 4並列で実行
```

---

## 可視化・プリポスト

* **ParaView**（VTKファイル）
* **REVOCAP\_PrePost**（メッシュ生成・前処理／実行連携）
* **FreeCAD FEM\_FrontISTR**（FrontISTR用Workbench）: [https://github.com/FrontISTR/FEM\_FrontISTR](https://github.com/FrontISTR/FEM_FrontISTR)

---

## バイナリ／配布

ダウンロードサイトはこちら  
https://www.frontistr.com/download/

---

## 貢献

* Issue / Merge Request(MR)は **[GitLab](https://gitlab.com/FrontISTR-Commons/FrontISTR)** 側で受け付けています。
* 変更提案の際は`CONTRIBUTING.md` と各種テンプレート（課題・MR）を参照してください。

---

## ライセンス

* 本体：**MIT License**（リポジトリ内 `License.txt` を参照）
* 連携ライブラリ（METIS/MUMPS/Trilinos-ML/MKL 等）は各ライセンスに従います。

---

## 謝辞

FrontISTRは、大学・研究機関・企業からなる開発コミュニティ（FrontISTR Commons）によって継続的に開発されています。関連する研究開発プロジェクト（革新的シミュレーションソフトウェア等）の支援に感謝します。  
GPUへの移植について東京大学情報基盤センターの支援に感謝します。
