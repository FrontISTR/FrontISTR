# Trilinos ML から MueLu への移行ガイド

## 概要

このドキュメントは、FrontISTRのTrilinos MLラッパーをMueLuに移行する際の手順と変更点について説明します。

## 主な変更点

### 1. ヘッダーファイルの変更

**ML (旧):**
```c
#include "ml_include.h"
#include "ml_config.h"
```

**MueLu (新):**
```cpp
#include "MueLu.hpp"
#include "MueLu_EpetraOperator.hpp"
#include "MueLu_CreateEpetraPreconditioner.hpp"
#include "MueLu_ParameterListInterpreter.hpp"
#include "Teuchos_ParameterList.hpp"
```

### 2. 言語の変更

- MLはCインターフェース
- MueLuはC++インターフェース（Cインターフェースも提供）
- ファイル拡張子を`.c`から`.cpp`に変更

### 3. API の変更

#### セットアップ

**ML (旧):**
```c
ML *ml_object;
ML_Create(&ml_object, N_grids);
ML_Init_Amatrix(ml_object, 0, nlocal, nlocal, id);
```

**MueLu (新):**
```cpp
Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::rcp(new Teuchos::ParameterList());
Teuchos::RCP<MueLu::EpetraOperator> muelu_prec = 
    MueLu::CreateEpetraPreconditioner(matrix, *paramList, nullspace);
```

#### パラメータ設定

**ML (旧):**
```c
ML_Gen_Smoother_Cheby(ml_object, level, ML_BOTH, 20.0, num_sweeps);
ML_Gen_Smoother_Amesos(ml_object, coarsest_level, ML_AMESOS_MUMPS, 1, 0.0);
```

**MueLu (新):**
```cpp
paramList->set("smoother: type", "CHEBYSHEV");
paramList->set("smoother: sweeps", num_sweeps);
paramList->set("coarse: type", "MUMPS");
```

### 4. 新機能と改善点

#### 4.1 改善されたパラメータ設定
- XML/YAML形式でのパラメータ設定サポート
- より直感的なパラメータ名
- 階層的なパラメータ構造

#### 4.2 新しいスムーザー
- ILUT スムーザー
- 改善されたChebyshevスムーザー
- ブロックスムーザーの強化

#### 4.3 新しい粗化スキーム
- HMIS (Heavy MIS) 粗化
- Coupled aggregation
- 改善されたパーティショニング

#### 4.4 パフォーマンス改善
- より効率的なメモリ使用
- 並列効率の向上
- 適応的パラメータ選択

## 移行手順

### ステップ 1: 依存関係の更新

#### CMakeLists.txt の変更
```cmake
# 旧
find_package(Trilinos REQUIRED COMPONENTS ML Amesos)

# 新
find_package(Trilinos REQUIRED COMPONENTS MueLu Amesos2 Tpetra)
```

#### コンパイルフラグの変更
```cmake
# 旧
target_compile_definitions(your_target PRIVATE HECMW_WITH_ML)

# 新
target_compile_definitions(your_target PRIVATE HECMW_WITH_MUELU)
```

### ステップ 2: ソースコードの移行

#### 2.1 ファイル名の変更
```bash
mv hecmw_ML_wrapper.c hecmw_MueLu_wrapper.cpp
```

#### 2.2 基本的な置換
```cpp
// 旧の関数名を新しい関数名に置換
hecmw_ML_wrapper_setup -> hecmw_MueLu_wrapper_setup
hecmw_ML_wrapper_apply -> hecmw_MueLu_wrapper_apply
hecmw_ML_wrapper_clear -> hecmw_MueLu_wrapper_clear
```

#### 2.3 データ構造の変更
```cpp
// 旧
struct ml_info {
    ML *ml_object;
    ML_Aggregate *agg_object;
    // ...
};

// 新
struct muelu_info {
    Teuchos::RCP<MueLu::EpetraOperator> muelu_prec;
    Teuchos::RCP<Epetra_CrsMatrix> matrix;
    // ...
};
```

### ステップ 3: パラメータマッピング

| ML パラメータ | MueLu パラメータ | 備考 |
|---------------|------------------|------|
| ML_AMESOS_MUMPS | "coarse: type": "MUMPS" | Amesos2使用 |
| ML_AMESOS_KLU | "coarse: type": "KLU2" | KLU2使用 |
| ML_MGV | "cycle type": "V" | V-cycle |
| ML_MGW | "cycle type": "W" | W-cycle |
| Cheby | "smoother: type": "CHEBYSHEV" | 改善版 |
| SymBlockGaussSeidel | "smoother: type": "RELAXATION" | + relaxation: type |

### ステップ 4: テスト

#### 4.1 基本動作テスト
```bash
# 既存のテストケースでの動作確認
cd tests/
./test.sh -t solver -s ml    # 旧ML
./test.sh -t solver -s muelu # 新MueLu
```

#### 4.2 性能比較
```bash
# ベンチマークテストの実行
cd benchmark/
./run_benchmark.sh --solver=ml
./run_benchmark.sh --solver=muelu
```

## 既知の問題と対策

### 1. メモリ使用量の増加
**問題:** MueLuはMLよりもメモリを多く使用する場合がある
**対策:** 
- `"repartition: enable": true` でメモリ分散
- `"aggregation: drop tol"` を調整

### 2. コンパイル時間の増加
**問題:** C++テンプレートによるコンパイル時間の増加
**対策:**
- 事前コンパイル済みヘッダーの使用
- 分割コンパイルの検討

### 3. デバッグ情報の変更
**問題:** MLとMueLuでデバッグ出力が異なる
**対策:**
- `"verbosity"` パラメータの調整
- ログレベルの統一

## 推奨設定

### 汎用設定（構造解析）
```cpp
paramList->set("verbosity", "low");
paramList->set("number of equations", 3);
paramList->set("max levels", 4);
paramList->set("coarse: type", "MUMPS");
paramList->set("smoother: type", "CHEBYSHEV");
paramList->set("smoother: sweeps", 2);
paramList->set("aggregation: type", "uncoupled");
```

### 大規模問題向け設定
```cpp
paramList->set("repartition: enable", true);
paramList->set("repartition: partitioner", "zoltan2");
paramList->set("repartition: min rows per proc", 800);
paramList->set("coarse: max size", 10000);
```

### メモリ制約がある場合
```cpp
paramList->set("aggregation: drop tol", 0.02);
paramList->set("coarse: type", "RELAXATION");
paramList->set("smoother: type", "RELAXATION");
```

## パフォーマンス最適化

### 1. スムーザーの選択
- **Chebyshev**: 一般的に最も効率的
- **Jacobi**: メモリ使用量が少ない
- **ILUT**: 収束性が良いが計算コストが高い

### 2. 粗化スキームの選択
- **Uncoupled**: 最も安定
- **HMIS**: 強結合問題に効果的
- **METIS/ParMETIS**: 不規則メッシュに適している

### 3. 並列化の最適化
- プロセス数に応じた `max levels` の調整
- `repartition` の有効活用
- `coarse: max size` の調整

## 参考資料

- [MueLu User's Guide](https://trilinos.github.io/muelu.html)
- [Trilinos Migration Guide](https://trilinos.github.io/migration.html)
- [MueLu Parameter Reference](https://trilinos.github.io/muelu_parameters.html)
