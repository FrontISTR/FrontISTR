# 修正プログラム一覧

本改修の主な目的は、MITC4シェル要素（要素タイプ `741`）を幾何学的非線形解析（`NLGEOM`）で扱えるようにすることである。
対象は、主に **4節点MITC4シェル要素、Total Lagrangian、弾性材料** の組み合わせである。

今回の修正では、シェル節点の回転を単なる回転自由度として扱うのではなく、節点ごとの姿勢情報として管理するようにした。具体的には、参照姿勢、現在姿勢、Newton反復中の試行姿勢を保持し、回転増分によって節点triadおよびdirectorを更新する。また、directorまわりのdrilling成分は物理的なdirectorの傾きとは分離し、別の状態量として扱う。

| No. | ファイル | 主な修正内容 |
| --- | --- | --- |
| 1 | `fistr1/src/lib/static_LIB_shell.f90` | MITC4シェル要素に対して、有限回転およびTotal Lagrangianに基づく剛性行列、内部力、応力評価の処理を追加した。 |
| 2 | `fistr1/src/analysis/static/fstr_Update.f90` | シェル節点の参照姿勢、現在姿勢、Newton反復中の試行姿勢を管理する処理を追加した。また、回転増分によるtriad更新と、収束時の状態確定処理を追加した。 |
| 3 | `fistr1/src/analysis/static/fstr_CreateMatrix.f90` | 接線剛性行列を計算する際に、シェルの参照姿勢と現在の試行姿勢を要素計算へ渡す処理を追加した。 |
| 4 | `fistr1/src/analysis/static/fstr_Residual.f90` | 残差評価時に、変位だけでなくシェル節点の回転状態も試行状態として反映する処理を追加した。 |
| 5 | `fistr1/src/analysis/static/fstr_AddBC.f90` | 回転境界条件に対して、指定された全体回転と現在回転の差から相対的な増分回転を作成する処理を追加した。 |
| 6 | `fistr1/src/analysis/static/fstr_Cutback.f90` | cutback発生時に、シェル節点の回転状態および厚み方向Gauss点の履歴を保存・復元する処理を追加した。 |
| 7 | `fistr1/src/analysis/static/fstr_NodalStress.f90` | シェル厚み方向の応力・ひずみ履歴を用いて、有限回転後のdirectorと整合する節点応力出力処理を追加した。 |
| 8 | `fistr1/src/analysis/static/fstr_solve_NLGEOM.f90` | `NLGEOM`ステップ開始時に、シェル節点の回転状態を保存し、Newton反復で使用する試行状態を初期化する処理を追加した。 |
| 9 | `fistr1/src/analysis/static/fstr_solve_NonLinear.f90` | Newton反復内で、線形ソルバから得られた回転増分を用いてシェル節点の試行姿勢を更新し、収束時にその状態を確定する処理を追加した。 |
| 10 | `fistr1/src/analysis/static/fstr_solve_QuasiNewton.f90` | Quasi-Newton解析でも、通常のNewton解析と同様にシェル節点の回転増分更新および収束時の状態確定を行うように修正した。 |
| 11 | `fistr1/src/common/fstr_setup.f90` | シェル節点の回転状態配列、および厚み方向Gauss点履歴を確保・初期化・解放する処理を追加した。 |
| 12 | `fistr1/src/lib/m_fstr.F90` | `fstr_solid` に、シェル節点の参照姿勢、現在姿勢、試行姿勢、drilling成分を保持するための配列を追加した。 |
| 13 | `fistr1/src/lib/physics/mechgauss.f90` | `tElement` に、シェル厚み方向Gauss点の履歴情報を保持する配列と、その保存・復元に関する補助処理を追加した。 |
