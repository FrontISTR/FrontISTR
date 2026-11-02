# FrontISTR のレビュー基準

Merge Request（MR）のレビューで確認することと、その進め方を定めます。

## 基本スタンス

- レビューは、価値のある変更をどうすれば本体に取り込めるか、提出者と一緒に考えるために行う。
- 修正を求める際は、理由が明らかな場合を除き、その理由を示す。
- マージ前に必要な対応と、後から改善できる課題を区別する。必要な品質を確保しつつ、そのMRに不要な作業を求めない。
- 見落としや残課題を理由に、レビュワー個人の責任を問わない。取り込み後に見つかった問題や、そのMRで残した課題は、新しいIssueとして扱う。
- レビュワーだけでは判断できない機能の採否や、取り込み後の機能の削除は、Software WGで最終判断する。

## レビューで確認すること

- 以下の観点でMRを確認し、対応が必要な問題があれば修正を求める。
- 確認した範囲と判断できなかった点を明らかにし、必要に応じて他のレビュワーと確認を分担する。
- 要求する修正の範囲は、変更の影響に見合ったものとする。修正要求の水準は、実際のレビュー例を参考にする。
- 修正要求では具体的な問題と影響を示し、今回の変更による問題と、もともと存在する問題を区別する。

| 観点 | 確認すること | 実際のMR・レビュー例 |
|---|---|---|
| 正しさ | 意図した動作や数理と整合しているか。関連する解析条件、並列実行、異常時の処理でも問題がないか。 | [同じ物理モデルで結果が変わる点を指摘](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/700#note_4080960348) ／ [内力・応力への座標系の渡し忘れを指摘](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/699#note_4068021673) ／ [一部ランクが通信を飛ばす問題を指摘](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/576#note_3824847816) |
| 検証 | 変更に見合う検証が行われているか。追加したテストが対象の機能を実行し、不具合や意図しない変更を検出できるか。 | [修正を戻しても通るテストを指摘](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/706#note_4080839899) ／ [追加テストがCIで実行されていない点を指摘](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/689#note_4081518373) ／ [圧縮機能を通らない圧縮テストを指摘](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/674#note_3942568647) |
| 設計・責任分担 | 各処理の責任分担や階層が保たれ、他の層が担うべき詳細が入り込んでいないか。同種の処理で使われる既存パターンと整合しているか。 | [時間増分の管理から接触固有の処理を分離](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/691#note_4041722311) ／ [アセンブルと要素計算を分離](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/593#note_3828414349) ／ [ミドルウェアへの解析条件の持込みを指摘](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/642#note_3885639989) |
| 入力・出力仕様 | 既存の項目名・書式・意味の使い分けに沿っているか。互換性や利用者への影響に問題がないか。 | [解析入力での機種指定を指摘](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/597#note_3751145248) ／ [減衰係数をオプションからデータ行へ移動](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/577) |
| 複雑さ・拡張性 | 過剰な分岐・重複・依存や、処理のまとまりの不明確さによって、理解や今後の変更を難しくしていないか。他の機能に不要な負担を増やしていないか。 | [剛性行列作成の二重管理を指摘](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/504#note_2845420845) ／ [要素ごとに条件分岐を散らさず集約するよう指摘](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/643#note_3827680479) ／ [材料処理の集約と責務の明確化（模範例）](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/706#note_4080839899) |
| 説明と実装の一致 | コメント、ログ、MR本文、マニュアルが実装と合っているか。対応範囲や既存入力への影響が説明されているか。 | [修正後も古いままの説明・計算結果を指摘](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/700#note_4094847272) ／ [ログとマニュアルの不一致を指摘](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/674#note_3942568647) |
| 知財 | コード・コメント・テストデータ等について、他ソフトの資材の不適切な利用など、第三者の著作権・特許権に関わる懸念がないか。疑いに気づいた場合に確認するものとし、網羅的な調査は求めない。 | [CSparse由来の逐語訳とLAPACKの出典表示漏れを指摘](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/674#note_3942568647) ／ [ライセンス表示の不足を再確認](https://gitlab.com/FrontISTR-Commons/FrontISTR/-/merge_requests/674#note_3946065466) |

### 設計・責任分担の補足

- 設計の基準は、既存コードやコメントから読み取れる構造・意図と、同種の処理で使われているパターンとする。設計書にまとまっていることは前提としない。
- 明確な既存パターンや記載がない詳細は、実装者・レビュワーの判断に委ねる。
- 新たな設計方針を決めた場合は、以後の実装・レビューでも読み取れるよう、コードの構成やコメントに残す。

## レビューの主な対象にしないこと

- **細かなコード作法やインデント等**  
  指摘や修正要求を妨げるものではないが、機械チェック等で統一することを目指す。
- **機能自体の有用性や優先度**  
  上記の基準を満たす機能については、レビュワー個人の賛同をマージ条件にしない。必要性に疑問がある場合は、気になった人が別途Issueを立てて議論することを原則とする。ただし、明らかにFrontISTRの対象範囲から外れるものは、レビュワーの判断で受入れを断ってよい。

## Software WGに判断を委ねる変更

レビュワーだけでは取り込みを判断できない場合や、以下のように共有する設計・仕様・互換性について新たな方針判断が必要な変更は、論点を整理してSoftware WGに判断を委ねる。既存の方針やパターンに沿う修正・拡張は、通常のレビューで進める。

- 共通部分の責任分担や依存関係について、既存の方針を変える変更
- 他の機能の追加・修正に大きな制約を生み、その制約を許容するか判断が必要な変更
- 既存入力・出力の意味や互換性を変える変更、既存パターンに沿わない新しいI/O形式の導入

WGに判断を委ねる場合も、上記のレビュー項目に沿った確認・修正要求を進め、論点や選択肢を具体化することは歓迎する。
