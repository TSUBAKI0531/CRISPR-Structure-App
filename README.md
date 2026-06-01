# CRISPR Structural Pipeline

**テーマ**: CAR-T細胞療法

## 概要

CAR-T細胞療法における非ウイルス的CARノックインで重要なHDR（相同組換え修復）とNHEJ（非相同末端結合）の修復効率を、CARコンストラクトのアミノ酸配列から予測するプロトタイプツールです。BioPythonで算出した物理化学的特性を用いたシミュレーションにより、HDR/NHEJ比率およびノックイン成功率スコアを算出します。結果はPass/Alert 2段階で判定し、設計の改善指針を提示します。

> **注意**: 現行実装はアミノ酸配列ベースのシミュレーションです。gRNA・核酸配列入力とHDR/NHEJ本格解析への対応は今後の拡張候補です。

## 入力

- アミノ酸配列（VH/VL または CDR3）
- 形式: FASTAまたはプレーンテキスト

## 出力

| スコア | 説明 |
|--------|------|
| Predicted HDR repair efficiency | HDR修復効率予測スコア |
| Predicted NHEJ repair efficiency | NHEJ修復効率予測スコア |
| Predicted Knock-in success rate | ノックイン成功率予測スコア |

判定: **Pass / Alert** 2段階（アラート閾値はサイドバーで調整可能）

## 使用技術

- Python
- Streamlit
- BioPython（ProteinAnalysis: MW / pI / Aromaticity / GRAVY 算出）

## ローカル起動方法

```bash
pip install -r requirements.txt
streamlit run app.py
```
