# WikiPathways Enrichment Analysis & Treemap Visualization

WikiPathways のエンリッチメント解析（Over-Representation Analysis）と階層的 Treemap 可視化を行う Jupyter Notebook です。

## 機能

1. **WikiPathways GMT ダウンロード** — 公式 GMT（Entrez ID）と Enrichr GMT（遺伝子シンボル）を自動取得
2. **エンリッチメント解析 (ORA)** — Fisher's exact test + Benjamini-Hochberg FDR 補正
3. **階層的 Treemap** — Plotly `go.Treemap` でカテゴリ別にパスウェイを可視化（-log₁₀ FDR 色付け）
4. **パスウェイ遺伝子マッピング** — WikiPathways の SVG ダイアグラム上にヒット遺伝子を赤色でハイライト

## セットアップ

```bash
pip install gseapy plotly matplotlib pandas requests jupyter
```

## 使い方

### 方法 1: 既存ノートブックを実行

```bash
jupyter notebook wikipathways_enrichment.ipynb
```

### 方法 2: ノートブックを再生成して実行

```bash
python create_notebook.py
jupyter notebook wikipathways_enrichment.ipynb
```

## カスタマイズ

- **遺伝子リスト**: セクション 3 の `gene_list` を自分のデータに差し替え
- **生物種**: セクション 2 の `SPECIES` と `ENRICHR_LIBRARY` を変更
  - マウス: `SPECIES = "Mus_musculus"`, `ENRICHR_LIBRARY = "WikiPathways_2024_Mouse"`
- **パスウェイ可視化**: セクション 5 の `pathway_index` で表示するパスウェイを選択

## 依存ライブラリ

- [GSEApy](https://github.com/zqfang/GSEApy) — エンリッチメント解析
- [Plotly](https://plotly.com/python/) — Treemap 可視化
- [Pandas](https://pandas.pydata.org/) — データ操作
- [Requests](https://requests.readthedocs.io/) — HTTP ダウンロード

## ライセンス

MIT
