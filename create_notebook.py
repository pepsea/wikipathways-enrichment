#!/usr/bin/env python3
"""Generate the WikiPathways enrichment analysis Jupyter Notebook."""

import json, os, uuid

cells = []

def md(source):
    cells.append({
        "cell_type": "markdown", "metadata": {},
        "id": uuid.uuid4().hex[:8],
        "source": [line + "\n" for line in source.split("\n")[:-1]] + [source.split("\n")[-1]]
    })

def code(source):
    lines = source.split("\n")
    cells.append({
        "cell_type": "code", "metadata": {}, "execution_count": None,
        "id": uuid.uuid4().hex[:8],
        "source": [line + "\n" for line in lines[:-1]] + [lines[-1]],
        "outputs": []
    })

# ===== Title =====
md("""# WikiPathways エンリッチメント解析 & Treemap可視化

このNotebookでは以下を行います：
1. WikiPathwaysのGMTファイルをローカルにダウンロード
2. Over-Representation Analysis (ORA) によるエンリッチメント解析
3. Treemapによる結果の可視化""")

# ===== Section 1 =====
md("""---
## 1. セットアップ & ライブラリのインポート""")

code("""# 必要なライブラリのインストール（未インストールの場合）
!pip install -q gseapy plotly matplotlib pandas requests""")

code("""import os
import requests
import re
import pandas as pd
import numpy as np
import gseapy as gp
import plotly.graph_objects as go
from collections import OrderedDict, defaultdict
import warnings
warnings.filterwarnings('ignore')

print(f"GSEApy version: {gp.__version__}")
print("セットアップ完了 ✓")""")

# ===== Section 2 =====
md("""---
## 2. WikiPathways GMTファイルのダウンロード

[WikiPathways](https://www.wikipathways.org/) から最新のGMTファイルをローカルにダウンロードします。

> **注意**: 公式GMTはEntrez Gene IDを使用しています。エンリッチメント解析には遺伝子シンボル対応のGMTも別途取得します。""")

code("""# --- 設定 ---
SPECIES = "Homo_sapiens"          # 対象の生物種（必要に応じて変更）
DATA_DIR = "data"
os.makedirs(DATA_DIR, exist_ok=True)

# ====================================
# 2-1. 公式WikiPathways GMTのダウンロード
# ====================================
GMT_BASE_URL = "https://data.wikipathways.org/current/gmt/"

def download_official_gmt(species: str, data_dir: str) -> str:
    \"\"\"WikiPathways公式サイトからGMTファイルをダウンロード（Entrez Gene ID版）\"\"\"
    print(f"WikiPathways GMTアーカイブを確認中...")
    index_resp = requests.get(GMT_BASE_URL, timeout=30)
    index_resp.raise_for_status()
    
    pattern = rf'(wikipathways-\d+-gmt-{species}\.gmt)'
    match = re.search(pattern, index_resp.text)
    if not match:
        raise FileNotFoundError(f"{species} のGMTファイルが見つかりません")
    
    filename = match.group(1)
    url = GMT_BASE_URL + filename
    filepath = os.path.join(data_dir, filename)
    
    if os.path.exists(filepath):
        print(f"既存ファイルを使用: {filepath}")
        return filepath
    
    print(f"ダウンロード中: {url}")
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    
    with open(filepath, 'w') as f:
        f.write(resp.text)
    
    print(f"✓ 保存完了: {filepath} ({len(resp.text):,} bytes)")
    return filepath

official_gmt = download_official_gmt(SPECIES, DATA_DIR)

# ====================================
# 2-2. Enrichr WikiPathways GMTのダウンロード（遺伝子シンボル版）
# ====================================
def download_enrichr_gmt(library_name: str, data_dir: str) -> str:
    \"\"\"Enrichrから遺伝子シンボル版GMTをダウンロード\"\"\"
    filepath = os.path.join(data_dir, f"{library_name}.gmt")
    
    if os.path.exists(filepath):
        print(f"既存ファイルを使用: {filepath}")
        return filepath
    
    print(f"Enrichrからダウンロード中: {library_name}")
    url = f"https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={library_name}"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    
    with open(filepath, 'w') as f:
        f.write(resp.text)
    
    print(f"✓ 保存完了: {filepath} ({len(resp.text):,} bytes)")
    return filepath

ENRICHR_LIBRARY = "WikiPathways_2024_Human"  # 遺伝子シンボル版
enrichr_gmt = download_enrichr_gmt(ENRICHR_LIBRARY, DATA_DIR)
print(f"\\n使用するGMTファイル:")
print(f"  公式 (Entrez ID): {official_gmt}")
print(f"  解析用 (Symbol):  {enrichr_gmt}")""")

code("""# GMTファイルの内容を確認
def parse_gmt(filepath: str) -> dict:
    \"\"\"GMTファイルをパースしてパスウェイ名 -> 遺伝子セットの辞書を返す\"\"\"
    pathways = OrderedDict()
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\\t')
            if len(parts) >= 3:
                pathway_name = parts[0]
                genes = [g for g in parts[2:] if g]  # 空文字除去
                if genes:
                    pathways[pathway_name] = genes
    return pathways

pathways = parse_gmt(enrichr_gmt)
all_genes = set(g for genes in pathways.values() for g in genes)

print(f"=== WikiPathways GMTファイル概要 ({ENRICHR_LIBRARY}) ===")
print(f"パスウェイ数: {len(pathways)}")
print(f"遺伝子数（ユニーク）: {len(all_genes):,}")
print()
print("--- パスウェイ例（先頭5件）---")
for i, (name, genes) in enumerate(pathways.items()):
    if i >= 5:
        break
    print(f"  {name}: {len(genes)} genes")
    print(f"    例: {', '.join(genes[:5])}...")""")

# ===== Section 3 =====
md("""---
## 3. エンリッチメント解析（Over-Representation Analysis）

### サンプル遺伝子リスト

以下はデモ用のサンプル遺伝子リストです。**実際のデータに差し替えて使用してください。**""")

code("""# =====================================================
# ★ ここにご自身の遺伝子リストを入力してください ★
# =====================================================
# 以下はデモ用のサンプルリスト（免疫・炎症関連遺伝子）

gene_list = [
    # 炎症性サイトカイン
    "TNF", "IL1B", "IL6", "IL8", "IFNG", "IL10", "IL17A",
    # NF-κBシグナル
    "NFKB1", "NFKB2", "RELA", "RELB", "IKBKB", "CHUK",
    # JAK-STATシグナル
    "JAK1", "JAK2", "JAK3", "STAT1", "STAT3", "STAT4", "STAT5A", "STAT6",
    # TLRシグナル
    "TLR2", "TLR4", "MYD88", "IRAK1", "IRAK4", "TRAF6",
    # アポトーシス関連
    "CASP3", "CASP8", "CASP9", "BCL2", "BAX", "BID",
    # MAPKシグナル
    "MAPK1", "MAPK3", "MAPK8", "MAPK14", "MAP2K1", "MAP3K7",
    # PI3K-Aktシグナル
    "PIK3CA", "PIK3CB", "AKT1", "MTOR", "PTEN",
    # Wntシグナル
    "CTNNB1", "GSK3B", "APC", "WNT1", "FZD1",
    # TP53関連
    "TP53", "MDM2", "CDKN1A", "RB1",
    # 転写因子
    "JUN", "FOS", "MYC", "HIF1A",
]

print(f"入力遺伝子数: {len(gene_list)}")
print(f"遺伝子リスト: {', '.join(gene_list[:10])}...")""")

md("""### ORA（Over-Representation Analysis）の実行

Fisher's exact testを用いてパスウェイの濃縮度を評価します。""")

code("""# GSEApyを使ったOver-Representation Analysis
enr = gp.enrich(
    gene_list=gene_list,          # 解析対象の遺伝子リスト
    gene_sets=enrichr_gmt,        # WikiPathways GMTファイル（遺伝子シンボル版）
    outdir=None,                  # 結果をファイルに保存しない
    no_plot=True,
    cutoff=0.05,                  # FDR閾値
)

# 結果をDataFrameとして取得
results_df = enr.results.copy()

# カラム名を日本語化して表示用に整形
display_df = results_df[[
    'Term', 'Overlap', 'P-value', 'Adjusted P-value',
    'Odds Ratio', 'Combined Score', 'Genes'
]].copy()
display_df.columns = [
    'パスウェイ', 'ヒット数/パスウェイサイズ', 'P値', 'FDR補正P値',
    'オッズ比', 'Combined Score', '遺伝子'
]

# FDR補正P値でソート
display_df = display_df.sort_values('FDR補正P値')

print(f"=== エンリッチメント解析結果 ===")
print(f"有意なパスウェイ数 (FDR < 0.05): {len(display_df)}")
print()

# 上位20件を表示
display_df.head(20).style.format({
    'P値': '{:.2e}',
    'FDR補正P値': '{:.2e}',
    'オッズ比': '{:.2f}',
    'Combined Score': '{:.2f}'
}).background_gradient(subset=['FDR補正P値'], cmap='RdYlGn')""")

# ===== Section 4 =====
md("""---
## 4. 階層的 Treemap可視化 (Plotly)

エンリッチメント解析結果を **Plotly go.Treemap** で可視化します。
- **階層**: Root → カテゴリ → パスウェイ（キーワードベースで自動分類）
- **面積**: ヒット遺伝子数に比例
- **色**: -log₁₀(FDR adjusted P-value)（子孫含む最大値）""")

code("""# === パスウェイをカテゴリに分類する関数 ===
def classify_pathway(name):
    \"\"\"パスウェイ名からカテゴリを推定する\"\"\"
    name_lower = name.lower()
    
    categories = [
        ('Cancer / Tumor', [
            'cancer', 'tumor', 'carcinoma', 'leukemia', 'lymphoma',
            'melanoma', 'glioblastoma', 'oncog', 'metasta',
        ]),
        ('Immune / Inflammation', [
            'immune', 'inflam', 'interleukin', 'il-', 'il1', 'il2', 'il4',
            'il6', 'il8', 'il10', 'il12', 'il13', 'il17', 'il18', 'il26',
            'interferon', 'toll', 'tlr', 'nf-kb', 'nfkb', 'nf-κb',
            'complement', 'cytokine', 'chemokine', 'phagocyt',
            't cell', 'b cell', 'mhc', 'antigen',
            'sepsis', 'infection', 'viral', 'virus',
        ]),
        ('Apoptosis / Cell Death', [
            'apoptosis', 'apoptotic', 'caspase', 'cell death',
            'autophagy', 'ferroptosis', 'necrosis', 'trail',
        ]),
        ('Signaling', [
            'signal', 'pathway', 'mapk', 'jak', 'stat',
            'pi3k', 'akt', 'mtor', 'wnt', 'notch', 'hedgehog',
            'vegf', 'egfr', 'erbb', 'pdgf', 'tgf', 'bmp',
            'ras', 'raf', 'receptor', 'kinase',
            'prolactin', 'rank', 'tnf', 'oncostatin',
            'bdnf', 'ngf', 'fgfr', 'gastrin',
        ]),
        ('Cell Cycle / DNA Repair', [
            'cell cycle', 'dna damage', 'dna repair', 'checkpoint',
            'mitotic', 'meiotic', 'tp53', 'p53', 'rb1',
            'senescence', 'telomer',
        ]),
        ('Metabolism', [
            'metabol', 'biosynth', 'glycol', 'oxida', 'lipid',
            'fatty acid', 'cholesterol', 'steroid', 'amino acid',
            'glucose', 'tca', 'citric', 'vitamin',
        ]),
        ('Development / Differentiation', [
            'develop', 'differentiat', 'stem cell', 'embryo',
            'morphogen', 'neurogenesis', 'hematopoie',
        ]),
    ]
    
    for category, keywords in categories:
        for kw in keywords:
            if kw in name_lower:
                return category
    
    return 'Other'

print("✓ classify_pathway 関数を定義しました")""")

code("""# === Treemap用データの構築（全パスウェイ表示） ===
df = results_df.sort_values('Adjusted P-value').copy()
top_n = len(df)  # 全パスウェイを使用

# データ整形
df['Hit_Count'] = df['Overlap'].apply(lambda x: int(x.split('/')[0]))
df['Pathway_Size'] = df['Overlap'].apply(lambda x: int(x.split('/')[1]))
df['-log10_P'] = -np.log10(df['Adjusted P-value'].clip(lower=1e-50))

# パスウェイ名を短縮
def shorten_name(name, max_len=40):
    for suffix in ['WP', ' Homo sapiens']:
        idx = name.find(suffix)
        if idx > 0:
            name = name[:idx].strip().rstrip('_').rstrip()
    if len(name) > max_len:
        return name[:max_len-2] + '..'
    return name

df['PathwayName'] = df['Term'].apply(shorten_name)
df['Category'] = df['Term'].apply(classify_pathway)

# ホバーテキスト作成
df['Hover'] = df.apply(
    lambda r: f"<b>{r['PathwayName']}</b><br>"
              f"Overlap: {r['Overlap']}<br>"
              f"P-value: {r['P-value']:.2e}<br>"
              f"FDR: {r['Adjusted P-value']:.2e}<br>"
              f"-log10(FDR): {r['-log10_P']:.1f}<br>"
              f"Genes: {r['Genes']}",
    axis=1
)

print(f"Top {len(df)} パスウェイを処理")
print(f"カテゴリ分布:")
for cat, cnt in df['Category'].value_counts().items():
    print(f"  {cat}: {cnt}")""")

code("""# === 階層構造の構築 (Root → Category → Pathway) ===
# Plotly go.Treemap用のデータ構造を作成

ids = []
parents = []
labels = []
values = []
hover_texts = []
log_p_values = []

# Root ノード（空文字列が親＝最上位）
ids.append('WikiPathways')
parents.append('')
labels.append('WikiPathways Enrichment')
values.append(0)  # 親ノードは0（子の合計が使われる）
hover_texts.append('WikiPathways Enrichment Analysis')
log_p_values.append(0)

# カテゴリノード
categories_in_data = df['Category'].unique()
for cat in categories_in_data:
    ids.append(cat)
    parents.append('WikiPathways')
    labels.append(cat)
    values.append(0)
    cat_df = df[df['Category'] == cat]
    hover_texts.append(f"{cat}<br>{len(cat_df)} pathways")
    log_p_values.append(0)  # 後で Max_Descendant_P で上書き

# パスウェイノード（リーフ）
for _, row in df.iterrows():
    pw_id = row['Term']  # ユニークID
    ids.append(pw_id)
    parents.append(row['Category'])  # 親はカテゴリ
    labels.append(row['PathwayName'])
    values.append(row['Hit_Count'])
    hover_texts.append(row['Hover'])
    log_p_values.append(row['-log10_P'])

# DataFrameに変換
df_tree = pd.DataFrame({
    'PathwayID': ids,
    'Parent': parents,
    'PathwayName': labels,
    'Gene_Count': values,
    'Hover': hover_texts,
    '-log10_P': log_p_values,
})

print(f"Treemapノード数: {len(df_tree)}")
print(f"  Root: 1, カテゴリ: {len(categories_in_data)}, パスウェイ: {len(df)}")""")

code("""# === Max_Descendant_P の計算 ===
# 子 → 親マッピングを作成
child_to_parent = dict(zip(df_tree['PathwayID'], df_tree['Parent']))

# 親 → 子の構造を作る
parent_to_children = defaultdict(list)
for child, parent in child_to_parent.items():
    parent_to_children[parent].append(child)

# 子孫全体を再帰的に収集
def get_descendants(node_id, tree):
    descendants = set()
    stack = [node_id]
    while stack:
        current = stack.pop()
        children = tree.get(current, [])
        for child in children:
            if child not in descendants:
                descendants.add(child)
                stack.append(child)
    return descendants

# 各ノードの子孫中の最大 -log10_P を計算
pval_dict = dict(zip(df_tree['PathwayID'], df_tree['-log10_P']))
max_pvals = []

for node_id in df_tree['PathwayID']:
    descendants = get_descendants(node_id, parent_to_children)
    all_relevant = [pval_dict[node_id]] + [pval_dict.get(d, 0) for d in descendants]
    max_pvals.append(max(all_relevant))

df_tree['Max_Descendant_P'] = max_pvals

print("✓ Max_Descendant_P を計算しました")
print(f"  最大値: {df_tree['Max_Descendant_P'].max():.1f}")
print(f"  最小値（>0）: {df_tree[df_tree['Max_Descendant_P'] > 0]['Max_Descendant_P'].min():.1f}")""")

code("""# === Plotly Treemap の描画 ===
fig = go.Figure(go.Treemap(
    ids=df_tree['PathwayID'],
    parents=df_tree['Parent'],
    values=df_tree['Gene_Count'],
    labels=df_tree['PathwayName'],
    text=df_tree['Hover'],
    textinfo='label',
    hovertext=df_tree['Hover'],
    hoverinfo='text',
    marker=dict(
        colors=df_tree['Max_Descendant_P'],
        colorscale='Blues',
        cmin=0,
        colorbar=dict(title='-log₁₀(FDR)<br>Max in descendants')
    )
))

fig.update_traces(marker_line_width=0.5, marker_line_color='black')

fig.update_layout(
    title=dict(
        text=f'WikiPathways エンリッチメント解析 — 階層的 Treemap (Top {top_n} pathways)',
        font=dict(size=16)
    ),
    width=800,
    height=550,
    margin=dict(t=50, l=0, r=0, b=0)
)

# HTMLファイルとして保存
fig.write_html('treemap_hierarchical.html')
print("✓ treemap_hierarchical.html に保存しました")

fig.show()""")

code("""# === カテゴリ別サマリー ===
print("=== カテゴリ別サマリー ===")
for cat in sorted(df['Category'].unique()):
    cat_df = df[df['Category'] == cat]
    total_hits = cat_df['Hit_Count'].sum()
    best_fdr = cat_df['Adjusted P-value'].min()
    print(f"  {cat}: {len(cat_df)} pathways, {total_hits} hits, best FDR={best_fdr:.2e}")""")

# ===== Section 5: Pathway Gene Mapping =====
md("""---
## 5. パスウェイ上への遺伝子マッピング

選択したパスウェイのダイアグラム（SVG）をWikiPathwaysからダウンロードし、
ヒットした遺伝子を**赤色**でハイライト表示します。""")

code("""import xml.etree.ElementTree as ET
from IPython.display import SVG, display, HTML

def extract_wp_id(term_name):
    \"\"\"パスウェイ名からWP IDを抽出 (例: 'Apoptosis WP254' -> 'WP254')\"\"\"
    match = re.search(r'(WP\\d+)', term_name)
    return match.group(1) if match else None


def download_pathway_svg(wp_id, save_dir='data/pathways'):
    \"\"\"WikiPathwaysからパスウェイSVGをダウンロード\"\"\"
    os.makedirs(save_dir, exist_ok=True)
    filepath = os.path.join(save_dir, f'{wp_id}.svg')
    
    if os.path.exists(filepath):
        with open(filepath, 'r', encoding='utf-8') as f:
            return f.read()
    
    url = f'https://www.wikipathways.org/wikipathways-assets/pathways/{wp_id}/{wp_id}.svg'
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(resp.text)
        print(f'✓ SVGダウンロード完了: {wp_id}')
        return resp.text
    except Exception as e:
        print(f'✗ SVGダウンロード失敗 ({wp_id}): {e}')
        return None


def highlight_genes_in_svg(svg_text, hit_genes, highlight_color='#FF6666'):
    \"\"\"SVG内のヒット遺伝子をハイライト（赤色に変更）\"\"\"
    hit_genes_upper = {g.upper() for g in hit_genes}
    highlighted = set()
    
    # SVGのXML名前空間を登録
    namespaces = {
        'svg': 'http://www.w3.org/2000/svg',
        'xlink': 'http://www.w3.org/1999/xlink',
    }
    for prefix, uri in namespaces.items():
        ET.register_namespace(prefix, uri)
    ET.register_namespace('', 'http://www.w3.org/2000/svg')
    
    root = ET.fromstring(svg_text)
    
    # class属性にGeneProductを含む要素を検索
    for elem in root.iter():
        cls = elem.get('class', '')
        if 'GeneProduct' not in cls:
            continue
        
        # クラス名からHGNC遺伝子シンボルを抽出
        class_parts = cls.split()
        gene_name = None
        for part in class_parts:
            if part.startswith('HGNC_'):
                gene_name = part.replace('HGNC_', '')
                break
        # HGNCが無い場合はクラス名から直接探す
        if not gene_name:
            for part in class_parts:
                if part.upper() in hit_genes_upper:
                    gene_name = part
                    break
        
        if gene_name and gene_name.upper() in hit_genes_upper:
            # この要素の子のrect要素のfillを変更
            for child in elem.iter():
                tag = child.tag.split('}')[-1] if '}' in child.tag else child.tag
                if tag == 'rect':
                    child.set('fill', highlight_color)
                    child.set('fill-opacity', '0.7')
            highlighted.add(gene_name)
    
    # SVG文字列に変換
    modified_svg = ET.tostring(root, encoding='unicode')
    return modified_svg, highlighted


def visualize_pathway(term_name, hit_gene_str, highlight_color='#FF6666'):
    \"\"\"パスウェイをダウンロードしてヒット遺伝子をハイライト表示\"\"\"
    wp_id = extract_wp_id(term_name)
    if not wp_id:
        print(f'✗ WP IDが見つかりません: {term_name}')
        return
    
    # ヒット遺伝子のリスト
    hit_genes = [g.strip() for g in hit_gene_str.split(';')]
    
    print(f'パスウェイ: {term_name}')
    print(f'WP ID: {wp_id}')
    genes_str = ", ".join(hit_genes)
    print(f'ヒット遺伝子 ({len(hit_genes)}): {genes_str}')
    print()
    
    # SVGダウンロード
    svg_text = download_pathway_svg(wp_id)
    if not svg_text:
        return
    
    # 遺伝子ハイライト
    modified_svg, highlighted = highlight_genes_in_svg(svg_text, hit_genes, highlight_color)
    
    hl_str = ", ".join(sorted(highlighted))
    print(f'ハイライトされた遺伝子 ({len(highlighted)}/{len(hit_genes)}): {hl_str}')
    not_found = set(g.upper() for g in hit_genes) - {g.upper() for g in highlighted}
    if not_found:
        nf_str = ", ".join(sorted(not_found))
        print(f'SVG内に見つからなかった遺伝子: {nf_str}')
    
    # ハイライト済みSVGを保存
    output_path = f'data/pathways/{wp_id}_highlighted.svg'
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(modified_svg)
    print(f'✓ 保存: {output_path}')
    
    # 凡例を表示
    display(HTML(f'''
    <div style="margin: 10px 0; padding: 8px; border: 1px solid #ddd; border-radius: 5px; background: #fafafa;">
        <b>凡例:</b>
        <span style="background: {highlight_color}; padding: 2px 8px; border-radius: 3px; margin-left: 10px;">■ ヒット遺伝子</span>
        <span style="margin-left: 15px;">□ 非ヒット遺伝子</span>
    </div>
    '''))
    
    # SVG表示
    display(SVG(data=modified_svg))

print('✓ パスウェイ可視化関数を定義しました')""")

md("""### パスウェイの選択と可視化

最も有意なパスウェイを自動選択して可視化します。  
別のパスウェイを可視化したい場合は、`pathway_index` を変更してください。""")

code("""# === 最も有意なパスウェイを可視化 ===
# 補正P値でソート済みの結果から選択
sorted_results = results_df.sort_values('Adjusted P-value')

# ★ 別のパスウェイを可視化する場合は index を変更 ★
pathway_index = 0  # 0 = 最も有意なパスウェイ

if len(sorted_results) > pathway_index:
    row = sorted_results.iloc[pathway_index]
    print(f"=== パスウェイ #{pathway_index + 1} (FDR = {row['Adjusted P-value']:.2e}) ===")
    visualize_pathway(row['Term'], row['Genes'])
else:
    print('有意なパスウェイがありません。')""")

code("""# === 2番目に有意なパスウェイも可視化（任意） ===
pathway_index = 1

if len(sorted_results) > pathway_index:
    row = sorted_results.iloc[pathway_index]
    print(f"=== パスウェイ #{pathway_index + 1} (FDR = {row['Adjusted P-value']:.2e}) ===")
    visualize_pathway(row['Term'], row['Genes'])
else:
    print('2番目のパスウェイがありません。')""")


# ===== Section 6: Summary =====
md("""---
## 6. 結果のサマリー & エクスポート""")

code("""# 結果をCSVにエクスポート
output_csv = "enrichment_results.csv"
results_df.to_csv(output_csv, index=False)
print(f"✓ 結果をCSVに保存しました: {output_csv}")

# サマリー統計
print(f"\\n=== 解析サマリー ===")
print(f"入力遺伝子数: {len(gene_list)}")
print(f"使用GMTファイル: {ENRICHR_LIBRARY}")
print(f"公式GMTファイル（ローカル保存）: {os.path.basename(official_gmt)}")
print(f"パスウェイデータベース: WikiPathways ({SPECIES})")
print(f"テスト済みパスウェイ数: {len(pathways)}")
print(f"有意なパスウェイ数 (FDR < 0.05): {len(results_df)}")
if not results_df.empty:
    print(f"\\n最も有意なパスウェイ:")
    print(f"  名前: {results_df.iloc[0]['Term']}")
    print(f"  FDR P-value: {results_df.iloc[0]['Adjusted P-value']:.2e}")
    print(f"  Overlap: {results_df.iloc[0]['Overlap']}")""")

md("""---
## 使い方ガイド

1. **遺伝子リストの変更**: セクション3の `gene_list` を自分のデータに差し替えてください
2. **生物種の変更**: セクション2の `SPECIES` と `ENRICHR_LIBRARY` を変更してください
   - マウス: `SPECIES = "Mus_musculus"`, `ENRICHR_LIBRARY = "WikiPathways_2024_Mouse"`
3. **閾値の変更**: `gp.enrich()` の `cutoff` パラメータでFDR閾値を調整できます
4. **Treemapの件数**: Treemapは全パスウェイを表示します
5. **パスウェイ可視化**: セクション5の `pathway_index` を変更して別のパスウェイを可視化できます
6. **背景遺伝子セット**: `gp.enrich()` に `background=20000` 等を追加して背景を指定できます""")

# ===== Write notebook =====
nb = {
    "nbformat": 4,
    "nbformat_minor": 5,
    "metadata": {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3"
        },
        "language_info": {
            "name": "python",
            "version": "3.11.0"
        }
    },
    "cells": cells
}

outpath = "/Users/yoshinorisatomi/Documents/antigravity/wikipathway/pj01/wikipathways_enrichment.ipynb"
with open(outpath, 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print(f"Notebook created: {outpath}")
