import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
from datetime import datetime
from modules.crispr_designer import CRISPRDesigner
from modules.protein_analyzer import ProteinStructureAnalyzer
from modules.report_generator import ReportGenerator

# ページ設定
st.set_page_config(page_title="CRISPR Analysis Pro", layout="wide")
st.title("🧬 Advanced CRISPR Screening Pipeline")

# --- 1. セッション状態の初期化 ---
if 'analyzed' not in st.session_state:
    st.session_state.analyzed = False
    st.session_state.results = None
    st.session_state.orig_prot = None
    st.session_state.p1_path = None

# サイドバー：入力
with st.sidebar:
    st.header("1. Parameters")
    default_dna = "ATGCAGATCTTCGTGAAGACCCTGACCGGCAAGACCATCACCCTGGAGGTGGAGCCCAGTGACACCATCGAGAATGTGAAGGCCAAGATCCAGGATAAGGAGGGCATCCCCCCTGACCAGCAGAGGCTGATCTTTGCCGGCAAGCAGCTGGAGGATGGCCGCACCCTGTCTGATTACAACATCCAGAAGGAGTCCACCCTGCACCTGGTCCTCCGTCTGAGGGGT"
    target_dna = st.text_area("Target DNA Sequence", default_dna, height=150)
    insert_seq = st.text_input("Insertion Sequence", "GGCAGCGGC")
    st.markdown("---")
    top_n = st.slider("Analyze top N candidates", 1, 10, 3)
    rmsd_threshold = st.number_input("RMSD Warning Threshold (Å)", value=2.0, step=0.5)

# --- 2. 解析実行ロジック ---
if st.button("🚀 Run Comprehensive Screening"):
    designer = CRISPRDesigner(target_dna)
    analyzer = ProteinStructureAnalyzer()
    candidates = designer.design_grna()
    
    if candidates:
        with st.spinner("Analyzing high-throughput data..."):
            target_candidates = candidates[:top_n]
            results = []
            orig_prot = designer.translate_to_protein(target_dna)
            p1 = analyzer.fetch_structure(orig_prot, "orig")
            
            for idx, cand in enumerate(target_candidates):
                edit_prot = designer.translate_to_protein(designer.simulate_knock_in(cand['cut_site'], insert_seq))
                p2 = analyzer.fetch_structure(edit_prot, f"edit_{idx}")
                rmsd = analyzer.calculate_rmsd(p1, p2) if p1 and p2 else None
                
                res_item = cand.copy()
                res_item["rmsd"] = rmsd
                res_item["edit_prot"] = edit_prot # レポート用に保持
                results.append(res_item)

            # 解析結果をセッションに保存
            st.session_state.results = results
            st.session_state.orig_prot = orig_prot
            st.session_state.p1_path = p1
            st.session_state.analyzed = True
    else:
        st.error("No candidates found.")

# --- 3. 解析完了後の表示ロジック ---
# ボタンの中ではなく、セッション状態を見て表示を切り替える
if st.session_state.analyzed:
    st.subheader("1. Screening Results & Ranking")
    
    df_res = pd.DataFrame(st.session_state.results).sort_values("rmsd")
    
    # 警告ステータス判定
    def get_status(row):
        if row['rmsd'] > rmsd_threshold: return "⚠️ High Impact"
        if row['off_target_sites'] > 2: return "❗ Off-target Risk"
        return "✅ Safe"

    # 前回のKeyErrorを修正した表示用DF作成
    display_df = (
        df_res.assign(
            rank=range(1, len(df_res) + 1),
            status=df_res.apply(get_status, axis=1)
        )
        [["rank", "status", "grna", "rmsd", "off_target_sites", "gc", "score"]]
        .set_index("rank")
    )
    st.dataframe(display_df, use_container_width=True)

    # 2. 詳細可視化 (ここを操作しても解析が消えなくなります)
    st.subheader("2. Structural Visualizer")
    selected_rank = st.selectbox(
        "Select candidate to visualize:", 
        range(len(df_res)), 
        format_func=lambda x: f"Rank {x+1}: {df_res.iloc[x]['grna']} (RMSD: {df_res.iloc[x]['rmsd']:.4f})"
    )
    
    selected_cand = df_res.iloc[selected_rank]
    analyzer = ProteinStructureAnalyzer()
    
    # 選択された候補の構造を読み込み
    p_edit = analyzer.fetch_structure(selected_cand['edit_prot'], f"edit_display_{selected_rank}")
    
    col1, col2 = st.columns([1, 2])
    with col1:
        st.metric("RMSD Value", f"{selected_cand['rmsd']:.4f} Å")
        if selected_cand['rmsd'] > rmsd_threshold:
            st.error("Structure highly affected.")
        
        # PDFレポート生成
        try:
            reporter = ReportGenerator()
            pdf = reporter.generate(selected_cand, st.session_state.orig_prot, selected_cand['edit_prot'], selected_cand['rmsd'])
            st.download_button(f"📥 Download Report (Rank {selected_rank+1})", pdf, f"Analysis_Rank{selected_rank+1}.pdf", "application/pdf")
        except:
            st.warning("PDF generator needs compatible characters.")

    with col2:
        # HTML/JSによる3D表示
        components.html(analyzer.get_html_view(st.session_state.p1_path, p_edit), height=550)