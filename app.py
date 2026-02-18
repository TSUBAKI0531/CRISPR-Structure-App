import streamlit as st
import streamlit.components.v1 as components
from modules.crispr_designer import CRISPRDesigner
from modules.protein_analyzer import ProteinStructureAnalyzer

# ページ設定
st.set_page_config(page_title="CRISPR Structure App", layout="wide")

st.title("🧬 CRISPR-to-Structure Pipeline")
st.markdown("Designed for **gRNA Design** & **Protein Structural Impact Analysis**")

# サイドバー：入力パラメータ
with st.sidebar:
    st.header("Input Parameters")
    target_dna = st.text_area("Target DNA Sequence", "ATGGCCCCAAACTGAGTACCCTAGGACCGGTTTTAGCCGATCGATCGATCGATCGATCGG", height=150)
    insert_seq = st.text_input("Insertion Sequence (Knock-in)", "GGCAGCGGC")
    
    st.markdown("---")
    st.info("Ensure the DNA sequence length is sufficient for gRNA search.")

# メイン処理
if st.button("🚀 Run Analysis"):
    # 1. gRNA設計
    st.subheader("1. gRNA Design & Selection")
    designer = CRISPRDesigner(target_dna)
    candidates = designer.design_grna()
    
    if not candidates:
        st.error("No suitable gRNA found. Please check PAM (NGG) or sequence length.")
    else:
        # 結果をデータフレームで表示
        st.dataframe(candidates)
        
        # スコアが一番高いものを自動選択
        best_grna = candidates[0]
        st.success(f"Best gRNA Selected: `{best_grna['grna']}` (Strand: {best_grna['strand']})")
        
        # 2. 編集シミュレーション
        st.subheader("2. Sequence Editing & Translation")
        edited_dna = designer.simulate_knock_in(best_grna['cut_site'], insert_seq)
        
        orig_prot = designer.translate_to_protein(designer.target_dna)
        edit_prot = designer.translate_to_protein(edited_dna)
        
        col1, col2 = st.columns(2)
        with col1:
            st.text_area("Original Protein", orig_prot, height=100)
        with col2:
            st.text_area("Edited Protein", edit_prot, height=100)
            
        # 3. 構造予測と解析
        st.subheader("3. Structural Prediction (ESMFold) & RMSD")
        analyzer = ProteinStructureAnalyzer()
        
        with st.spinner("Predicting 3D structures via ESMFold API... (This may take a minute)"):
            p1 = analyzer.fetch_structure(orig_prot, "original_prot")
            p2 = analyzer.fetch_structure(edit_prot, "edited_prot")
        
        if p1 and p2:
            rmsd = analyzer.calculate_rmsd(p1, p2)
            st.metric("RMSD Value", f"{rmsd:.4f} Å", delta_color="inverse")
            
            if rmsd < 2.0:
                st.caption("✅ Minor structural change detected.")
            else:
                st.caption("⚠️ Significant structural change detected.")
            
            # 3D可視化
            st.markdown("### 3D Visualization")
            mol_html = analyzer.render_mol_html(p1, p2)
            components.html(mol_html, height=600, width=800)
        else:
            st.error("Failed to fetch structures from ESMFold API.")