from fpdf import FPDF
from datetime import datetime

class ReportGenerator:
    def __init__(self):
        self.pdf = FPDF()
        self.pdf.add_page()
        
    def generate(self, grna_info, orig_prot, edit_prot, rmsd):
        self.pdf.set_font("Arial", 'B', 16)
        self.pdf.cell(200, 10, txt="CRISPR Advanced Analysis Report", ln=True, align='C')
        
        self.pdf.set_font("Arial", 'B', 12)
        self.pdf.ln(10)
        self.pdf.cell(200, 10, txt="1. Selection Metrics", ln=True)
        self.pdf.set_font("Arial", size=10)
        self.pdf.multi_cell(0, 8, txt=(
            f"- gRNA: {grna_info['grna']}\n"
            f"- RMSD: {rmsd:.4f} Angstrom\n"
            f"- Predicted Off-target Sites: {grna_info['off_target_sites']}\n"
            f"- GC Content: {grna_info['gc']}\n"
            f"- Composite Score: {grna_info['score']}"
        ))
        
        # 配列情報などは以前と同様
        self.pdf.ln(5)
        self.pdf.set_font("Arial", 'B', 12)
        self.pdf.cell(200, 10, txt="2. Protein Sequences", ln=True)
        self.pdf.set_font("Courier", size=8)
        self.pdf.multi_cell(0, 5, txt=f"Original:\n{orig_prot}\n\nEdited:\n{edit_prot}")

        return self.pdf.output(dest='S').encode('latin-1')