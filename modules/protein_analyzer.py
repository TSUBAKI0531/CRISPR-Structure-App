import requests
import os
from Bio.PDB import PDBParser, Superimposer
import py3dmol

class ProteinStructureAnalyzer:
    """ESMFold APIを用いた構造予測と可視化クラス"""
    
    def __init__(self, output_dir: str = "data"):
        self.api_url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
        self.output_dir = output_dir
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    def fetch_structure(self, sequence: str, filename: str) -> str:
        file_path = os.path.join(self.output_dir, f"{filename}.pdb")
        # キャッシュがあればそれを使う
        if os.path.exists(file_path):
            return file_path
        
        response = requests.post(self.api_url, data=sequence)
        if response.status_code == 200:
            with open(file_path, "w") as f:
                f.write(response.text)
            return file_path
        else:
            return None

    def calculate_rmsd(self, pdb1: str, pdb2: str) -> float:
        try:
            parser = PDBParser(QUIET=True)
            s1 = parser.get_structure("s1", pdb1)
            s2 = parser.get_structure("s2", pdb2)

            atoms1 = [res["CA"] for m in s1 for c in m for res in c if "CA" in res]
            atoms2 = [res["CA"] for m in s2 for c in m for res in c if "CA" in res]

            min_len = min(len(atoms1), len(atoms2))
            sup = Superimposer()
            sup.set_atoms(atoms1[:min_len], atoms2[:min_len])
            return sup.rms
        except Exception:
            return -1.0

    def render_mol_html(self, pdb1: str, pdb2: str) -> str:
        """Streamlitで表示するためのHTML文字列を生成する"""
        with open(pdb1, 'r') as f: d1 = f.read()
        with open(pdb2, 'r') as f: d2 = f.read()

        view = py3dmol.view(width=700, height=500)
        view.addModel(d1, 'pdb')
        view.setStyle({'model': 0}, {'cartoon': {'color': 'blue', 'opacity': 0.6}})
        view.addLabel("Original (Blue)", {'position': {'x':0, 'y':0, 'z':0}, 'useScreen': True, 'backgroundColor': 'blue'})
        
        view.addModel(d2, 'pdb')
        view.setStyle({'model': 1}, {'cartoon': {'color': 'red'}})
        view.addLabel("Edited (Red)", {'position': {'x':0, 'y':20, 'z':0}, 'useScreen': True, 'backgroundColor': 'red'})
        
        view.align({'model': 1}, {'model': 0})
        view.zoomTo()
        
        # HTMLとして書き出し
        return view._make_html()