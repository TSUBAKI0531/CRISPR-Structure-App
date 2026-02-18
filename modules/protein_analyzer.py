import requests
import os

class ProteinStructureAnalyzer:
    def __init__(self, output_dir: str = "data"):
        self.api_url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
        self.output_dir = output_dir
        if not os.path.exists(output_dir): os.makedirs(output_dir)

    def fetch_structure(self, sequence: str, filename: str) -> str:
        file_path = os.path.join(self.output_dir, f"{filename}.pdb")
        # 停止コドン '*' が含まれているとESMFoldがエラーを出すため除去
        clean_seq = sequence.replace("*", "")
        response = requests.post(self.api_url, data=clean_seq, verify=False)
        if response.status_code == 200:
            with open(file_path, "w") as f: f.write(response.text)
            return file_path
        return None

    def calculate_rmsd(self, pdb1: str, pdb2: str) -> float:
        from Bio.PDB import PDBParser, Superimposer
        try:
            parser = PDBParser(QUIET=True)
            s1, s2 = parser.get_structure("s1", pdb1), parser.get_structure("s2", pdb2)
            atoms1 = [res["CA"] for m in s1 for c in m for res in c if "CA" in res]
            atoms2 = [res["CA"] for m in s2 for c in m for res in c if "CA" in res]
            min_len = min(len(atoms1), len(atoms2))
            sup = Superimposer()
            sup.set_atoms(atoms1[:min_len], atoms2[:min_len])
            return sup.rms
        except: return -1.0

    def get_html_view(self, pdb1: str, pdb2: str) -> str:
        """完全自立型JavaScript描画ロジック"""
        with open(pdb1, 'r') as f: d1 = f.read().replace('\n', '\\n').replace("'", "\\'")
        with open(pdb2, 'r') as f: d2 = f.read().replace('\n', '\\n').replace("'", "\\'")

        return f"""
        <div id="container" style="height: 600px; width: 100%; border: 1px solid #ddd; border-radius: 8px;"></div>
        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
        <script>
            function initViewer() {{
                if (typeof $3Dmol === 'undefined') {{ setTimeout(initViewer, 100); return; }}
                const element = document.getElementById('container');
                const viewer = $3Dmol.createViewer(element, {{ backgroundColor: 'white' }});
                viewer.addModel('{d1}', "pdb");
                viewer.setStyle({{model: 0}}, {{cartoon: {{color: 'blue', opacity: 0.7}}}});
                viewer.addModel('{d2}', "pdb");
                viewer.setStyle({{model: 1}}, {{cartoon: {{color: 'red'}}}});
                viewer.zoomTo();
                viewer.render();
            }}
            initViewer();
        </script>
        """