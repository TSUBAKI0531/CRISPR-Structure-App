import re
from Bio.Seq import Seq

class CRISPRDesigner:
    def __init__(self, target_dna: str):
        self.target_dna = target_dna.upper()

    def design_grna(self):
        candidates = []
        seq_len = len(self.target_dna)

        # センス鎖 (+)
        for i in range(seq_len - 23):
            sub = self.target_dna[i:i+23]
            if sub.endswith("GG"):
                grna = sub[:20]
                candidates.append(self._create_candidate(grna, sub[20:], i + 17, "Sense (+)"))

        # アンチセンス鎖 (-)
        for i in range(seq_len - 23):
            sub = self.target_dna[i:i+23]
            if sub.startswith("CC"):
                grna = str(Seq(sub[3:]).reverse_complement())
                candidates.append(self._create_candidate(grna, sub[:3], i + 6, "Antisense (-)"))
        
        return sorted(candidates, key=lambda x: x['score'], reverse=True)

    def _create_candidate(self, grna, pam, cut_site, strand):
        gc = (grna.count('G') + grna.count('C')) / 20 * 100
        off_target_count = self._predict_off_target(grna) # オフターゲット予測
        
        return {
            "strand": strand,
            "grna": grna,
            "pam": pam,
            "cut_site": cut_site,
            "gc": f"{gc:.1f}%",
            "off_target_sites": off_target_count,
            "score": self._calculate_score(gc, off_target_count)
        }

    def _predict_off_target(self, grna):
        """
        ゲノム全体への類似配列検索をシミュレーション。
        実際にはBLAST等のツールが必要ですが、ここでは配列の連続性から
        リスク(0~5箇所のヒット)を擬似的に算出します。
        """
        import random
        # 特定のモチーフ(TTTTなど)があるとオフターゲットリスクが上がると仮定
        risk_base = grna.count("AAAA") + grna.count("TTTT")
        return min(5, risk_base + random.randint(0, 1))

    def _calculate_score(self, gc, off_target):
        # GC40-60%が理想、かつオフターゲットが少ないほど高スコア
        base_score = 100 - abs(50 - gc) * 2
        penalty = off_target * 15
        return max(0, int(base_score - penalty))

    def simulate_knock_in(self, cut_site, insert_seq):
        return self.target_dna[:cut_site] + insert_seq.upper() + self.target_dna[cut_site:]

    def translate_to_protein(self, dna_seq):
        trim_len = (len(dna_seq) // 3) * 3
        return str(Seq(dna_seq[:trim_len]).translate(to_stop=False))