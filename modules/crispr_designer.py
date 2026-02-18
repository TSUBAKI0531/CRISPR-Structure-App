from typing import List, Dict, Tuple, Optional
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

class CRISPRDesigner:
    """CRISPR-Cas9 gRNA設計および編集シミュレーションクラス"""
    
    def __init__(self, target_dna: str, exon_regions: List[Tuple[int, int]] = None, genome_db: Optional[List[str]] = None):
        self.target_dna = Seq(target_dna.upper())
        self.exon_regions = exon_regions if exon_regions else []
        self.genome_db = genome_db if genome_db else []

    def _is_in_exon(self, pos: int) -> Tuple[bool, str]:
        for start, end in self.exon_regions:
            if start <= pos <= end:
                return True, "Exon"
        return False, "Intron/Intergenic"

    def _calculate_off_target_risk(self, grna_seq: str) -> int:
        hits = 0
        for ref_seq in self.genome_db:
            mismatches = sum(1 for a, b in zip(grna_seq, ref_seq) if a != b)
            if mismatches <= 3:
                hits += 1
        return hits

    def design_grna(self) -> List[Dict]:
        candidates = []
        target_len = len(self.target_dna)
        
        # センス鎖(+)とアンチセンス鎖(-)の探索
        search_configs = [
            ("Sense (+)", range(20, target_len - 2), True),
            ("Antisense (-)", range(0, target_len - 22), False)
        ]

        for strand, loop_range, is_sense in search_configs:
            for i in loop_range:
                # PAM判定 (NGG or CCN)
                if (is_sense and self.target_dna[i+1:i+3] == "GG") or \
                   (not is_sense and self.target_dna[i:i+2] == "CC"):
                    
                    if is_sense:
                        grna_seq = self.target_dna[i-20:i]
                        cut_site = i - 3
                        pam = str(self.target_dna[i:i+3])
                    else:
                        grna_seq = self.target_dna[i+3:i+23].reverse_complement()
                        cut_site = i + 3
                        pam = str(self.target_dna[i:i+3])

                    gc = gc_fraction(grna_seq)
                    if 0.4 <= gc <= 0.6 and "TTTT" not in grna_seq:
                        in_exon, region = self._is_in_exon(cut_site)
                        off_hits = self._calculate_off_target_risk(str(grna_seq))
                        
                        candidates.append({
                            "strand": strand,
                            "grna": str(grna_seq),
                            "pam": pam,
                            "cut_site": cut_site,
                            "region": region,
                            "off_target_risk": off_hits,
                            "score": 100 - (off_hits * 20)
                        })
        return sorted(candidates, key=lambda x: x['score'], reverse=True)

    def simulate_knock_in(self, cut_site: int, insert_seq: str) -> Seq:
        return self.target_dna[:cut_site] + Seq(insert_seq.upper()) + self.target_dna[cut_site:]

    @staticmethod
    def translate_to_protein(dna_seq: Seq) -> str:
        # 簡易的な翻訳（ORF探索なし、ベタ翻訳）
        over = len(dna_seq) % 3
        valid_dna = dna_seq[:-over] if over != 0 else dna_seq
        return str(valid_dna.translate())