"""
Pore model to convert DNA sequence to ONT raw signal
"""
import numpy as np
import pandas as pd

class PoreModel:
    """
    ONT's rudimentary pore model to convert DNA sequence to a raw signal using mean pore voltage per k-mer
    
    Example:
        # mkdir -p runs/data/pore_models/legacy_r9.4_180mv_450bps_6mer
        # wget https://raw.githubusercontent.com/nanoporetech/kmer_models/4e56daed7fbb79b538f58e41262d5c54b07356ea/legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model -O runs/data/pore_models/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model
        pore_filename = Path("runs/data/pore_models/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model")
        pore_model = PoreModel(pore_filename, signals_per_bp=max(1, round(4000 / 450)))
        pore_model.to_raw("AAAACCTGGGC")
        
    Args:
        pore_filename: path to the pore model file
        signals_per_bp: how many signals per bp to generate, e.g., max(1, round(4000 / bp_speed))
    """
    def __init__(self, pore_filename, signals_per_bp=9):
        df = pd.read_csv(pore_filename, sep="\t")
        # could also take level_stdv, sd_mean, sd_stdv into account, as described here: https://github.com/nanoporetech/kmer_models/tree/4e56daed7fbb79b538f58e41262d5c54b07356ea/legacy
        self.kmer_to_mean = dict(zip(df["kmer"], df["level_mean"]))
        
        self.signals_per_bp = signals_per_bp
        self.k = len(df.iloc[0]["kmer"])
        
    # @classmethod
    # def from_pore_filename(cls, pore_filename):
    #     pore_params = cls.get_pore_params(pore_filename)
    #     return cls(pore_filename, **pore_params)
    # @staticmethod
    # def get_pore_params(pore_filename):
    #     pore_filename = Path(pore_filename)
    #     _, _, bias, bp_speed, k = pore_filename.parent.name.split("_")
    #     bias = float(bias.replace("mv", ""))
    #     bp_speed = float(bp_speed.replace("bps", ""))
    #     k = int(k.replace("mer", ""))
    #     return {"bp_speed": bp_speed, "k": k}
    
    def to_raw(self, seq: str):
        """
        Convert a sequence to raw signal
        
        Since each bp corresponds to several signals, its value is duplicated.
        The voltage bias is not added because methods usually normalize the signal. It is also not zscore normalized.
        """
        return np.array([self.kmer_to_mean[seq[i:i+self.k]] for i in range(len(seq) - self.k + 1) for _ in range(self.signals_per_bp)])
    