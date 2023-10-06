from simreaduntil.simulator.pore_model import PoreModel

def test_pore_model(shared_datadir):
    seq = "ACCCTTTGGG"
    k = 6
    signals_per_bp = 2
    pore_filename = shared_datadir / "dummy_pore_model.csv"
    
    pore_model = PoreModel(pore_filename, signals_per_bp=signals_per_bp)
    raw_signal = pore_model.to_raw(seq)
    assert len(raw_signal) == (len(seq) - k + 1) * signals_per_bp
    raw_signal