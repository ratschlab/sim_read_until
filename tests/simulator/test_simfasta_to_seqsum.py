from pathlib import Path
from textwrap import dedent
import pandas as pd
from simreaduntil.seqsum_tools.seqsum_preprocessing import sort_and_clean_seqsum_df
from simreaduntil.simulator.channel_element import ChunkedRead
from simreaduntil.simulator.simfasta_to_seqsum import SequencingSummaryWriter, convert_simfasta_dir_to_seqsum, convert_simfasta_to_seqsum, write_seqsum_header, write_seqsum_record_line
from Bio import SeqIO
from Bio.Seq import Seq

def get_dummy_record():
    # >chr20_36784526_aligned_proc0:16m_R_0_2226_0 full_seqlen=13552 t_start=0.1752480000000105 t_end=6.115326881408691 t_delay=0.09375 ended=user_unblocked tags= full_read_id=chr20_36772816_aligned_proc0:16_R_0_13936_0 ch=ch138
    # GTGCAATTTATACTCATGGCCAGTGTACAGTGACTCATGCCTGTACCCCACTTTAGGAGA
    description = "chr20_36784526_aligned_proc0:16m_R_0_2226_0 full_seqlen=13552 t_start=0.1752480000000105 t_end=6.115326881408691 t_delay=0.09375 ended=user_unblocked tags= full_read_id=chr20_36772816_aligned_proc0:16_R_0_13936_0 ch=ch138"
    read_id = description.split(" ")[0]
    seq = "GTGCAATTTATACTCATGGCCAGTGTACAGTGACTCATGCCTGTACCCCACTTTAGGAGA"
    return SeqIO.SeqRecord(id=read_id, description=description, seq=Seq(seq))
    
# to match with get_dummy_record
expected_seqsum_file_content=dedent(f"""\
read_id\tchannel\tmux\tstart_time\tduration\tpasses_filtering\ttemplate_start\ttemplate_duration\tsequence_length_template\tend_reason\tnb_ref_bps_full\tstopped_receiving\tnever_requested
chr20_36784526_aligned_proc0:16m_R_0_2226_0\tch138\t1\t0.1752480000000105\t5.940078881408681\tTrue\t0.2689980000000105\t5.846328881408681\t{len(get_dummy_record().seq)}\tdata_service_unblock_mux_change\t13936\tFalse\tFalse
""")

def test_write_seqsum_record_line(tmp_path):
    sequencing_summary_filename = tmp_path / "seqsummary_filename_simple1.txt"
    with open(sequencing_summary_filename, mode="w") as seqsummary_fh:
        write_seqsum_header(seqsummary_fh)
        
        write_seqsum_record_line(get_dummy_record(), seqsummary_fh)
        
    assert sequencing_summary_filename.read_text() == expected_seqsum_file_content
    
def test_SequencingSummaryWriter(tmp_path):
    sequencing_summary_filename = tmp_path / "seqsummary_filename_simple2.txt"
    with open(sequencing_summary_filename, mode="w") as seqsummary_fh:
        with SequencingSummaryWriter(seqsummary_fh) as seqsum_writer:
            seqsum_writer.write_read(get_dummy_record())
        
    assert sequencing_summary_filename.read_text() == expected_seqsum_file_content

def test_seqsum_line_with_chunked_read(tmp_path):
    chunked_read = ChunkedRead("read1", "111112222222222333333", 10.1, read_speed=10, min_chunk_size=4)
    chunked_read.get_new_samples(10.1+0.5)
    chunked_read.stop_receiving()
    seq_record = chunked_read.finish()
    seq_record.description += f" ch=ch1" # added by channel on top
    
    sequencing_summary_filename = tmp_path / "seqsummary_filename_simple3.txt"
    with open(sequencing_summary_filename, mode="w") as seqsummary_fh:
        write_seqsum_record_line(seq_record, seqsummary_fh, read_id=seq_record.id)

    # test with SequencingSummaryWriter
    sequencing_summary_filename = tmp_path / "seqsummary_filename_simple4.txt"
    with open(sequencing_summary_filename, mode="w") as seqsummary_fh:
        with SequencingSummaryWriter(seqsummary_fh) as seqsum_writer:
            seqsum_writer.write_read(seq_record)
    
    expected_filecontent = dedent("""\
        read_id\tchannel\tmux\tstart_time\tduration\tpasses_filtering\ttemplate_start\ttemplate_duration\tsequence_length_template\tend_reason\tnb_ref_bps_full\tstopped_receiving\tnever_requested
        read1\tch1\t1\t10.1\t2.0999999999999996\tTrue\t10.1\t2.0999999999999996\t21\tsignal_positive\tnan\tTrue\tFalse
        """)
    assert Path(sequencing_summary_filename).read_text() == expected_filecontent

def test_convert_simfasta_to_seqsum(shared_datadir, tmp_path):
    number_of_lines_in_file = lambda filename: sum(1 for _ in open(filename))
    number_of_lines_starting_with = lambda filename, val: sum(x.startswith(val) for x in open(filename))
    
    sequencing_summary_file = tmp_path / "seqsummary_filename1.txt"
    convert_simfasta_to_seqsum(shared_datadir / "sim_reads" / "reads1.fasta", sequencing_summary_file)
    assert number_of_lines_in_file(sequencing_summary_file) - 1 == number_of_lines_starting_with(shared_datadir / "sim_reads" / "reads1.fasta", ">")
    
    # test it can be read and processed
    df_read = pd.read_csv(sequencing_summary_file, sep="\t")
    seqsum_df = sort_and_clean_seqsum_df(df_read, min_gap_duration=0.05)
    
    # should overwrite sequencing_summary_file
    convert_simfasta_dir_to_seqsum(shared_datadir / "sim_reads", sequencing_summary_file)
    assert number_of_lines_in_file(sequencing_summary_file) - 1 == (
        number_of_lines_starting_with(shared_datadir / "sim_reads" / "reads1.fasta", ">")
        + number_of_lines_starting_with(shared_datadir / "sim_reads" / "reads2.fasta", ">")
    )
    # test it can be read and processed
    df_read = pd.read_csv(sequencing_summary_file, sep="\t")
    seqsum_df = sort_and_clean_seqsum_df(df_read, min_gap_duration=0.05)
    