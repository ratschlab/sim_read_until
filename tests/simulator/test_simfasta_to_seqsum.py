import pandas as pd
from simreaduntil.seqsum_tools.seqsum_preprocessing import sort_and_clean_seqsum_df
from simreaduntil.simulator.simfasta_to_seqsum import convert_simfasta_dir_to_seqsum, convert_simfasta_to_seqsum


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
    