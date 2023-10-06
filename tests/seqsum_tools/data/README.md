# Creating data

From a simulated run (see usecase)
```{python}
# alternatively run with less channels
from Bio import SeqIO
from simreaduntil.simulator.channel_element import ReadDescriptionParser, ReadTags, end_reason_to_ont_map

# filter sequencing summary to a few channels
full_reads_filename = "runs/enrich_usecase/simulator_run/reads/reads.fasta"
channels = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# keep all fasta records coming from a channel in channels
with open("tests/param_extractor/data/simreads.fasta", "w") as f:
    SeqIO.write(
        (record for record in SeqIO.parse(full_reads_filename, "fasta") if int(ReadDescriptionParser(record.description.split(" ", maxsplit=1)[1]).ch) in channels),
        f, "fasta"
    )
!simfasta_to_seqsum tests/param_extractor/data/simreads.fasta --seqsummary_filename tests/param_extractor/data/sim_sequencing_summary.txt
```

From a real run
```{python}
sequencing_summary_file = "/Users/maximilianmordig/Desktop/ont_sequencing/data/sequencing_summaries/uncalled/20190809_zymo_seqsum.txt" # from UNCALLED
df_read = pd.read_csv(sequencing_summary_file, sep="\t")
sub_seqsum_df = df_read[(df_read["start_time"] < 13000) & (df_read["channel"] < 10)]
sub_seqsum_df["end_reason"] = "signal_positive"
sub_seqsum_df.to_csv("tests/param_extractor/data/zymo_short_seqsum.txt", sep="\t", index=False)
```
