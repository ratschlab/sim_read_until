# Data preparation

Take data from an existing run and put it into the data directory with
```{bash}
mkdir -p tests/simulator/data/sim_reads
sed -n '1,10007p;10008q' runs/replicate_run/reads.fasta > tests/simulator/data/sim_reads/reads1.fasta
sed -n '10008,20006p;20007q' runs/replicate_run/reads.fasta > tests/simulator/data/sim_reads/reads2.fasta
```
