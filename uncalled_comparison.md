# Comparison to UNCALLED

We tried running UNCALLED with the following commands, but it segfaults:

```{bash}
# download a tar.gz containing sequencing data, use guppy to basecall it, then input the fast5 directory together with the sequencing summary and the uncalled sequencing summary files to the uncalled simulator
# in big_data dir
mkdir uncalled_sim_exp
cd uncalled_sim_exp

# wget https://sra-download.be-md.ncbi.nlm.nih.gov/vast/sra01/SRZ/013127/SRR13127888/20190809_zymo_RUfull.tar.gz
tar -xvzf $CUSTOM_APPS/data/uncalled_data/20190809_zymo_RUfull.tar.gz
# Ctrl-C after a few files + delete incomplete file

# get readuntil UNCALLED files (the ones they link to in their README)
wget https://labshare.cshl.edu/shares/schatzlab/www-data/UNCALLED/simulator_files/20190809_zymo_seqsum.txt.gz
gunzip --keep 20190809_zymo_seqsum.txt.gz
wget https://labshare.cshl.edu/shares/schatzlab/www-data/UNCALLED/simulator_files/20190809_zymo_uncalled.paf.gz
gunzip --keep 20190809_zymo_uncalled.paf.gz

# get Zymomock ref index from https://files.zymoresearch.com/protocols/_d6322_zymobiomics_hmw_dna_standard.pdf
wget https://s3.amazonaws.com/zymo-files/BioPool/D6322.refseq.zip
unzip D6322.refseq.zip

# basecall raw reads with guppy to get a sequencing summary
srun --job-name interactive --cpus-per-task 6 --mem 8G --partition gpu --gres=gpu:1 --time 04:00:00 --pty /bin/zsh
guppy_basecaller --disable_pings --verbose_logs --input_path 20190809_zymo_RUfull/fast5 --save_path basecalled_reads_fast -c dna_r9.4.1_450bps_fast.cfg -x 'cuda:0'
guppy_basecaller --disable_pings --verbose_logs --input_path 20190809_zymo_RUfull/fast5 --save_path basecalled_reads_hac -c dna_r9.4.1_450bps_hac.cfg -x 'cuda:0'

# build with BWA for so it does not need to be rebuilt
mamba activate ont_project_310
uncalled index -o D6322_yeast D6322.refseq/Genomes/Saccharomyces_cerevisiae_draft_genome.fasta

#Format: uncalled sim E.coli.fasta /path/to/control/fast5s --ctl-seqsum /path/to/control/sequencing_summary.txt --unc-seqsum /path/to/uncalled/sequencing_summary.txt --unc-paf /path/to/uncalled/uncalled_out.paf -t 16 --enrich -c 3 --sim-speed 0.25 > uncalled_out.paf 2> uncalled_err.txt
uncalled sim D6322_yeast 20190809_zymo_RUfull/fast5 --ctl-seqsum basecalled_reads_fast/sequencing_summary.txt --unc-seqsum 20190809_zymo_seqsum.txt --unc-paf 20190809_zymo_uncalled.paf -t 16 --enrich -c 3 --sim-speed 0.25 > sim_uncalled_out.paf 2> sim_uncalled_err.txt

# this gives the following segmentation fault!
# [1]    2754126 segmentation fault (core dumped)  uncalled sim D6322_yeast 20190809_zymo_RUfull/fast5 --ctl-seqsum  --unc-seqsu
```

No precise error message is given and people have reported similar issues on UNCALLED's github.