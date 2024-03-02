#!/usr/bin/env bash 

# note: aligned and unaligned reads are not shuffled here because they are shuffled by the ReadPoolFromFile
# more efficient to run with 2 cpus and run many in parallel

##CONDOR request_cpus=2
##CONDOR request_memory=6000
##CONDOR request_disk=2G
##CONDOR +JobBatchName = "ont_enrich_usecase"
##CONDOR log = /home/mmordig/joblogs/job-$(ClusterId)-$(ProcId).log
##CONDOR output = /home/mmordig/joblogs/job-$(ClusterId)-$(ProcId).out
##CONDOR error = /home/mmordig/joblogs/job-$(ClusterId)-$(ProcId).err
# seems to be broken/filesystem very slow
##CONDOR Requirements = (Machine != "g110.internal.cluster.is.localnet")

# bash ~/ont_project_all/ont_project/usecases/generate_nanosim_reads.sh 1000 3
# conda activate nanosim
# sbatch ~/ont_project_all/ont_project/usecases/generate_nanosim_reads.sh 100000 3
# or use ##CONDOR queue 100 together with $(Item)
# for seed in range(1, 100+1):
#     print(f"sbatch ~/ont_project_all/ont_project/usecases/generate_nanosim_reads.sh 100000 {seed}")
#
# initially
# mkdir -p runs/nanosim_models
# echo tar -xvzf external/ont_nanosim/pre-trained_models/human_NA12878_DNA_FAB49712_guppy.tar.gz -C "runs/nanosim_models"

# set -x
# source ~/.bashrc
# conda hell: conda not found, not sure why, so hardcoding python executable from conda env

# for python, samtools; conda activate not really working (need to copy entire env)
export PATH=/home/mmordig/tools/mambaforge/envs/nanosim/bin:$PATH
# export PATH=/home/mmordig/miniforge3/envs/nanosim/bin:$PATH

set -eux

cd ~/ont_project_all/ont_project/

num_reads=$1
seed=$2

# output_dir=runs/data/nanosim_reads/human_genome_med15000
output_dir=runs/data/nanosim_reads/human_genome_med15000_alignedrate2
# output_dir=runs/data/nanosim_reads/human_genome
# output_dir=runs/data/nanosim_reads/human_genome_with_flanking
mkdir -p "$output_dir"
num_procs=$(nproc)
# ((num_procs--)) # 1 manager process, not really needed
# num_procs=1 #todo

# genome=runs/data/random_genome.fasta # see below for how to generate
genome="runs/enrich_usecase/data/chm13v2.0_normalized.fa.gz"
# aligned_rate="100%"
aligned_rate="2"

echo "nanosim read generation: generating ${num_reads} using seed $seed using $num_procs threads from genome '$genome' with aligned_rate '$aligned_rate' into output_dir '$output_dir'"

echo "Generating reads"
# conda slow to run, instead use "conda activate nanosim" once and then launch the script several times
# conda run -n nanosim python \
rm "$output_dir/reads_seed$seed"* || true
# in NanoSim, replaced by uniform distribution now because median length was unreliable, sd=6.9=ln(1000) is the std of the lognormal (lognormal = distribution whose log is normally distributed with stddeviation std)
python \
    "external/ont_nanosim/src/simulator.py" genome \
    --model_prefix "runs/nanosim_models/human_NA12878_DNA_FAB49712_guppy/training" \
    --ref_g "$genome" \
    -dna_type linear \
    -med 15000  -max 20000 -min 400 -sd 2000 \
    --output "$output_dir/reads_seed$seed" \
    --number "${num_reads}" \
    --seed "$seed" \
    --strandness 0.5 \
    --basecaller guppy \
    --aligned_rate "$aligned_rate" \
    --num_threads "$num_procs" \
    --no_flanking \
    --no_error_profile

echo "Merging files"
# merge 'reads_seed1_aligned_reads.fasta', 'reads_seed1_unaligned_reads.fasta' into 'reads_seed1_merged_reads.fasta'
files=$(ls "$output_dir/reads_seed${seed}_"*)
merged_file="$output_dir/reads_seed${seed}_merged_reads.fasta"
# note: aligned and unaligned reads are not shuffled here because they are shuffled by the ReadPoolFromFile
cat $files > "$merged_file"
rm $files

echo "Generating .fai file"
samtools faidx "$merged_file"

echo "nanosim read generation: done with seed $seed"


# # generate small fake genome and write to file using pysam
# from simreaduntil.shared_utils.dna import get_random_DNA_seq
# from Bio import SeqIO
# from Bio.Seq import Seq
# with open("runs/data/random_genome.fasta", "w") as fasta:
#     SeqIO.write((SeqIO.SeqRecord(id=f"fakechr_{i}", seq=Seq(get_random_DNA_seq(1_000_000))) for i in range(20)), fasta, "fasta")