#!/usr/bin/env bash 

# it seems 2 CPUs are fine based on condor log average resource usage, simforward,muxscan,stopthread
##CONDOR request_cpus=4
# takes about 8GB of memory
##CONDOR request_memory=64000
##CONDOR request_disk=100G
##CONDOR +JobBatchName = "ont_enrich_usecase"
##CONDOR log = /home/mmordig/joblogs/job-$(ClusterId)-$(ProcId).log
##CONDOR output = /home/mmordig/joblogs/job-$(ClusterId)-$(ProcId).out
##CONDOR error = /home/mmordig/joblogs/job-$(ClusterId)-$(ProcId).err

#SBATCH --job-name=enrich_usecase-%j
#SBATCH --error=/cluster/home/mmordig/joblogs/job-%j.err
#SBATCH --output=/cluster/home/mmordig/joblogs/job-%j.out
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
## #SBATCH --time=00:10:00
#SBATCH --time=2:00:00
# not avail #SBATCH --tmp=10G
#SBATCH --partition=compute


# launch_condor_job 20 --- ~/ont_project_all/ont_project/usecases/enrich_usecase_submission.sh
# sbatch ~/ont_project_all/ont_project/usecases/enrich_usecase_submission.sh

echo "Content of job ad file $_CONDOR_JOB_AD:"; cat "$_CONDOR_JOB_AD"
echo "Starting job with args: " "$@"
echo "Cwd: $(pwd)"
source ~/.bashrc

cd ~/ont_project_all/ont_project/

source ~/ont_project_all/ont_project_venv/bin/activate
set -ex
export PATH=~/ont_project_all/tools/bin:$PATH && which minimap2

output_dir=${1:-full_genome_run_sampler_per_window}
config_rel_dir=${2:-sampler_per_window}

cd runs/enrich_usecase
rm -rf "$output_dir"
mkdir -p "$output_dir"
ln -s "$(pwd)"/data "$output_dir"
cp -rL configs/full_genome_run/"${config_rel_dir}" "$output_dir"/configs # -L: expand symlinks
cd "$output_dir"
pwd
python ~/ont_project_all/ont_project/usecases/enrich_usecase.py

# symlink job output files to directory
# parse stderr from job ad file
# Err = "/home/mmordig/joblogs/job-13951206-0.err"
if [ -n "$_CONDOR_JOB_AD" ]; then
    grep -oP '(?<=Err = ").*(?=")' "$_CONDOR_JOB_AD" | xargs -I {} cp {} .
    grep -oP '(?<=Out = ").*(?=")' "$_CONDOR_JOB_AD" | xargs -I {} cp {} .
    grep -oP '(?<=Log = ").*(?=")' "$_CONDOR_JOB_AD" | xargs -I {} cp {} .
    # grep -oP '(?<=Err = ").*(?=")' "$_CONDOR_JOB_AD" | xargs -I {} ln -s {} .
    # grep -oP '(?<=Out = ").*(?=")' "$_CONDOR_JOB_AD" | xargs -I {} ln -s {} .
    # grep -oP '(?<=Log = ").*(?=")' "$_CONDOR_JOB_AD" | xargs -I {} ln -s {} .
fi

echo "Done with job, pwd $(pwd)"
