#!/usr/bin/env bash 

##CONDOR request_cpus=1
# takes very little memory, 100 MB
##CONDOR request_memory=16G
##CONDOR request_disk=20G
##CONDOR log = /home/mmordig/joblogs/job-$(ClusterId)-$(ProcId).log
##CONDOR output = /home/mmordig/joblogs/job-$(ClusterId)-$(ProcId).out
##CONDOR error = /home/mmordig/joblogs/job-$(ClusterId)-$(ProcId).err

# <run_dir>: the directory containing the seqsum to plot
# launch_condor_job 15 --- ~/ont_project_all/ont_project/usecases/plot_existing_seqsum_submission.sh <run_dir>

echo "Content of job ad file $_CONDOR_JOB_AD:"; cat "$_CONDOR_JOB_AD"
echo "Starting job with args: " "$@"
echo "Cwd: $(pwd)"
source ~/.bashrc

source ~/ont_project_all/ont_project_venv/bin/activate
export PATH=~/ont_project_all/tools/bin:$PATH && which minimap2

set -ex
cd "$1" # if empty, it is cwd

python ~/ont_project_all/ont_project/usecases/plot_existing_seqsum.py

echo "Done with job, pwd $(pwd)"
