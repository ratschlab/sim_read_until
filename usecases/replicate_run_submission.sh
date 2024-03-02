#!/usr/bin/env bash 

##CONDOR request_cpus=2
# takes about 8GB of memory
##CONDOR request_memory=32G
##CONDOR request_disk=100G
##CONDOR +JobBatchName = "ont_replicate_run"
##CONDOR log = /home/mmordig/joblogs/job-$(ClusterId)-$(ProcId).log
##CONDOR output = /home/mmordig/joblogs/job-$(ClusterId)-$(ProcId).out
##CONDOR error = /home/mmordig/joblogs/job-$(ClusterId)-$(ProcId).err

# launch_condor_job 20 --- ~/ont_project_all/ont_project/usecases/replicate_run_submission.sh "constant_gaps"
# launch_condor_job 20 --- ~/ont_project_all/ont_project/usecases/replicate_run_submission.sh "replication"
# launch_condor_job 20 --- ~/ont_project_all/ont_project/usecases/replicate_run_submission.sh "sampler_per_window"
# launch_condor_job 20 --- ~/ont_project_all/ont_project/usecases/replicate_run_submission.sh "sampler_per_rolling_window_channel"

# check if environment var exists
if [ -n "$_CONDOR_JOB_AD" ]; then
    echo "Content of job ad file $_CONDOR_JOB_AD:"; cat "$_CONDOR_JOB_AD"
fi

echo "Starting job with args: " "$@"
echo "Cwd: $(pwd)"
source ~/.bashrc

cd ~/ont_project_all/ont_project/

source ~/ont_project_all/ont_project_venv/bin/activate
set -ex
export PATH=~/ont_project_all/tools/bin:$PATH && which minimap2

cd runs/run_replication

method=$1

# for method in "constant_gaps" "replication" "sampler_per_window" "sampler_per_rolling_window_channel"; do
echo "Method: $method"

rm -rf "$method" && mkdir "$method"
cd "$method"
ln -s ../data .
ln -s ../configs/"$method" configs
python ~/ont_project_all/ont_project/usecases/replicate_run.py
# profile_filename=scalene_$(date +%s).html; echo "Filename: $profile_filename"; scalene --off --html --program-path ~/ont_project_all/ont_project --outfile "$profile_filename" --profile-all ~/ont_project_all/ont_project/usecases/replicate_run.py
cd ..
# done

echo "Done with job, pwd $(pwd)"
