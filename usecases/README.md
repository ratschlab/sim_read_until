# Usecases

Config files are located in the directory `configs` which contain configs for the simulator (and ReadFish).
The paths in these scripts may need to be adapted to your local setup.

Activate the environment and start jupyter lab:
```{bash}
source ~/ont_project_all/ont_project_venv/bin/activate
cd ~/ont_project_all/ont_project
jupyter lab
```

## Setup

### TLDR

You can get the data and install the dependencies with:
```{bash}
cd ~/ont_project_all/ont_project
curl https://public.bmi.inf.ethz.ch/user/mmordig/ont_project/prepackaged/usecase_data.tar.gz -O
tar -xvzf usecase_data.tar.gz
cd runs
mkdir enrich_usecase && cd enrich_usecase
ln -s ~/ont_project_all/ont_project/usecases/configs/enrich_usecase configs
ln -s ../data data
cd ..
mkdir run_replication && cd run_replication
ln -s ~/ont_project_all/ont_project/usecases/configs/run_replication configs
ln -s ../data data
cd ..
cd ..

# install ReadFish
git submodule update --init --depth 1 external/ont_readfish
source ~/ont_project_all/ont_project_venv/bin/activate
pip install -e './[readfish]'

# optional: install NanoSim and minimap2
# git submodule update --init --depth 1 external/ont_nanosim
# bash usecases/install_usecase_deps.sh

# adapt config in ~/ont_project_all/ont_project/usecases/configs/ to fit your setup
~/ont_project_all/ont_project/usecases/replicate_run_submission.sh "sampler_per_rolling_window_channel"
~/ont_project_all/ont_project/usecases/enrich_usecase_submission.sh
```

You may want to check out the Docker image as described in the `README.md` in the repo root, if it does not work.

### Explanation

We provide more details:

```{bash}
cd ~/ont_project_all/ont_project
mkdir runs && cd runs
mkdir data && cd data
curl https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz -O
wget https://labshare.cshl.edu/shares/schatzlab/www-data/UNCALLED/simulator_files/20190809_zymo_seqsum.txt.gz
# gz file slow to read and write, so we uncompress it
gunzip 20190809_zymo_seqsum.txt.gz

source ~/ont_project_all/ont_project_venv/bin/activate
normalize_fasta "chm13v2.0.fa.gz" "chm13v2.0_normalized3chroms.fa.gz" --seq_names "chr20,chr21,chr22"
normalize_fasta "chm13v2.0.fa.gz" "chm13v2.0_normalized.fa.gz"
cd ..

mkdir enrich_usecase && cd enrich_usecase
ln -s ~/ont_project_all/ont_project/usecases/configs/enrich_usecase configs
ln -s ../data data
cd ..

mkdir run_replication && cd run_replication
ln -s ~/ont_project_all/ont_project/usecases/configs/run_replication configs
ln -s ../data data
```

If you want to run the usecase, install the dependencies with:
```{bash}
source ~/ont_project_all/ont_project_venv/bin/activate
pip install -e './[readfish]'
git submodule update --init --depth 1 external/ont_readfish
```

Optionally, to install minimap2 and NanoSim as well:
```{bash}
# you need to have conda or mamba installed
git submodule update --init --depth 1 external/ont_nanosim
bash usecases/install_usecase_deps.sh
```

## Parameter Extraction from an Existing Run

This experiment learns parameters from an existing run and then simulated a run with the learned parameters.
The gaps are learnt from the run.
Since the input to the simulator is in terms of basecalled reads, reads consisting of random letters are generated to match the read durations during the real run as closely as possible.

Files:
- `replicate_run.py`: end-to-end script that runs the parameter extraction and simulation, see the instructions in that file
- `replicate_run_submission.sh`: condor submission script, can also be run locally
- `compare_replication_methods.ipynb`: merge plots into one (for the paper) to compare the different parameter extraction methods

Parameter extraction takes some time (few minutes) since the sequencing summary file must be loaded. Therefore, we implement a caching mechanism in the function `create_simparams_if_inexistent`:
```{bash}
sequencing summary from an existing run 
--> sequencing summary with mux scans removed (prefixed with 'no_mux_scans_') 
--> cleaned sequencing summary (prefixed with 'cleaned_'): some reads of a channel overlap, so we shift them
--> sim params (per channel, saved in one '.npz' file)
```
We check from the back whether any of the files already exists and use the existing file if possible.
Make sure to delete the appropriate files when changing relevant parameters.
When running several configurations in parallel and some cache files do not exist, make sure they don't get created by different processes at the same time.

## Enrichment with ReadFish

This combines `ReadFish` with the `SimReadUntil` simulator. Reads are generated on the fly by sampling read start position, length and strand. With this ground-truth information, the alignment step in ReadFish can be accelerated which is useful when we accelerate the run by a factor of 10 where minimap2 alignment may otherwise become a bottleneck.
Alternatively, a reads file can be provided by adding `reads_file = <path>` to the simulator config.
If the read ids are NanoSim ids with ground-truth alignment information, `minimap2` is not needed. Alternatively, a minimap2 reference index can be provided by adding it in the ReadFish config file. The index itself can be created by uncommenting the line `create_minimap_index_if_inexistent` in the usecase file.

Files:
- `enrich_usecase.py`: end-to-end script that runs an enrichment with ReadFish connected to the simulator, see the instructions in that file
- `enrich_usecase_submission.sh`: condor submission script, can also be run locally
- `install_usecase_deps.sh`: to install `minimap2` and `NanoSim` (optional), launch it from the repo root
- `create_nanosim_reads.ipynb`: notebook to create NanoSim reads that can be fed into the simulator by modifying the config file

## Other Files

These files are for our own reference and may not work for you out of the box:
- `analyze_readfish_outputs.py`: to check whether ReadFish is mapping reads correctly by parsing the ground-truth from the read id
- `plot_existing_seqsum.py`: to plot an existing sequencing summary file, e.g., from a real run
- `remove_mux_scans.ipynb`: notebook showing how mux scans are removed (you don't need to run this, this is done automatically in the usecases)
- `prepare_small_refgenome.py`: to create a small reference genome for the usecase
- `results_preparation.md`: commands to create the results in the paper
- `plot_existing_seqsum_submission.sh`: condor submission script, can also be run locally
