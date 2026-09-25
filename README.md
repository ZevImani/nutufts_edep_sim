# NuTufts Edep-Sim Data set Maker

This is a repository of code to make toy datasets for ML/AI experiments for LArTPCs.

See below for example commands to complete the initial setup and build the programs you need.

Also see the folder `troubleshooting` for a record of tricky errors that we've run into.

## Initial Setup and Making a first test dataset

```
# go into a worker node 
# i dont expect this to matter too much, but maybe there are 'rules' imposed about talking/writing to databases if you're on the login noe?
srun --pty -p wongjiradlab --time 8:00:00 bash

# setup singularity
module load singularity/3.5.3

# go to your copy of the repo. (replace relevant info in between the [] brackets.)
cd /cluster/tufts/wongjiradlabnu/[your-username-here]/[path-of-folders]/nutufts_edepsim_dataset_maker/

# if you want to get a copy of the recpo
cd /cluster/tufts/wongjiradlabnu/[your-username-here]/
mkdir [path-of-folders]
git clone https://github.com/NuTufts/nutufts_edepsim_dataset_maker.git
# this repo is currently private, so ask Taritree to add you to NuTufts
# though presumably you couldnt read this readme without already being addded

# start the container
singularity shell -B /tmp:/tmp,/cluster:/cluster /cluster/tufts/wongjiradlabnu/larbys/larbys-container/simpledet_petastorm_u20.04_geant4.10.06.p03.sif bash

# start a bash shell inside the container
source setenv_py3.sh

# if you were starting from a fresh container, load the submodules and then build
# skippable if you already have built everything
git submodule init
git submodule update
# if you pulled and updated the repo, to make sure the version of the submodules are current, run 'git submodule update'
# build
source build_edep_sim.sh
source build_simpledet.sh

# go into simpledet/test where some scripts live
#  if you need to make an example set of events from the simulation
# will simulate 5 electrons by default
edep-sim -g ../../shower_classifier_clue/larblock.gdml -o test.root -e 5 ../../shower_classifier_clue/electron.mac

# convert simulation information into images+metadata and store into petastorm db
python3 test_save2petastormdb.py --input-edepsim test.root -ow --petastorm-db-folder /tmp/test_db/ --tag electron_test --pdgcode 11 --runid 0

# Hack: visulaize -- to avoid petastorm db 
python3 test_save2petastormdb.py -v true --input-edepsim test.root -ow --petastorm-db-folder /tmp/test_db/ --tag electron_test --pdgcode 11 --runid 0

# the above command made a pyspark database that now lives in /tmp/test_db
# look inside the folder
ls /tmp/test_db/
ls /tmp/test_db/
 _SUCCESS   _common_metadata  'partition=pdg11_run0000_electron_test'

# the entries you made are in the partition older: 'partition=pdg11_run0000_electron_test'
# the '-ow' argument in principle should overwrite this partition if you use the same '--tag' '--pdgcode' '--runid' combination. 
# the partition is jsut a way to have flexibility to make events and organize them. 
# not 100% sure, but i've just removed partition folders without any issues.
# also, you can image that providing a 'train' and 'valid' label allows one to load certain partitions during training/testing

# let's test that the entries in the database are good
python3 test_dataloader.py --petastorm-db-folder /tmp/test_db
```