analyse_h5
============

Analysis tools for h5 files
- writing out and plotting with root

# Initialization
Depends on..
- python3
- and root 

Environment is set to ```source /opt/exp_software/limadou/set_env_standalone.sh```

Clone the git repository:
~~~
git clone https://github.com/CoralieNeubueser/analyse_h5.git .
~~~

# Directory structure
--working_directory  <br />
----|--python  <br />
--------|  <br />
--------|--run.py  <br /> 
--------|--utils.py  <br />
--------|--readH5.py  <br />
----|--root  <br />
--------|  <br />
----|--plots  <br />
--------|  <br />

# Run analysis
use run.py scipt in python/

~~~
python3 python/run.py --numRuns 1 --hepd --test
~~~

- it finds h5 files in ```/storage/gpfs_data/limadou/data/flight_data/L3_test``` if flag ```--hepd``` is raised. To run on LEOS data, specify ```--hepp``` and the files are found in ```/storage/gpfs_data/limadou/data/cses_data/HEPP_LEOS/```. 

- the analysis is run from the top of the list, in case that the corresponding root file already exists in ```{workDir}/root/```, the analysis is run on the next file on the list. 

- add ```--merge``` to merge the output root trees to one.
- add ```--ana``` to get energy bin vise flux distributions and L-pitch maps

# Examples
to work on the root trees and run higher level analysis, try e.g.:
~~~
python3 python/analyseRoot.py --inputFile ...
~~~
this example analysis divides the tree into energy bins, and fills histograms. It also implements a cut on the SAA and the correction of the flux for different geometrical factor, see [talk on Ryver](https://contattafiles.s3.us-west-1.amazonaws.com/tnt11944/xPI899xCKI3S9AZ/20200513_LIMADOU_analysis_L3_status.pdf)


