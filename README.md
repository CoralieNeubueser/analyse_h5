analyse_h5
============

Analysis tools for h5 files
- writing out and plotting with root

# Initialization
Depends on..
- python3
- and root 

Environment is set to 
```source /opt/exp_software/limadou/set_env_standalone.sh```

Clone the git repository:
~~~
git clone https://github.com/CoralieNeubueser/analyse_h5.git 
~~~

# Directory structure
--working_directory  <br />
----|--python  <br />
--------|  <br />
--------|--run.py  <br /> 
--------|--utils.py  <br />
--------|--readH5.py  <br />
----|--data  <br />
--------|  <br />
--------|--root <br />
--------|--averages <br />
--------|--bad_runs <br />
--------|--pdf <br />
----|--plots  <br />
--------|  <br />

# Read h5 files, and produce ouput root files
use run.py script in python/ directory:

~~~
python3 python/run.py --numRuns 1 --hepd/hepp
~~~

- it finds h5 files in ```/storage/gpfs_data/limadou/data/flight_data/L3_test``` if flag ```--hepd``` is raised. To run on HEPP-L data, specify ```--hepp``` and the files are found in ```/home/LIMADOU/cneubueser/public/HEPP_august_2018/```. 

- the analysis is run from the top of the list, in case that the corresponding root file already exists in ```data/root/``` ```{limadou_directory}=/storage/gpfs_data/limadou/data/flight_data/analysis/```, the analysis is run on the next file on the list. 

OPTIONS:
- ```-q```: run within minimal printed statements
- ```--test```: in combination with ```--hepd``` run on data at the 1.-5. August 2018
- ```--submit```: will submit each single job via condor to the batch farm on cnaf

# Determine the daily average fluxes
This is a locally run analysis.
1. first merge the root files within a certain time range, either use ```hadd -f -k``` or add ```--merge``` to the command, e.g. for the test period of 1.-5. August 2018:
~~~
python3 python/run.py --hepd --test --merge
~~~
- the output of this example will be stored locally in your working directory: ```{limadou_directory}/data/root/hepd/``` as ```L3_test/L3h5_orig/all.root```

2. run the analysis on the file with all half-orbits, here e.g. for 1.-5. day of August 2018:
~~~
python3 python/writeDayAverages.py --data hepd --inputFile {limadou_directory}/data/root/hepd/L3_test/L3h5_orig/all.root 
~~~
The input file can be varied, and the type of data (either HEPD of HEPP-L) needs to be specified.
- the mean and rms of the flux distributions get stored for equatorial pitch angle / L-shell value and every energy bin, only if entries of the histograms are larger than a threshold (default=100). The threshold can be changed to any integer x, using the flag ```--thr x```. 
The flux averages will be stored in ```{limadou_directory}/data/averages/hepd/{yyyymmdd}_min_{thr}ev.txt```

# Create root file of potential particle bursts
a new root file with fluxes larger than the mean+(3x rms) can be created in the next step, with:
~~~
python3 python/findHighFluxes.py --data hepd --inputFile {limadou_directory}/data/root/hepd/L3_test/L3h5_orig/all.root
~~~
The input file can be varied, and the type of data (either HEPD of HEPP-L) needs to be specified.
You will find the output in the directory of the ```inputFile``` as ```all_highFluxes.root```

# Run the Clustering algorithm on high fluxes for EQ correlation
The previously created 'high flux' data is used as an input for the clustering algorithm that consists of 3 steps:
1. determine (again) the highest fluxes by a cut to be used as threshold for the seeds of the clustering. Here you have the choice between '99perc','z_score_more2','z_score_more3','dummy_cut' (```--cut str```). The definitions can be found [here](https://agenda.infn.it/event/27436/contributions/139005/attachments/82451/108325/clustering_20210614.pdf). These cuts are determined per day, per L shell (in int steps), and per alpha bin. The values are stored in pickle files.
2.  If you want to run the clustering per month, in the second step the accessible data over that defined month (```--month YYYYMM```) are merged into a sinlge root file.
3.  The actual clustering is performed, here there are 2 possible paramters to be chosen: window size in seconds (```--window x```) and the number of seeds that have to cross the threshold to be found wihtin that window in order to build a cluster (```--seeds y```).

An example execution in the command line looks like that:
~~~
python3 python/run.py --cluster --hepp_l --channel wide --day 20180801 --numRuns 1
~~~
The output of step 1. is stored in: ```{limadou_directory}/data/thresholds/{SW_version}/{det}/```

The root file with the fluxes, and the cluster info is stored in: ```{limadou_directory}/data/root/{SW_version}/{det}/clustered/{type_of_cut}/{windowSize}s_window/{seeds}_seeds/```

# Example analysis
to work on the root trees and run higher level analysis, try e.g.:
~~~
python3 python/analyseRoot.py --inputFile ...
~~~
this example analysis divides the tree into energy bins, and fills histograms. It also implements a cut on the SAA and the correction of the flux for different geometrical factor, see [talk on Ryver](https://contattafiles.s3.us-west-1.amazonaws.com/tnt11944/xPI899xCKI3S9AZ/20200513_LIMADOU_analysis_L3_status.pdf)


