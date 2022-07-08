# CaFe Online Replay
This repository, `cafe_online_replay`, is a version of the official `hallc_replay` repository with additions specific to the CaFe experiment. This online version has replay scripts for data online monitoring of the data, as well as specialized shell scripts for analysis of special runs such as target_boiling, H(e,e'p) elastics, production, detector calibrations, etc. 

Please refer to the following useful link:
[https://hallcweb.jlab.org/wiki/index.php/CaFe\_Experiment](https://hallcweb.jlab.org/wiki/index.php/CaFe_Experiment)

where you can find helpful information related to the experiment. 
# Instructions for Shift Takers
 NOTE: In the instructions below, the pound sign, `# `, denotes comments, 
 and the `$` denotes the user terminal command-line input

### General Instructions:
Go to the main `cafe_online_replay` directory on the cdaq cluster: <br>

`step 1:` Assuming you are on a desktop [hcdesk] machine in the Counting House, open a terminal: 

```sh
# login to cdaql[1,2 3, 5 or 6]
$ ssh cdaql3 
$ cd cafe-2022/cafe_online_replay 
```

This is the top-level directory where the relevant scripts for online data replay and output directories are symbolically linked. 


### HMS/SHMS 50k Replay:
To monitor the detector response in both spectrometers, a sample of 50,000 event
replay is done on the HMS and SHMS separately. This monitoring is done on the on-going live run and its purpose is to access the quality of the detectors and quickly determine whethere there are any issues that might need to be resolved. 
The 50k replay can be done simultaneosly, provided the user has multiple tabs
opend on the terminal.

```sh
# execute 50k HMS replay script
$ ./run_coin_hms.sh

# execute 50k SHMS replay script
$ ./run_coin_shms.sh
```

The script will automatically replay the most recent on-going run for up to 50,0000 events. However, if the user wants to replay a different run with a particular number of events:

```sh
# execute 50k HMS replay script for specific run
$ ./run_coin_hms.sh <run_number> <event_number>

# execute 50k SHMS replay script for specific run
$ ./run_coin_shms.sh <run_number> <event_number> 
```

There is an additional option which is used to set a run as a **golden run**, and that can simply be done via:

```sh
# execute 50k HMS replay script for specific run
$ ./run_coin_hms.sh <run_number> <event_number> set-golden

# execute 50k SHMS replay script for specific run
$ ./run_coin_shms.sh <run_number> <event_number> set-golden

# numerical example: ./run_coin_shms 3288 75000 set-golden
```
If the `set-golden` option is used, then it will make a copy of the output ROOTfile and rename it with the '_golden.root' suffix in the name. The first production run will need to be set as the **golden run** so that subsequent runs can be compared to the golden one. <br>

The purpose of the **golden run** is that it be used as a benchmark to compare the detector pedestal values of the golden run with the present run. That way, a monitoring of the pedestals can be done for each detector, and one can assess the stability of the pedestals.

After 	50k replay event finishes, an onlineGUI will appear with the monitorin plots, and once the user cycles through all the histograms (and click Exit GUI in red square) the data will be ouputted to the following directories:

* `ROOTfiles/(s)hms50k/`
* `REPORT_OUTPUT/(s)hms50k/`
* `HISTOGRAMS/(s)hms50k/PDF/`

### CaFe Sample/Production Replay:
If the 50k replay monitoring plots are OK, the user can proceed to executing a sample physics analysis of the data as follows:

```sh
# execute cafe sample replay script
./run_cafe_sample.sh <run_number> <kin_type> <event_number>
```  
The required arguments are the `<run_number>` and `<kin_type>`. If no event number is specified, it will default to 100k events. The `<kin_type>` specifies the type of study to be done and can be one of the following: <br>

* `<bcm_calib>`
* `<lumi>`
* `<optics>`
* `<heep_singles>`
* `<heep_coin>` 
* `<MF>`
* `<SRC>`

Each of these `<kin_type>` arguments triggers very specific analysis to be done (hidden from the user)