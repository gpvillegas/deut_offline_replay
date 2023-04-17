# Overview
Useful Docuentation: [https://hallcweb.jlab.org/DocDB/0010/001032/001/analysis_notes.pdf]() <br>
Analysis E-LOG: [https://hallcweb.jlab.org/elogs/Deuteron+Analysis/]()
### directory structure: <br>

`scripts/set_reftimes.C` main script to check/set reference times or detector time window cuts

`scripts/set_reftimes.h` header file to define variables used in the main script

`set_reftime.sh`: shell script that (1) replays data and (2) checks ref. time of replayed data

`set_timewin.sh`: shell script that (1) replays data and (2) checks detector time windows of replayed data

**NOTE**: Since checking/setting ref. times/detector time windows requires specific leaf variables to be in a ROOTfile, the shell script automatically replays the data with the  required variables and thne analyzes the output ROOTfile to determine the existing reference time/detector time windows currently
set and in the case of hodoscopes/calorimeters, it suggests new time window cuts as well, since these are high-density detectors with many channels for a user to set manually.

## how-to guide:

### executing the scripts
To replay + analyze data in a single step
(e.g. running HCANA on a raw .dat file to generate a .root file + analyze the .root file to generate the relevant sub-directories with analysis plots, execute the following:

```sh
# for setting / checking the reference time cuts
./set_reftime.sh <run_number> <evt_number> replay

# NOTE: user need to check if ref. times are set properly, as it
# may be necessary to update the reference time .param before
# running the script to set time window cuts

# for setting / checking detector time window cuts
/set_timewin.sh <run_number> <evt_number> replay
```

If the .root file already exists, and user would like to re-analyze (e.g. with different ranges in cuts, histograms, etc.) then execute the following:

```sh
# for setting / checking the reference time cuts
./set_reftime.sh <run_number> <evt_number> 

# NOTE: user need to check if ref. times are set properly, as it
# may be necessary to update the reference time .param before
# running the script to set time window cuts

# for setting / checking detector time window cuts
/set_timewin.sh <run_number> <evt_number>
```

### analysis outputs
1) When the user executes `./set_reftimes.sh`,  the following output directories will be generated: <br>

`Time_cuts_refTime<run_number>/` with the contents: <br>

`HMS/`: sub-directory with plots showing the existing HMS reference times <br>
`SHMS/`: sub-directory with the plots showing the existing SHMS reference times <br>
`.root` a ROOTfile with the relevant histogram objects that were used to make the plots.

At this point, the user can decide whether the existing reference time cuts (blue lines) are satisfactory or whether they will need to be modified before executing the `./set_timewin.sh` script


2)  When the user executes `./set_timewin.sh`,  the following output directories will be generated: <br>

`Time_cuts_tWinSet<run_number>/` with contents: <br>

`param_files`: sub-directory with both `_existing.param` or `_new.param` files <br>
`HMS/`: sub-directory with detector plots showing the existing/new HMS time window cuts <br>
`SHMS/`: sub-directory with detector plots showing the existing/new SHMS time window cuts <br>
`.root` a ROOTfile with the relevant histogram objects that were used to make the plots.

There may also be an existing `.pdf` file with plots of trigger times, but this is subject to which variables the user will set a cut on the `tcoin.param` file (and can always be changed in the code, depending on which variables from `tcoin.param` the user will examine)

