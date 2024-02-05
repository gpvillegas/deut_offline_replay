This directory is used to check calibrations for the
various detectors

SCRIPTS                     
---------
--> replay_deut_checkCalib.sh

description: This shell script replays the raw (.dat) files 
       and outputs the raw (.root) files for calibration
       purposes. The leaf variables to be written are 
       specified under: DEF-files/deut_checkCalib.def 
       and may be modified.


--> checkCalib.C (checkCalib.h)

description: This script is used to make calibration plots for each
             detector of each spectrometer of a replayed data file
	     to check if its calibrated properly.


---> replay_all.sh

description: This shell script facilitates the user's life
             by taking a run list to check the calibrations
	     as input, and running the two scripts above to
	     1) generate a raw .root file and 2) analyze that
	     .root file to produce the relevant histograms
