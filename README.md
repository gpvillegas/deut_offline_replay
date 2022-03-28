# CaFe Online Replay
This repository, `cafe_online_replay`, is a version of the official `hallc_replay` repository with additions specific to the CaFe experiment. This online version has replay scripts for data online monitoring of the data, as well as specialized shell scripts for analysis of special runs such as target_boiling, H(e,e'p) elastics, production, detector calibrations, etc. 

Please refer to the following useful link:
[https://hallcweb.jlab.org/wiki/index.php/CaFe_Experiment](https://hallcweb.jlab.org/wiki/index.php/CaFe_Experiment)

where you can find helpful information related to the experiment. 
# How-To Guide for Starters
This guide is intended for *novice level* on hall c analysis, using the CaFe experiment as an example.

### Exercise #1:
Set up and get familiar with the general hall c analysis replay structure on ifarm. <br>

`step 1:` Assuming you have a Jefferson Lab account, log-in to ifarm: 
>ssh -Y *user*@login.jlab.org <br>
>ssh -Y ifarm <br>
>source /site/12gev_phys/softenv.csh 2.5 # setup necessary environment variables

`step 2:` Go to the relevant work directory and setup the relevant Hall C Analyzer repository.
>  # Change directories. For more info see [https://hallcweb.jlab.org/wiki/index.php/CaFe_Disk_Space](https://hallcweb.jlab.org/wiki/index.php/CaFe_Disk_Space) <br>
>cd /w/hallc-scshelf2102/c-cafe-2022 <br>
>
>  # If you don't have a directory, make one and change to it !<br>
> mkdir *user* ; cd *user*<br>
> 
>  # Clone and properly setup the Hall C analyzer source code <br>
>     >> git clone https://github.com/JeffersonLab/hcana <br>
>     >> cd hcana <br> 
>    # 3 separate commands (usually only required 1st time afte cloning hcana, but if there are updates to podd, then it need to be done.) <br>
>     >> git submodule init; git submodule  sync; git submodule update <br>
>     >> source setup.csh  # setup hcana environment variables <br> 
>   # the command below will compile the analyzer and generate an executable, *hcana*<br>
>     >> scons -jN # scons is a compiler and -j specifies the number N of machine (try -j4) cores so it compiles faster

`step 3:` Go up one directory (cd .. ),  and clone and setup the Hall C replay repository <br>
>  # Clone and properly setup the Hall C data analysis replay <br>
>     >> 