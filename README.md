# CaFe Online Replay
This repository, `cafe_online_replay`, is a version of the official `hallc_replay` repository with additions specific to the CaFe experiment. This online version has replay scripts for data online monitoring of the data, as well as specialized shell scripts for analysis of special runs such as target_boiling, H(e,e'p) elastics, production, detector calibrations, etc. 

Please refer to the following useful link:
[https://hallcweb.jlab.org/wiki/index.php/CaFe\_Experiment](https://hallcweb.jlab.org/wiki/index.php/CaFe_Experiment)

where you can find helpful information related to the experiment. 
# How-To Guide for Starters
This guide is intended for *novice level* on hall c analysis, using the CaFe experiment as an example. NOTE: In the exercises below, the semi-colons ';' denote a comment, and the '>>' denotes the command the user MUST input into the command-line terminal.

### Exercise #1:
Set up and get familiar with the general hall c analysis replay structure on ifarm. <br>

`step 1:` Assuming you have a Jefferson Lab account, log-in to ifarm: 
> ; login to ifarm<br>
>     >> ssh -Y *user*@login.jlab.org <br>
>     >> ssh -Y ifarm <br>
>; setup necessary environment variables <br>
>     >> source /site/12gev_phys/softenv.csh 2.5 

`step 2:` Go to the relevant work directory and setup the relevant Hall C Analyzer repository. 
>  ; Create symbolic link to the CaFe work directory. For more info see [https://hallcweb.jlab.org/wiki/index.php/CaFe\_Disk\_Space](https://hallcweb.jlab.org/wiki/index.php/CaFe_Disk_Space) <br>
>     >> ln -s /w/hallc-scshelf2102/c-cafe-2022 cafe\_work <br>
>     >> cd cafe\_work <br>
>
>  ; If you don't have a user directory, make one and *cd* to it !<br>
>     >> mkdir *user* <br> 
>     >> cd *user*<br>
> 
>  ; Clone and properly setup the Hall C analyzer source code <br>
>     >> git clone https://github.com/JeffersonLab/hcana <br>
>     >> cd hcana <br> 
> 
>    ; execute the three commands below separately <br>
>    ; (usually only required after cloning hcana, but if there are updates to Hall A analyzer (podd), then it will need to be updated.) <br>
>     >> git submodule init; git submodule  sync; git submodule update <br>
> 
>  ; setup hcana environment variables <br>
>     >> source setup.csh   <br> 
> 
>   ; the command below will compile the analyzer and generate an executable, *hcana*<br>
>     >> scons -jN # scons is a compiler and -j specifies the number N of machine (try -j4) cores so it compiles faster

`step 3:` Go up one directory (cd .. ),  and clone and setup the Hall C replay repository <br>
>  ; Clone and properly setup the Hall C data analysis replay <br>
>  ; Alternatively, you can fork a copy of the repository remotely and the clone it directly from your github account. <br>
>     >> git clone [https://github.com/Yero1990/cafe\_online\_replay](https://github.com/Yero1990/cafe_online_replay) <br>
>     >> cd cafe\_online\_replay <br>

>  ; execute this script to create the necessary sybmolic links required by the replay script <br>
>     >> ./cafe\_setup.sh 

`step 4:` Try to run your first replay on a sample data file (you might get errors. This is a learning process with a majority of the time spent on de-bugging code. Don't feel bad about it.)
>  ; Execute a replay script (reads a raw data file and generates an output ROOTfile and REPORT\_FILE) <br>
>  ; Follow the command-line requests to enter run number, number of events as well as <br>
>  ; an additional option, depending on the analysis you want to carry out. The last option, <br>
>  ; I added specially for CaFe, to faciliate our analysis (on a normal hallc_replay data analysis, <br> 
>  ; you will typically only be asked for a run and event number) <br>
>     >> ./hcana replay\_cafe.C <br>
 
# How-To Guide for Active Contributors
This guide is intended for users who would like to actively contribute to the development of this repository. <br>
See [https://docs.github.com/en/get-started/quickstart/github-flow](https://docs.github.com/en/get-started/quickstart/github-flow) for helpful GitHub documentaion 

`step 1:` Go to [https://github.com/Yero1990/cafe\_online\_replay](https://github.com/Yero1990/cafe_online_replay) and fork a copy to your own repo.<br>

`step 2:` Clone the forked repository from your own repo and set up upstream track to keep track of any updates made to the official  repository. <br> 
> ; clone the repository either locally or remotely (wherever you plan to work)<br>
>     >> git clone https://github.com/your\_github\_username/cafe\_online\_replay <br>
>     >> cd cafe\_online\_replay <br>

> ; track the upstram branch is needed to keep your forked copy and local machine copy up-to-date <br>
>     >> git remote add --track master upstream https://github.com/Yero1990/cafe\_online\_replay <br>
>     HINT: if you do: git remote -v,  it will specify the track to the remote branches origin and upstream<br>

`step 3:` To make any contributions locally and save it to the official [https://github.com/Yero1990/cafe\_online\_replay](https://github.com/Yero1990/cafe_online_replay) repository:
> ; change to a new working branch on your local copy of the repo, for example: <br>
>     >> git checkout -b work\_branch <br>

> ; make changes/contributions while on the work\_branch and then add the changes to a staging area <br>
> ; (i.e., putting all your changes into a common area, for saving later) <br>
>     >> git add file1.ext file2.ext . . .  ; where fileN.ext represnets the files added or modified <br>

> ; save (i.e. commit) the changes locally on your work\_branch, the -m flag it for adding a commit message (obligatory) <br>
>     >> git commit -m "message stating a summary of your commits. " <br>

> ; save the changes to your remote forked copy of the cafe\_online\_replay repo <br>
>     >> git push origin work\_branch

> ; make a pull request to incorporate the changes into the official cafe\_online\_replay <br>

>  a) Go to:  https://github.com/your\_github\_username/cafe\_online\_replay and select the work\_branch (top left) <br>

>  b) Click on 'Contribute' and select 'Open pull request' <br>

>  c) Add a brief description of your work on the commnet box and click 'Create pull request' <br>

>  d) The person with admin priviliges will receive a notification and will check the work, and if there is no conflict with other files, then the 'Pull request' will be accepted and merged onto the official branch <br>

>  e) To update the changes locally,  and remotely do the following: <br>
>      >> git checkout master ; change back to the 'master' branch <br> 
>      >> git pull upstream master ; to pull the changes just merged to upstream down to your master branch <br>
>      >> git push origin master ; to push the changes from your local 'master' branch up to your remote origin 'master' branch <br>
> 
AND THE CYCLE COMPLETES ! ! ! <br>

`step 4:` Now you can safely remove the work branch locally and remotely and repeat the cycle in step 3. <br>
>   ; remove local copy of branch (optional) <br>
>     >> git branch -d work\_branch
