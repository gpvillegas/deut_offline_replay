# Deuteron Online Replay
This repository, `deut_online_replay`, is a version of the official `hallc_replay` repository with additions specific to the CaFe experiment. This online version has replay scripts for data online monitoring of the data, as well as specialized shell scripts for analysis of special runs such as target_boiling, H(e,e'p) elastics, production, detector calibrations, etc. 

Please refer to the following useful link:
[https://hallcweb.jlab.org/wiki/index.php/CaFe\_Experiment](https://hallcweb.jlab.org/wiki/index.php/CaFe_Experiment)

where you can find helpful information related to the experiment. 
# How-To Guide for Starters
This guide is intended for *novice level* on hall c analysis, using the CaFe experiment as an example. NOTE: In the exercises below, the pound sign, `# `, denotes comments, and the `$`
denotes the user terminal command-line input

### Exercise #1:
Set up and get familiar with the general hall c analysis replay structure on ifarm. <br>

`step 1:` Assuming you have a Jefferson Lab account, log-in to ifarm: 

```sh
# login to ifarm
$ ssh -Y user@login.jlab.org 
$ ssh -Y ifarm 
```

`step 2:` Go to the relevant work directory and setup the relevant Hall C Analyzer repository. 
##### (ignore the `git clone` command if using this repo as a submodule of `CaFe-Online`)

```sh
# Create symbolic link to the CaFe work directory. 
# For more info see [https://hallcweb.jlab.org/wiki/index.php/CaFe_Disk_Space]
$ ln -s /work/hallc/c-cafe-2022/ cafe_work 
$ cd cafe_work 

# If you don't have a user directory, make one and cd to it !
$ mkdir <user> 
$ cd <user>

# Clone and properly setup the Hall C analyzer source code <br>
$ git clone https://github.com/JeffersonLab/hcana
$ cd hcana 
 
# initialize and update the Hall A Analyzer (podd), which is itself a submodule of hcana
# usually only requires ONCE, but if podd is updated, then submodule will need to be 
# updated as well (extremely rare)
$ git submodule update --init --recursive  

# setup hcana environment variables 
# tip: add this command to your (.bashrc, .cshrc, etc.) so its done automatically at the 
# start of each terminal session
$ source setup.csh   
 
# the command below will compile the analyzer and generate an executable, hcana
# scons is a compiler and -j specifies the number N of machine (try -j4) cores so it compiles faster
$ scons -jN 
```

`step 3:` Go up one directory (cd .. ),  and clone and setup the Hall C replay repository <br>
##### (ignore the `git clone` command if using this repo as a submodule of `CaFe-Online`)
```sh
# Clone and properly setup the Hall C data analysis replay 
# Alternatively, you can fork a copy of the repository remotely and the clone it directly 
# from your github account.
$ git clone https://github.com/Yero1990/deut_online_replay 
$ cd cafe_online_replay 

# execute this script to create the necessary sybmolic links required by the replay script 
$ ./cafe_setup.sh 
```

```sh
# You will need to set-up the following direcotries, either via symbolic link or normal "mkdir" command
# for example, if on ifarm, I've create the following directories, and made a symbolic link to this repo.

mkdir /volatile/hallc/c-deuteron/cyero
mkdir /volatile/hallc/c-deuteron/cyero/ROOTfiles
mkdir /volatile/hallc/c-deuteron/cyero/REPORT_OUTPUT
mkdir /volatile/hallc/c-deuteron/cyero/DEUT_OUTPUT
mkdir /volatile/hallc/c-deuteron/cyero/DEUT_OUTPUT/ROOT
mkdir /volatile/hallc/c-deuteron/cyero/DEUT_OUTPUT/REPORT
mkdir /volatile/hallc/c-deuteron/cyero/DEUT_OUTPUT/PDF

# then for symbolic links, we only need to symbolic link the following:
ln -sf /volatile/hallc/c-deuteron/cyero/worksim
ln -sf /volatile/hallc/c-deuteron/cyero/ROOTfiles
ln -sf /volatile/hallc/c-deuteron/cyero/REPORT_OUTPUT 
ln -sf /volatile/hallc/c-deuteron/cyero/DEUT_OUTPUT

# in volatile for testing my analysis/initial offline analysis,         
# but keep in mind files are not backed up in this directory if the limit is exceeded. 
# you can make these directories anywhere you like, provided there is space / permissions
# also make sure to point to the directory where the .raw data files are located for replay
# I usually put a symbolic link to the .raw data files under the CACHE_LINKS/ directory
```


`step 4:` Try to run your first replay on a sample data file (you might get errors. This is a learning process with a majority of the time spent on de-bugging code. Don't feel bad about it.)

```sh
# Execute a replay script (reads a raw data file and generates an output ROOTfile and REPORT_FILE)
# Follow the command-line requests to enter run number, number of events as well as 
# an additional option, depending on the analysis you want to carry out. The last option, 
# I added specially for deuteron, to faciliate our analysis. To test it, execute:
$ ./hcana SCRIPTS/COIN/PRODUCTION/replay_deut.C  

# the script will ask you to enter a specific run number, event number and analysis type to use. Based on the
# analysis type you choose, it will read the corresponding definition file located in: DEF-files/ . Each definition
# file specifies which leaf variables to add. For specific detector calibrations / studies, I've added specific definition
# files which will only use the required leaf variables in order to keep the TTree ntuple small as possible so that it may
# run faster. 
# 
```
 
# How-To Guide for Active Contributors
This guide is intended for users who would like to actively contribute to the development of this repository. <br>
See [https://docs.github.com/en/get-started/quickstart/github-flow](https://docs.github.com/en/get-started/quickstart/github-flow) for helpful GitHub documentaion 

`step 1:` Go to [https://github.com/Yero1990/cafe\_online\_replay](https://github.com/Yero1990/cafe_online_replay) and fork a copy to your own repo.<br>

`step 2:` Clone the forked repository from your own repo and set up upstream track to keep track of any updates made to the official  repository. <br> 

```sh
# clone the repository either locally or remotely (wherever you plan to work)
$ git clone https://github.com/your_github_username/cafe_online_replay 
$ cd cafe_online_replay 

# track the upstram branch is needed to keep your forked copy and local machine copy up-to-date 
$ git remote add --track master upstream https://github.com/Yero1990/cafe_online_replay 

# check that the track to the remote branches origin and upstream are specified
$ git remote -v  
```

`step 3:` To make any contributions locally and save it to the official [https://github.com/Yero1990/cafe\_online\_replay](https://github.com/Yero1990/cafe_online_replay) repository, do:

```sh
# change to a new working branch on your local copy of the repo, for example:
$ git checkout -b work_branch 

# make changes/contributions while on the work_branch and then add the changes to a staging area 
# (i.e., putting all your changes into a common area, for saving later) 
$ git add file1.ext file2.ext . . .  # where fileN.ext represnets the files added or modified 

# save (i.e. commit) the changes locally on your work_branch, the -m flag 
# is for adding a commit message (obligatory) 
$ git commit -m "message stating a summary of your commits. "

# save the changes to your remote forked copy of the cafe_online_replay repo 
$ git push origin work_branch

# make a pull request to incorporate the changes into the official cafe_online_replay

   a) Go to:  https://github.com/your_github_username/cafe_online_replay and select the work_branch

   b) Click on 'Contribute' and select 'Open pull request'. Alternatively, click on 'Compare & pull request'

   c) Add a brief description of your work on the commnet box and click 'Create pull request' 

   d) The person with admin priviliges will receive a notification and will check the work, and if there is 
      no conflict with other files, then the 'Pull request' will be accepted and merged onto the official 
      branch (at this point you will be given the option to delete the branch remotely) 

   e) To update the changes locally,  and remotely do the following: 
      $ git checkout master       # change back to the 'master' branch  
      $ git pull upstream master  # pull the changes just merged to upstream down to your master branch
      $ git push origin master    # push the changes to your remote origin 'master' branch  
                                  # AND THE CYCLE COMPLETES ! ! ! 
```

`step 4:` Now you can safely remove the work branch locally and remotely and repeat the cycle in step 3. <br>

```sh
# remove local copy of branch (optional)
$ git branch -d work_branch

# remove the remote copy of the branch (optional, if you haven't done so) 
$ git push origin --delete work_branch 
```
