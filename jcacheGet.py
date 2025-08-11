import os

run_nums = open("./run-list.dat","r")
run_list = [r.strip() for r in run_nums]
#print(run_list)
#run_nums = [3283,3284,3285,3286,3287]
#run_nums = [20840,20841,20846,20851,20858,20861,20868,20869]
#run_nums = [20886]
spec = 'shms'       # shms, hms, or coin
exp = 'c-deuteron'  # c-deuteron, c-cafe-2022, spring17 (deuteron 2018), ...

i = 0
for run in run_list:
    tapeDIR = f'/mss/hallc/{exp}/raw/{spec}_all_{run}.dat'   

    if i == 0:
        jcacheCMD = f'jcache get {tapeDIR} -e gvill@jlab.org -D 60'
        i+=1
    elif i == (len(run_list)-1):
        jcacheCMD = f'jcache get {tapeDIR} -e gvill@jlab.org -D 60'
    else:
        jcacheCMD = f'jcache get {tapeDIR} -D 60'
        i+=1    
    
    #print(jcacheCMD)
    os.system(jcacheCMD)
run_nums.close()

