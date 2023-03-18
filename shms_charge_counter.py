#!/usr/bin/env python
import time,sys,os,argparse,atexit,subprocess,math


XVALS,YVALS=[],[]
START=0
NCHAN=28


def getScalars():
    zvals=[]
    pzvals=[]
    chan=[]
    overflow=[]
    pchan=[]
    XVALS[:] = []
    YVALS[:] = []
    
    #file = open("output.txt",'r')
    #row = file.readlines()
    lines = subprocess.check_output(["getscalers","hcvme08"])
    for line in lines.splitlines():
        #print line
    #for line in row:
        info = line.split()
        info[0] = info[0].replace(':', '')
        if info[0].isdigit():
           chan.append(int (info[0])+1)
           temp  = info[1].rpartition('(')[0]
           temp2 =  info[1].rpartition('(')[2]
           temp3 =  temp2.rpartition(')')[0]
           overflow.append(int(temp3))
           zvals.append(float(temp))
           
    pzvals.append(zvals[255])
    pzvals.append(zvals[249])
    pzvals.append(zvals[192])
    pzvals.append(zvals[193])
    pzvals.append(zvals[194])
    return [pzvals,overflow]


def calcRates(rates,rates2):
    left,right,maximum=0,0,0
    for ii in range(28):#len(XVALS)):
        xx,yy,zz,zz2=XVALS[ii],YVALS[ii],rates[ii],rates2[ii]
        if xx==1:left+=(zz2-zz)/10
        else : right+=(zz2-zz)/10
        if (zz2-zz)/10>maximum: maximum=(zz2-zz)/10
    return [left,right,maximum]


def makeTime():
    datime=TDatime()
    return '%d/%d/%d %.2d:%.2d:%.2d'%(datime.GetDay(),
             datime.GetMonth(),
             datime.GetYear(),
             datime.GetHour(),
             datime.GetMinute(),
             datime.GetSecond())





def main():


    timeCul = 0
    timeCulCh = 0
    chargeCul = 0
    aveI = 0
    FracTime = 0
    while True:

        [zvals,overflow1] = getScalars()
        time.sleep(5)
        [zvals2,overflow2] = getScalars()
        diff = zvals2[0]-zvals[0]
        if diff>0 :
          time1 = (zvals2[0]-zvals[0])/1000000.
          current = ((zvals2[1]-zvals[1])/time1+3290)/9235.
          t1rate=  (zvals2[2]-zvals[2])/time1/1000.
          t2rate=  (zvals2[3]-zvals[3])/time1/1000.
          t3rate=  (zvals2[4]-zvals[4])/time1/1000.
          charge = 0
          timeCul = timeCul + time1/60.
          print('Current = {0:10.5f} '.format(current))
          if current>5. : 
              charge = current*time1
              chargeCul = chargeCul + charge/1000.
              timeCulCh = timeCulCh + time1/60.
          if timeCulCh>0: aveI=chargeCul/(timeCulCh*60)*1000
          if timeCul>0: FracTime=timeCulCh/timeCul*100
          print('Total Time = {0:6.2f} min, Beam on time = {1:6.2f} min.,  Frac = {4:7.2f} %  , Total Charge = {2:6.2f} mC, Ave I = {3:6.2f} uA'.format(timeCul,timeCulCh,chargeCul,aveI,FracTime))
          print('SHMS T1 rate kHz = {0:10.5f} , SHMS T2 Rate kHz = {1:10.5f}, SHMS T3 Rate kHz = {2:10.5f}'.format(t1rate,t2rate,t3rate))
        else :
          if diff==0:
            print(' Time diff = {0} {1}, Assume between runs and reset counters'.format(zvals2[0],zvals[0]))
            timeCul = 0
            chargeCul = 0
            timeCulCh = 0

if __name__=='__main__': main()


