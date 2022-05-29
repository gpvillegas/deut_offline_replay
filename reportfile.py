#! /usr/bin/python

# Stephen JD Kay - 19/08/21 - University of Regina
# Script to grab info from a report file, for use in the run list shell script
# This script scans every line of the report file to find the correct info

# C. Yero - May 26, 2022
# Modified to meet the needs of CaFe 2022 Experiment

# Import relevant packages
import sys, math, os, subprocess

sys.path.insert(0, 'python/')

if len(sys.argv)-1!=2:
    print("Invalid number of arguments. \n 
    e.g., python reportfile.py <path/to/cafe_report_run.txt> <entry_type> \n
    <entry_type> = gen_run_info, trig_info, eff_info, good_evt_info" )
    sys.exit(1)


# user input
cafe_report_path = sys.argv[1]
entry_type = sys.argv[2]  # <entry_type> = "gen_run", "trig_info", "eff_info", "good_evt_info"

cafe_report = open(cafe_report_path)


# general run info
run_num=0
kin_type=0
run_len=0
evt_num=0
beam_e=0
tgt_name=0
tgt_mass=0
hms_p=0
hms_angle=0
shms_p=0
shms_angle=0
bcm_thrs=0
bcm_current=0
bcm_charge=0

# good events counts
heep_singles      =-1
heep_singles_rate =-1
heep_real         =-1
heep_real_rate    =-1
MF_real           =-1
MF_real_rate      =-1
SRC_real          =-1
SRC_real_rate     =-1

# trigger info (only enabled triggers, i.e PS# != -1 will be written to kin file)
PS1=-1    # SHMS 3/4
PS2=-1    # SHMS EL-REAL
PS3=-1    # HMS 3/4
PS5=-1    # SHMS EL-REAL x HMS 3/4 ( this needs to be implemented after pionLT run ends in August 2022 )
PS6=-1    # HMS 3/4 x SHMS 3/4

T1_scaler=-1
T2_scaler=-1
T3_scaler=-1
T5_scaler=-1
T6_scaler=-1

T1_scaler_rate=-1
T2_scaler_rate=-1
T3_scaler_rate=-1
T5_scaler_rate=-1
T6_scaler_rate=-1

T1_accp=-1
T2_accp=-1
T3_accp=-1
T5_accp=-1
T6_accp=-1

T1_accp_rate=-1
T2_accp_rate=-1
T3_accp_rate=-1
T5_accp_rate=-1
T6_accp_rate=-1

# tracking efficiency
hms_trk_eff=-1
shms_trk_eff=-1

# daq live time
T1_cpuLT=-1
T1_tLT=-1

T2_cpuLT=-1
T2_tLT=-1

T3_cpuLT=-1
T3_tLT=-1

T5_cpuLT=-1
T5_tLT=-1

T6_cpuLT=-1
T6_tLT=-1

TestVar = 0 # Counter to check the right number of variables have been set, 
for line in cafe_report:
    if (line[0]=="#") continue;

    # general run information
    if "run_number" in line :
        run_num = int((line.split(":")[1]).strip())
        TestVar+=1
        print(run_num)
    if "kin_type" in line :
        kin_type = (line.split(":")[1]).strip()
        TestVar+=1
        print(kin_type)
    if "daq_run_length" in line :
        run_len = float((line.split(":")[1]).strip())
        TestVar+=1
        print(run_len)
    if "events_replayed" in line :
        evt_num = int((line.split(":")[1]).strip())
        TestVar+=1
        print(evt_num)
    if "beam_energy" in line :
        beam_e = float((line.split(":")[1]).strip())
        TestVar+=1
        print(beam_e)
    if "target_name:" in line :
        tgt_name = (line.split(":")[1]).strip()
        TestVar+=1
        print(tgt_name)
    if "target_amu" in line :
        tgt_mass = float((line.split(":")[1]).strip())
        TestVar+=1
        print(tgt_mass)
    if "hms_h_momentum" in line :
        hms_p = float((line.split(":")[1]).strip())
        TestVar+=1
        print(hms_p)
    if "hms_h_angle" in line :
        hms_angle = float((line.split(":")[1]).strip())
        TestVar+=1
        print(hms_angle)
    if "shms_e_momentum" in line :
        shms_p = float((line.split(":")[1]).strip())
        TestVar+=1
        print(shms_p)
    if "shms_e_angle" in line :
        shms_angle = float((line.split(":")[1]).strip())
        TestVar+=1
        print(shms_angle)
    if "_Current_Threshold" in line :
        bcm_thrs = (line.split(":")[1]).strip()
        TestVar+=1
        print(bcm_thrs)
    if "_Average_Current" in line :
        bcm_current = float((line.split(":")[1]).strip())
        TestVar+=1
        print(bcm_current)
    if "_Charge" in line :
        bcm_charge = float((line.split(":")[1]).strip())
        TestVar+=1
        print(bcm_charge)
        
    # good elastic, MF or SRC counts/rates
    if "heep_total_singles_counts" in line :
        heep_singles = float((line.split(":")[1]).strip())
        TestVar+=1
        print(heep_singles)
    if "heep_total_singles_rate" in line :
        heep_singles_rate = float((line.split(":")[1]).strip())
        TestVar+=1
        print(heep_singles_rate)
    if "heep_real_counts" in line :
        heep_real = float((line.split(":")[1]).strip())
        TestVar+=1
        print(heep_real)
    if "heep_real_rate" in line :
        heep_real_rate = float((line.split(":")[1]).strip())
        TestVar+=1
        print(heep_real_rate)
    if "MF_real_counts" in line :
        MF_real = float((line.split(":")[1]).strip())
        TestVar+=1
        print(MF_real)
    if "MF_real_rate" in line :
        MF_real_rate = float((line.split(":")[1]).strip())
        TestVar+=1
        print(MF_real_rate)
    if "SRC_real_counts" in line :
        SRC_real = float((line.split(":")[1]).strip())
        TestVar+=1
        print(SRC_real)
    if "SRC_real_rate" in line :
        SRC_real_rate = float((line.split(":")[1]).strip())
        TestVar+=1
        print(SRC_real_rate)

    # trigger info
    if "Ps1_factor" in line :
        PS1 = int((line.split(":")[1]).strip())
        TestVar+=1
        print(PS1)
    if "Ps2_factor" in line :
        PS2 = int((line.split(":")[1]).strip())
        TestVar+=1
        print(PS2)
    if "Ps3_factor" in line :
        PS3 = int((line.split(":")[1]).strip())
        TestVar+=1
        print(PS3)
    if "Ps5_factor" in line :
        PS5 = int((line.split(":")[1]).strip())
        TestVar+=1
        print(PS5)
    if "Ps6_factor" in line :
        PS6 = int((line.split(":")[1]).strip())
        TestVar+=1
        print(PS6)
        
    # scaler counts
    if "T1_scaler" in line :
        T1_scaler = int((line.split(":")[1]).split("[")[0].strip())
        TestVar+=1
        print(T1_scaler)
    if "T2_scaler" in line :
        T2_scaler = int((line.split(":")[1]).split("[")[0].strip())
        TestVar+=1
        print(T2_scaler)
    if "T3_scaler" in line :
        T3_scaler = int((line.split(":")[1]).split("[")[0].strip())
        TestVar+=1
        print(T3_scaler)
    if "T5_scaler" in line :
        T5_scaler = int((line.split(":")[1]).split("[")[0].strip())
        TestVar+=1
        print(T5_scaler)
    if "T6_scaler" in line :
        T6_scaler = int((line.split(":")[1]).split("[")[0].strip())
        TestVar+=1
        print(T6_scaler)
        
    # scaler rates (same line, but get str within '[,]' brackets
    if "T1_scaler" in line :
        T1_scaler_rate = float((line.split("[")[1]).split("kHz")[0].strip())
        TestVar+=1
        print(T1_scaler_rate)
    if "T2_scaler" in line :
        T2_scaler_rate = float((line.split("[")[1]).split("kHz")[0].strip())
        TestVar+=1
        print(T2_scaler_rate)
    if "T3_scaler" in line :
        T3_scaler_rate = float((line.split("[")[1]).split("kHz")[0].strip())
        TestVar+=1
        print(T3_scaler_rate)
    if "T5_scaler" in line :
        T5_scaler_rate = float((line.split("[")[1]).split("kHz")[0].strip())
        TestVar+=1
        print(T5_scaler_rate)
    if "T6_scaler" in line :
        T6_scaler_rate = float((line.split("[")[1]).split("kHz")[0].strip())
        TestVar+=1
        print(T6_scaler_rate)

    # accepted counts
    if "T1_accepted" in line :
        T1_accp = int((line.split(":")[1]).split("[")[0].strip())
        TestVar+=1
        print(T1_accp)
    if "T2_accepted" in line :
        T2_accp = int((line.split(":")[1]).split("[")[0].strip())
        TestVar+=1
        print(T2_accp)
    if "T3_accepted" in line :
        T3_accp = int((line.split(":")[1]).split("[")[0].strip())
        TestVar+=1
        print(T3_accp)
    if "T5_accepted" in line :
        T5_accp = int((line.split(":")[1]).split("[")[0].strip())
        TestVar+=1
        print(T5_accp)
    if "T6_accepted" in line :
        T6_accp = int((line.split(":")[1]).split("[")[0].strip())
        TestVar+=1
        print(T6_accp)
        
    # accepted rates (same line, but get str within '[,]' brackets
    if "T1_accepted" in line :
        T1_accp_rate = float((line.split("[")[1]).split("kHz")[0].strip())
        TestVar+=1
        print(T1_accp_rate)
    if "T2_accepted" in line :
        T2_accp_rate = float((line.split("[")[1]).split("kHz")[0].strip())
        TestVar+=1
        print(T2_accp_rate)
    if "T3_accepted" in line :
        T3_accp_rate = float((line.split("[")[1]).split("kHz")[0].strip())
        TestVar+=1
        print(T3_accp_rate)
    if "T5_accepted" in line :
        T5_accp_rate = float((line.split("[")[1]).split("kHz")[0].strip())
        TestVar+=1
        print(T5_accp_rate)
    if "T6_accepted" in line :
        T6_accp_rate = float((line.split("[")[1]).split("kHz")[0].strip())
        TestVar+=1
        print(T6_accp_rate)

    # tracking efficiency
    if "hms_had_track_eff" in line :
        hms_trk_eff = float(line.split(":")[1].split("+")[0].strip())
        print(hms_trk_eff)
    if "shms_elec_track_eff" in line :
        shms_trk_eff = float(line.split(":")[1].split("+")[0].strip())
        print(shms_trk_eff)
        
    # live times
    if "T1_cpuLT" in line :
        T1_cpuLT = float(line.split(":")[1].split("+")[0].strip())
        print(T1_cpuLT)
    if "T1_tLT" in line :
        T1_tLT = float(line.split(":")[1].split("+")[0].strip())
        print(T1_tLT)

    if "T2_cpuLT" in line :
        T2_cpuLT = float(line.split(":")[1].split("+")[0].strip())
        print(T2_cpuLT)
    if "T2_tLT" in line :
        T2_tLT = float(line.split(":")[1].split("+")[0].strip())
        print(T2_tLT)

    if "T3_cpuLT" in line :
        T3_cpuLT = float(line.split(":")[1].split("+")[0].strip())
        print(T3_cpuLT)
    if "T3_tLT" in line :
        T3_tLT = float(line.split(":")[1].split("+")[0].strip())
        print(T3_tLT)

    if "T5_cpuLT" in line :
        T5_cpuLT = float(line.split(":")[1].split("+")[0].strip())
        print(T5_cpuLT)
    if "T5_tLT" in line :
        T5_tLT = float(line.split(":")[1].split("+")[0].strip())
        print(T5_tLT)

    if "T6_cpuLT" in line :
        T6_cpuLT = float(line.split(":")[1].split("+")[0].strip())
        print(T6_cpuLT)
    if "T6_tLT" in line :
        T6_tLT = float(line.split(":")[1].split("+")[0].strip())
        print(T6_tLT)

# general run entry list
gen_run_info = "%i       %s        %.3f     %i       %.4f     %s        %.6f      %.4f   %.3f       %.4f    %.3f        %s        %.3f         %.3f       " % \
               (run_num, kin_type, run_len, evt_num, beam_e,  tgt_name, tgt_mass, hms_p, hms_angle, shms_p, shms_angle, bcm_thrs, bcm_current, bcm_charge)

# trigger info
trig_info = "%i   %i   %i   %i   %i   %.3f            %.3f            %.3f            %.3f            %.3f            %.3f          %.3f          %.3f           %.3f           %.3f                        "  % \
            (PS1, PS2, PS3, PS5, PS6, T1_scaler_rate, T2_scaler_rate, T3_scaler_rate, T5_scaler_rate, T6_scaler_rate, T1_accp_rate, T2_accp_rate, T3_accp_rate,  T5_accp_rate,  T6_accp_rate                          )

# live time and trk_eff info
efficiency_info = "%.3f    %.3f    %.3f    %.3f    %.3f    %.3f         %.3f        " % \
           (T1_tLT, T2_tLT, T3_tLT, T5_tLT, T6_tLT, hms_trk_eff, shms_trk_eff )

# good event count info
good_evt_info = "%.2f           %.3f               %.2f       %.3f            %.2f     %.3f          %.2f      %.3f    " % \
                (heep_singles,  heep_singles_rate, heep_real, heep_real_rate, MF_real, MF_real_rate, SRC_real, SRC_real_rate )

if(entry_type == "gen_run_info"):
    print(gen_run_info)
elif(entry_type == "trig_info"):
    print(trig_info)
elif(entry_type == "eff_info"):
    print(efficiency_info)
elif(entry_type == "good_evt_info"):
    print(good_evt_info)
else:
    print("Invalid <entry_type> \n 
    e.g., python reportfile.py <path/to/cafe_report_run.txt> <entry_type> \n
    <entry_type> = gen_run_info, trig_info, eff_info, good_evt_info")

cafe_report.close()
