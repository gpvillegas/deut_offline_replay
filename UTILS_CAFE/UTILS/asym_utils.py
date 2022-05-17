#Various python functions to aid in asymmetry calculations
#for beam-spin asymmetry analysis of Kaon LT (2018)

import math
from scipy import optimize
import LT.box as B
from LT.datafile import dfile
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import numpy.ma as ma
import sys                                     
import os                                                                                                       
from sys import argv  
import matplotlib
from matplotlib import rc
from matplotlib import *
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)


#Use latex commands (e.g. \textit ot \textbf)
rc('text', usetex=True)
#Set default font to times new roman
font = {'family' : 'Times New Roman',
        'weight' : 'normal',
        'size'   : 12
}
plt.rc('font', **font)

#Set font
csfont = {'fontname':'Times New Roman'}

def calc_asym(particle="Kaon",  kin_group="KIN-1"):

    # Brief: Read +/- helicity datafiles and calculate asymmetry 'per KIN group'
    # An output file is created with the relevant bin information and asymmetry
    # for 'Kaons' or 'Pions'

    # Code Usage: calc_asym("Pion", "KIN-5")
    # The code will produce the output file: asymmetry_Pion_KIN-5.txt
    
    #The central values are: thxq_cm = [2.5, 7.5, 12.5 17.5, 22.5, 27.5, 32.5, 37.5] +/- 2.5 deg
        
    #Read datafiles (extracted from 2D bin information of thxq_CM_vs_phxq)
    fname_pos = '../BeamPolAsymmetry/analysis_scripts/OUTPUT/%s/pos_hel_%s_combined.txt'%(kin_group, particle) 
    fname_neg = '../BeamPolAsymmetry/analysis_scripts/OUTPUT/%s/neg_hel_%s_combined.txt'%(kin_group, particle)  

    asy_pos = dfile(fname_pos)
    asy_neg = dfile(fname_neg)

    #-----Read general binning information (same for +/- helicities)------
    
    ib = np.array(asy_pos['ib'])                          #2d bin
    xb = np.array(asy_pos['xb'])                          #central bin value (shouls be the same for + or - hel)
    yb = np.array(asy_pos['yb'])                          #central bin value (shouls be the same for + or - hel)

    xlow =  np.array(asy_pos['xlow'])
    xup =  np.array(asy_pos['xup'])
    
    ylow =  np.array(asy_pos['ylow'])
    yup =  np.array(asy_pos['yup'])
        
    th_xqCM_c = np.array(asy_pos['y0'])                #central bin value (shouls be the same for + or - hel)
    ph_xq_c = np.array(asy_pos['x0'])                  #central bin value (shouls be the same for + or - hel)

    xbins = len(ph_xq_c[th_xqCM_c==th_xqCM_c[0]])
    ybins = len(th_xqCM_c[ph_xq_c==ph_xq_c[0]])

    xbin_width = xup[0] - xlow[0]
    ybin_width = yup[0] - ylow[0]
    
    #---Read positive helicity/error---
    N_pos = np.array(asy_pos['zcont'])
    N_pos_err = np.array(asy_pos['zcont_err'])
    
    #---Read negative helicity/error---
    N_neg = np.array(asy_neg['zcont'])
    N_neg_err = np.array(asy_neg['zcont_err'])

    #---Create output file to write helicities/asymmetry---
    fout_name = "asymmetry_%s_%s.txt" % (particle, kin_group)
    header = """# %s Beam-Spin Asymmetry (%s)
#
# histogram parameters:
# ybins      =  %i    | th_xq_CM [deg]
# xbins      =  %i    | phi_xq   [deg]
# ybin width =  %.3f  (or +/- %.3f deg)
# xbin width =  %.3f  (or +/- %.3f deg)
#
# header definitions:
# ib:         (x,y)  bin number 
# xb:         x-axis bin number 
# yb:         y-axis bin number 
# x0:         x-axis central bin value 
# xlow:       x-axis low-edge bin value 
# xup:        x-axis up-edge bin value 
# y0:         y-axis central bin value 
# ylow:       y-axis low-edge bin value 
# yup:        y-axis up-edge bin value 
# Np:         pos (+) helicity counts
# Np_err:     pos (+) helicity counts uncertainty
# Nn:         neg (-) helicity counts
# Nn_err:     neg (-) helicity counts uncertainty
# Asy:        Beam-Spin Asymmetry
# Asy_err:    Beam-Spin Asymmetry uncertainty
#
#! ib[i,0]/  xb[i,1]/  yb[i,2]/  x0[f,3]/  xlow[f,4]/  xup[f,5]/  y0[f,6]/  ylow[f,7]/  yup[f,8]/  Np[f,9]/  Np_err[f,10]/  Nn[f,11]/  Nn_err[f,12]/  Asy[f,13]/  Asy_err[f,14]/
""" % (particle, kin_group, ybins, xbins, ybin_width, ybin_width/2., xbin_width, xbin_width/2.)

    f = open(fout_name, "w")
    f.write(header)
        
    #Calculate the beam asymmetry and error per (th_xq_cm, ph_xq) bins
    for i, i_b in enumerate(ib):
    
        N = N_pos[i] + N_neg[i]  
        Asy = (N_pos[i] - N_neg[i]) / N 
        Asy_err2 = (2*N_neg[i] / N**2)**2 * N_pos_err[i]**2 + (2*N_pos[i]/N**2)**2 * N_neg_err[i]**2
        Asy_err = np.sqrt(Asy_err2)
        f.write("%i  %i  %i  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f \n"  % (i_b,  xb[i],  yb[i], ph_xq_c[i], xlow[i], xup[i], th_xqCM_c[i], ylow[i], yup[i], N_pos[i], N_pos_err[i], N_neg[i], N_neg_err[i], Asy, Asy_err) )
    f.close()

def calc_asym_combined(kin_group = np.array([]), particle="", group=""):

    # Brief: Read multiple kinematic files (KIN-1, KIN-2, . . .)
    # to combine different shms angles of same (Q2, W) kinematic setting
    # to gain better statistics. The output combined filename is set by the user.

    # Code Usage: calc_asym_combined(np.array(["KIN-1", "KIN-2", "KIN-3"]), "Kaon", "group-1")
    # The code will produce the output file (example): asymmetry_Kaon_group-1.txt
    
    # Helicity Counters (running sum over different KIN settings of the same (W,Q2))
    N_pos_tot     = 0. 
    N_pos_err_tot = 0.
    N_neg_tot     = 0.
    N_neg_err_tot = 0.
    
    #Loop over each KIN setting of kin_group to be combined
    for i, idx in enumerate(kin_group):

        #Get filename of ith KIN setting
        fname_pos = '../OUTPUT_Pass1/%s/pos_hel_%s_combined.txt'%(kin_group[i], particle) 
        fname_neg = '../OUTPUT_Pass1/%s/neg_hel_%s_combined.txt'%(kin_group[i], particle)  

        asy_pos = dfile(fname_pos)
        asy_neg = dfile(fname_neg)

        #Read  bin number 
        ib = np.array(asy_pos['ib'])                          #2d bin
        xb = np.array(asy_pos['xb'])                          #entral bin value (shouls be the same for + or - hel)
        yb = np.array(asy_pos['yb'])                          #entral bin value (shouls be the same for + or - hel)
        
        #Read central bin value 
        th_xqCM_c = np.array(asy_pos['y0'])                          #entral bin value (shouls be the same for + or - hel)
        ph_xq_c = np.array(asy_pos['x0'])                            #entral bin value (shouls be the same for + or - hel)
    
        #Read low/high bin edges
        xlow =  np.array(asy_pos['xlow'])
        xup =  np.array(asy_pos['xup'])
        
        ylow =  np.array(asy_pos['ylow'])
        yup =  np.array(asy_pos['yup'])

        xbins = len(ph_xq_c[th_xqCM_c==th_xqCM_c[0]])
        ybins = len(th_xqCM_c[ph_xq_c==ph_xq_c[0]])
        
        xbin_width = xup[0] - xlow[0]
        ybin_width = yup[0] - ylow[0]
        
        #Read positive helicity/error
        N_pos = np.array(asy_pos['zcont'])
        N_pos_err = np.array(asy_pos['zcont_err'])
        
        #Read negative helicity/error
        N_neg = np.array(asy_neg['zcont'])
        N_neg_err = np.array(asy_neg['zcont_err'])

        #Running Sum of +/- helicity states and standard error propagation of sums
        N_pos_tot = N_pos_tot + N_pos
        N_pos_err_tot = np.sqrt( N_pos_err_tot**2 + N_pos_err**2)

        N_neg_tot = N_neg_tot + N_neg
        N_neg_err_tot = np.sqrt( N_neg_err_tot**2 + N_neg_err**2)


    #Calculate Asymmetry/Error
    N = N_pos_tot + N_neg_tot
    Asy = (N_pos_tot - N_neg_tot) / N 
    Asy_err2 = (2*N_neg_tot/N**2)**2 * N_pos_err_tot**2 + (2*N_pos_tot/N**2)**2 * N_neg_err_tot**2
    Asy_err = np.sqrt(Asy_err2)
    
    #Write combined +/- helicities and calculated asymmetry to file
    fout_name = "asymmetry_%s_%s.txt" % (particle, group)
    header = """# %s Beam-Spin Asymmetry (%s)
#
# histogram parameters:
# ybins      =  %i    | th_xq_CM [deg]
# xbins      =  %i    | phi_xq   [deg]
# ybin width =  %.3f  (or +/- %.3f deg)
# xbin width =  %.3f  (or +/- %.3f deg)
#
# header definitions:
# ib:         (x,y)  bin number 
# xb:         x-axis bin number 
# yb:         y-axis bin number 
# x0:         x-axis central bin value 
# xlow:       x-axis low-edge bin value 
# xup:        x-axis up-edge bin value 
# y0:         y-axis central bin value 
# ylow:       y-axis low-edge bin value 
# yup:        y-axis up-edge bin value 
# Np:         pos (+) helicity counts
# Np_err:     pos (+) helicity counts uncertainty
# Nn:         neg (-) helicity counts
# Nn_err:     neg (-) helicity counts uncertainty
# Asy:        Beam-Spin Asymmetry
# Asy_err:    Beam-Spin Asymmetry uncertainty
#
#! ib[i,0]/  xb[i,1]/  yb[i,2]/  x0[f,3]/  xlow[f,4]/  xup[f,5]/  y0[f,6]/  ylow[f,7]/  yup[f,8]/  Np[f,9]/  Np_err[f,10]/  Nn[f,11]/  Nn_err[f,12]/  Asy[f,13]/  Asy_err[f,14]/
""" % (particle, kin_group, ybins, xbins, ybin_width, ybin_width/2., xbin_width, xbin_width/2.)
    
    fo = open(fout_name, "w")
    fo.write(header)
      
    for i, b in enumerate(ib):
        fo.write("%i  %i  %i  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f \n"  % (b,  xb[i],  yb[i], ph_xq_c[i], xlow[i], xup[i], th_xqCM_c[i], ylow[i], yup[i], N_pos_tot[i], N_pos_err_tot[i], N_neg_tot[i], N_neg_err_tot[i], Asy[i], Asy_err[i]) )
    fo.close()

def combine_thcm_bins(particle="", group=""):

    # Brief: Code to integrate (or combine) over all th_xq_cm bins.
    # The integration can be done either per KIN setting or per group (combined KIN settings)
    # where a group constitutes multiple KIN settings that have already been combined.

    # Code Usage (examples): combine_thcm_bins("Kaon", "KIN-1") or
    #                        combine_thcm_bins("Pion", "group-3")
    # The code will produce the output file (example):
    # asymmetry_Kaon_KIN-1.txt OR asymmetry_Pion_group-3.txt
    # These output files will have the intergated asymmetries over all th_xq_cm bins
    
    #Read in datafile (either per KIN or per group)
    fname = "asymmetry_%s_%s.txt" % (particle, group)

    f = dfile(fname)

    #-----Read general binning information (same for +/- helicities)------
    
    ib = np.array(f['ib'])                          #2d bin
    xb = np.array(f['xb'])                          #central bin value (should be the same for + or - hel)
    yb = np.array(f['yb'])                          #central bin value (should be the same for + or - hel)

    xlow =  np.array(f['xlow'])
    xup =  np.array(f['xup'])
    
    ylow =  np.array(f['ylow'])
    yup =  np.array(f['yup'])
        
    th_xqCM_c = np.array(f['y0'])                #central bin value (should be the same for + or - hel)
    ph_xq_c = np.array(f['x0'])                  #central bin value (should be the same for + or - hel)

    xbins = len(ph_xq_c[th_xqCM_c==th_xqCM_c[0]])    
    xbin_width = xup[0] - xlow[0]
        
    #---Read positive helicity/error---
    N_pos = np.array(f['Np'])
    N_pos_err = np.array(f['Np_err'])
    
    #---Read negative helicity/error---
    N_neg = np.array(f['Nn'])
    N_neg_err = np.array(f['Nn_err'])

    #Find pattern in th_xq_cm and ph_xq arrays by requiring the first element of the correlated array
    th_xq_cm_pat = th_xqCM_c[ph_xq_c==ph_xq_c[0]]
    ph_xq_pat = ph_xq_c[th_xqCM_c==th_xqCM_c[0]]

    #Write combined +/- helicities and calculated asymmetry to file
    fout_name = "asymmetry_%s_%s_allthCM.txt" % (particle, group)
    header = """# %s Beam-Spin Asymmetry (%s) :: Integrated over all th_xq_cm bins
#
# histogram parameters:
# xbins      =  %i    | phi_xq   [deg]
# xbin width =  %.3f  (or +/- %.3f deg)
#
# header definitions:
# ph_xq:      out-of-plane angle between detected particle (X) and q-vector [deg] 
# Np:         pos (+) helicity counts
# Np_err:     pos (+) helicity counts uncertainty
# Nn:         neg (-) helicity counts
# Nn_err:     neg (-) helicity counts uncertainty
# Asy:        Beam-Spin Asymmetry
# Asy_err:    Beam-Spin Asymmetry uncertainty
#
#! ph_xq[f,0]/  Np[f,1]/  Np_err[f,2]/  Nn[f,3]/  Nn_err[f,4]/  Asy[f,5]/  Asy_err[f,6]/
""" % (particle, group, xbins, xbin_width, xbin_width/2.)

    f = open(fout_name, "w")
    f.write(header)
    
    
    #To combine asymmetry for a single bin in phi_xq for multiple bins in (th_xq_cm),
    #we must loop over the pattern in phi_xq: [-180, -172, . . . 0, 180], and for each element,
    #sum over the + (or -) helicity states of each th_xq_cm bin and calculate the corresponding errors 
    for i, iph_xq in enumerate(ph_xq_pat):


        #Sum the helicities of all th_xq_cm bins  per ph_xq bin 
        N_pos_tot = np.sum( N_pos[ph_xq_c==iph_xq] )
        N_pos_tot_err = np.sqrt(np.sum( N_pos_err[ph_xq_c==iph_xq]**2 )) #standard error propagation for a sum

        N_neg_tot = np.sum( N_neg[ph_xq_c==iph_xq] )
        N_neg_tot_err = np.sqrt(np.sum( N_neg_err[ph_xq_c==iph_xq]**2 )) #standard error propagation for a sum

        #Calculate the total asymmetry from the combined +/- helicities integrated over th_xq_cm bins
        N = N_pos_tot + N_neg_tot
        Asy_tot = (N_pos_tot - N_neg_tot) / N
        Asy_err2 = (2*N_neg_tot/N**2)**2 * N_pos_tot_err**2 + (2*N_pos_tot/N**2)**2 * N_neg_tot_err**2
        Asy_tot_err = np.sqrt(Asy_err2)
        
        f.write("%.2f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f \n"  % (iph_xq,  N_pos_tot, N_pos_tot_err, N_neg_tot, N_neg_tot_err, Asy_tot, Asy_tot_err) )
    f.close()
    
    
def plot_asym(kin_group="", particle="", all_thCM=0):

    #Brief: Function to plot asymmetry per KIN- setting or group (combined KIN- settings)
    #       In this instance: kin_group can be:   KIN-1, KIN-2, etc.  | OR | group-1, group-2, etc.

    #Code Usage: example 1. plot_asym("KIN-1", "Kaon", 0)   :: This will plot the Kaon Asymmetry vs. ph_xq for every th_xqCM bin
    #            example 2. plot_asym("group-2", "Pion", 1) :: This will plot the Pion Asymmetry vs. ph_xq integrated over all th_xqCM bins

    
    thxq_val = [2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5]  #central values of in-plane theta_xq_cm angles [deg]


    #Set General fontsize parameters
    xlfs = 20   #x-label font-size
    ylfs = 20   #y-label font-size
    tfs  = 22   #title font-size
    tks  = 20   #general tick size (for both axes)
    
    #Plot a single asymmetry graph (integrated over all th_cm bins)
    if (all_thCM):
        fig = plt.figure(figsize=(16,8))
    #Plot a series of subplots to plot asymmetry per each th_cm bin
    else:
        fig, axs = plt.subplots(2, 4, sharex=True, sharey=True, figsize=(16,8))

        
    #Set strings for common plot titles
    if particle=="Kaon":
        title=r'Beam-Spin Asymmetry (%s): $^{1}H(e,e^{\prime}K^{+})\Lambda$' % (kin_group)
        title2=r'Beam-Spin Asymmetry (%s, $\int\theta^{\mathrm{cm}}_{K^{+},q}$): $^{1}H(e,e^{\prime}K^{+})\Lambda$' % (kin_group)
        xlabel=r'$\phi_{K^{+},q}$ [deg]'
        ylabel='Asymmetry'
    if particle=="Pion":
        title=r'Beam-Spin Asymmetry (%s): $^{1}H(e,e^{\prime}\pi^{+})n$' % (kin_group)
        title2=r'Beam-Spin Asymmetry (%s, $\int\theta^{\mathrm{cm}}_{\pi^{+},q}$): $^{1}H(e,e^{\prime}\pi^{+})n$' % (kin_group)
        xlabel=r'$\phi_{\pi^{+},q}$ [deg]'
        ylabel='Asymmetry'



    #plot asymmetry integrated over all th_xqCM bins
    if(all_thCM):
        #Read datafiles (extracted from 2D bin information of thxq_CM_vs_phxq)
        fname = './asymmetry_%s_%s_allthCM.txt'%(particle, kin_group) 
        asy_file = dfile(fname)

        ph_xq = np.array(asy_file['ph_xq'])
        Asy = np.array(asy_file['Asy'])
        Asy_err = np.array(asy_file['Asy_err'])
        
        B.plot_exp(ph_xq, Asy, Asy_err, marker='o', color='r', markerfacecolor='white')
        B.pl.xlabel(xlabel, fontsize=xlfs)
        B.pl.ylabel(ylabel, fontsize=ylfs)
        B.pl.title(title2, fontsize=tfs)
        B.pl.xticks(fontsize=tks)
        B.pl.yticks(fontsize=tks)
        B.pl.savefig('%s_asymmetry_%s_allthCM.pdf' % (particle, kin_group))

    #plot asymmetry binned in th_xqCM
    else:
        #Read datafiles (extracted from 2D bin information of thxq_CM_vs_phxq)
        fname = './asymmetry_%s_%s.txt'%(particle, kin_group) 
        asy_file = dfile(fname)

       
        #Set Common Title
        fig.text(0.5, 0.05, xlabel, ha='center', fontsize=xlfs)
        fig.text(0.06, 0.5, ylabel, va='center', rotation='vertical', fontsize=ylfs)
        plt.suptitle(title, fontsize=tfs)
        
        for i, ax in enumerate(axs.flatten()):

            #set tick params
            ax.tick_params(axis='both', labelsize=tks)
        
            #Read array values
            ith_xq = thxq_val[i]
            th_xqCM_c = np.array(asy_file['y0'])
            Asy = np.array(asy_file['Asy'])[th_xqCM_c==ith_xq]
            Asy_err = np.array(asy_file['Asy_err'])[th_xqCM_c==ith_xq]
            ph_xq_c = np.array(asy_file['x0'])[th_xqCM_c==ith_xq]
            #print(ith_xq, Asy)
            ax.errorbar(ph_xq_c, Asy, Asy_err, marker='o', markerfacecolor='white', color='r', ls='none')
            #set y-axis limit per subplot
            if particle=="Pion":
                ax.set_ylim(-0.25,0.25)
                ax.set_title(r'$<\theta^{\mathrm{cm}}_{\pi^{+},q}>$=%.1f $\pm 2.5^{\circ}$'%(ith_xq), fontsize=xlfs)

            if particle=="Kaon":
                ax.set_ylim(-0.5,0.5)
                ax.set_title(r'$<\theta^{\mathrm{cm}}_{K^{+},q}>$=%.1f $\pm 2.5^{\circ}$'%(ith_xq), fontsize=xlfs)

            
        #plt.show()
        plt.savefig('%s_asymmetry_%s.pdf'%(particle, kin_group))



    
def plot_report(kin_group=""):
        
    #Brief: Funtion to plot report output of each kinematic group, i.e.,
    #       tracking efficiencies (+/- helicities), live time, rates, etc.

    #Code Usage: plot_report("KIN-1")
    
    #Set General fontsize parameters
    xlfs = 16   #x-label font-size
    ylfs = 16   #y-label font-size
    tfs  = 20   #title font-size
    tks  = 20   #general tick size (for both axes)

    #Read report file
    fname_report = '../BeamPolAsymmetry/analysis_scripts/OUTPUT/%s/basic_report.txt'%(kin_group) 
    print(fname_report)
    f = dfile(fname_report)
    
    #read the headers
    run = np.array(f['Run'])

    #hadron arm (SHMS) tracking eff.
    htrk = np.array(f['hTrkEff'])
    htrk_err = np.array(f['hTrkEff_err'])
    
    #+ helicity
    htrk_pos = np.array(f['hTrkEff_pos'])
    htrk_pos_err = np.array(f['hTrkEff_err_pos'])

    #- helicity
    htrk_neg = np.array(f['hTrkEff_neg'])
    htrk_neg_err = np.array(f['hTrkEff_err_neg'])

    #electron arm (HMS) tracking eff.
    etrk = np.array(f['pTrkEff'])
    etrk_err = np.array(f['pTrkEff_err'])
    
    #+ helicity
    etrk_pos = np.array(f['pTrkEff_pos'])
    etrk_pos_err = np.array(f['pTrkEff_err_pos'])

    #- helicity
    etrk_neg = np.array(f['pTrkEff_neg'])
    etrk_neg_err = np.array(f['pTrkEff_err_neg'])

    # total live time
    tLT = np.array(f['tLT'])
    cpuLT = np.array(f['cpuLT'])

    # HMS / SHMS  rates
    hms_rates = np.array(f['hms_Ps3_elreal_rate'])
    shms_rates = np.array(f['shms_Ps1_3of4_rate'])
    # coin rates
    coin_rates = np.array(f['coin_Ps5_rate'])
    
    fig, axes = plt.subplots(2,2, figsize=(14, 12))

    #PLOT SHMS TRK EFF VS. RUN
    l1 = axes[0, 0].errorbar(run, htrk, htrk_err, marker='o', markerfacecolor='white', color='k', ls='none', label='ALL')
    l2 = axes[0, 0].errorbar(run, htrk_pos, htrk_pos_err, marker='o', markerfacecolor='white', color='b', ls='none', label='+ hel')
    l3 = axes[0, 0].errorbar(run, htrk_neg, htrk_neg_err, marker='o', markerfacecolor='white', color='r', ls='none', label='= neg')
    axes[0,0].legend([l1, l2, l3], ['ALL', 'hel (pos)', 'hel (neg)' ])
    axes[0,0].set_ylim(htrk.min()-0.01, htrk.max()+0.01)
    axes[0,0].set_title('Hadron (SHMS) Tracking Efficiencies', fontsize=tfs)
    axes[0,0].set_xlabel('Run', fontsize=xlfs)
    axes[0,0].set_ylabel(r'Hadron Tracking Efficiency', fontsize=ylfs)

    #PLOT HMS TRK EFF. VS. RUN
    l1 = axes[0, 1].errorbar(run, etrk, etrk_err, marker='o', markerfacecolor='white', color='k', ls='none', label='ALL')
    l2 = axes[0, 1].errorbar(run, etrk_pos, etrk_pos_err, marker='o', markerfacecolor='white', color='b', ls='none', label='+ hel')
    l3 = axes[0, 1].errorbar(run, etrk_neg, etrk_neg_err, marker='o', markerfacecolor='white', color='r', ls='none', label='= neg')
    axes[0,1].legend([l1, l2, l3], ['ALL', 'hel (pos)', 'hel (neg)' ])
    axes[0,1].set_ylim(etrk.min()-0.01, etrk.max()+0.01)
    axes[0,1].set_title('Electron (HMS) Tracking Efficiencies', fontsize=tfs)
    axes[0,1].set_xlabel('Run', fontsize=xlfs)
    axes[0,1].set_ylabel(r'Electron Tracking Efficiency', fontsize=ylfs)
    
    #PLOT TOTAL LIVE TIME VS. RUN
    l1, = axes[1, 0].plot(run, tLT,  marker='o', markerfacecolor='white', color='k', ls='none', label='EDTM Live Time')
    l2, = axes[1, 0].plot(run, cpuLT,  marker='o', markerfacecolor='white', color='r', ls='none', label='CPU Live Time')   
    axes[1,0].legend((l1,l2), ('Total (EDTM) Live Time', 'Computer Live Time'))
    axes[1,0].set_ylim(tLT.min()-0.01, cpuLT.max()+0.01)
    axes[1,0].set_title('Live Time', fontsize=tfs)
    axes[1,0].set_xlabel('Run', fontsize=xlfs)
    axes[1,0].set_ylabel(r'Live Time', fontsize=ylfs)

    
    #PLOT RATES VS. RUN
    l1, = axes[1, 1].plot(run, shms_rates,  marker='o', markerfacecolor='white', color='r', ls='none')   
    l2, = axes[1, 1].plot(run, hms_rates,  marker='o', markerfacecolor='white', color='b', ls='none')
    l3, = axes[1, 1].plot(run, coin_rates,  marker='o', markerfacecolor='white', color='k', ls='none')
    axes[1,1].legend((l1,l2,l3), ('SHMS Rates', 'HMS Rates', 'Coincidence Rates'), loc='best')
    axes[1,1].set_ylim(coin_rates.min()-10*coin_rates.min(), shms_rates.max()+0.25*shms_rates.max())
    axes[1,1].set_title('Trigger Rates', fontsize=tfs)
    axes[1,1].set_xlabel('Run', fontsize=xlfs)
    axes[1,1].set_ylabel(r'Trigger Rates [kHz]', fontsize=ylfs) 
                     
    #B.pl.show()

    B.pl.savefig('%s_report.pdf' % (kin_group))
    


    
def main():
    print("Plotting Asym")

    kin_set = "KIN-16"

    
    # calculate asymmetry for a specified kin. setting
    # (output .txt file with numerical asymmetry is created)
    calc_asym("Kaon", kin_set)
    calc_asym("Pion", kin_set)

    # integrate over all th_xq_cm bins for specified kinematic
    # (output .txt file with integrated numerical asymmetry is created)
    combine_thcm_bins(particle="Kaon", group=kin_set)
    combine_thcm_bins(particle="Pion", group=kin_set)

    # plot asymmetry for specified kin. setting, particle, binned (0) or integrated (1) over th_xq_cm bins
    plot_asym(kin_set, "Kaon", 0)
    plot_asym(kin_set, "Kaon", 1)

    plot_asym(kin_set, "Pion", 0)
    plot_asym(kin_set, "Pion", 1)
   
    plot_report(kin_set)
    
    
     
    '''
    #combined multiple Kaon LT kinematics with the same HMS setting (and varying SHMS angles)
    calc_asym_combined(np.array(["KIN-1", "KIN-2", "KIN-3"]),              "Kaon", "group-1")
    calc_asym_combined(np.array(["KIN-4", "KIN-5", "KIN-6"]),              "Kaon", "group-2")
    calc_asym_combined(np.array(["KIN-7", "KIN-8", "KIN-9"]),              "Kaon", "group-3")
    calc_asym_combined(np.array(["KIN-10", "KIN-11", "KIN-12"]),           "Kaon", "group-4")
    calc_asym_combined(np.array(["KIN-13", "KIN-14", "KIN-15", "KIN-16"]), "Kaon", "group-5")

    calc_asym_combined(np.array(["KIN-1", "KIN-2", "KIN-3"]),              "Pion", "group-1")
    calc_asym_combined(np.array(["KIN-4", "KIN-5", "KIN-6"]),              "Pion", "group-2")
    calc_asym_combined(np.array(["KIN-7", "KIN-8", "KIN-9"]),              "Pion", "group-3")
    calc_asym_combined(np.array(["KIN-10", "KIN-11", "KIN-12"]),           "Pion", "group-4")
    calc_asym_combined(np.array(["KIN-13", "KIN-14", "KIN-15", "KIN-16"]), "Pion", "group-5")
    
    
    #calculate asymmetries integrated over all th_cm_bins
    combine_thcm_bins("Kaon", "group-1")
    combine_thcm_bins("Kaon", "group-2")
    combine_thcm_bins("Kaon", "group-3")
    combine_thcm_bins("Kaon", "group-4")
    combine_thcm_bins("Kaon", "group-5")
    
    combine_thcm_bins("Pion", "group-1")
    combine_thcm_bins("Pion", "group-2")
    combine_thcm_bins("Pion", "group-3")
    combine_thcm_bins("Pion", "group-4")
    combine_thcm_bins("Pion", "group-5")   
    '''
    
    #plot asymmetries binned in th_cm
    #plot_asym("group-1", "Kaon", 0)
    #plot_asym("group-2", "Kaon", 0)
    #plot_asym("group-3", "Kaon", 0)
    #plot_asym("group-4", "Kaon", 0)
    #plot_asym("group-5", "Kaon", 0)
    
    
    #plot_asym("group-1", "Pion", 0)
    #plot_asym("group-2", "Pion", 0)
    #plot_asym("group-3", "Pion", 0)
    #plot_asym("group-4", "Pion", 0)
    #plot_asym("group-5", "Pion", 0)
    
    #plot integrated asymmetries
    #plot_asym("group-1", "Kaon", 1)
    #plot_asym("group-2", "Kaon", 1)
    #plot_asym("group-3", "Kaon", 1)
    #plot_asym("group-4", "Kaon", 1)
    #plot_asym("group-5", "Kaon", 1)

    
    #plot_asym("group-1", "Pion", 1)
    #plot_asym("group-2", "Pion", 1)
    #plot_asym("group-3", "Pion", 1)
    #plot_asym("group-4", "Pion", 1)
    #plot_asym("group-5", "Pion", 1)
    
    
    
if __name__ == "__main__":
    main()
