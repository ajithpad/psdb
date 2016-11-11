# script to analyze batch data : firing rates, synchronicity and regularity

import numpy as np
from numpy import array
import pylab as pl
import scipy.io as sio
import os
import common.analyze_data_psd as adata
#reload(adata)
import itertools
import socket
import sys
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import LogNorm
import pdb

labelsize = 13
ticksize = 13
def get_sim_3d(match_pars):
    ''' get simulation object corresponding to (g,eta,inh2_ratio)
    match_pars is a dictionary, e.g. {'g':5.6,'eta':7.6, 'inh2_ratio':0.5}
    CAUTION: check order of keys'''
    #batch_pars ={'pg_rate_exc':exc_range,'num_spikes_burst':spb_range,'inh2_ratio':ratio_range,'g':g_range}
    batch_pars ={'pg_rate_exc':exc_range,'num_spikes_burst':spb_range,'inh2_ratio':ratio_range,'g':g_range,'extra_inh':inh_range}
    comb_pars = [[value for (key, value) in zip(batch_pars, values)] for values in itertools.product(*batch_pars.values())] 
    sim_id = -1
    keys_sort = batch_pars.keys() 
    match_pars_lis = [match_pars[k] for k in keys_sort]
    for ii,pars in enumerate(comb_pars):
	if match_pars_lis==np.array(pars).round(2).tolist():
	    sim_id = ii	  	    
    if sim_id>-1:
	return sim_id
    else:
	print 'No %s found'%match_pars.keys()
	return comb_pars,batch_pars.keys()


# simulation name

#prefix0 = 'mps_izh_tc_2'
prefix0 = 'mps_izh_rs_fs_ib-2'
prefix1 = sys.argv[1]
prefix0 = sys.argv[2]
prefix0 = prefix1+'_'+prefix0

col0 = np.array([51,102,204])/255.

# rate, ff and CV path
path1 = '/media/DISK_B/results/levcomp/res_analysis/'
path2 = '/media/DISK_B/scripts_backup/levcomp/batch/'


contents = os.listdir(path1)
max_sim_id = 1319#1583#1799#2474#2639#1055#5279
rate = np.load(path1+prefix0+'_rate_0_%s.npy'%max_sim_id)
ff = np.load(path1+prefix0+'_ff_0_%s.npy'%max_sim_id)
spec = np.load(path1+prefix0+'_spec_0_%s.npy'%max_sim_id)
pars = eval(open(path2+prefix1+'.log').read())
nr = 0
nr2 = 0
rate2 = np.zeros((len(rate),))
ff2 = np.zeros((len(ff),))
spec2 = np.zeros((len(spec),))

for ii in np.arange(len(rate2)):   
	rate2[nr] = rate[ii]
	ff2[nr] = ff[ii]
	spec2[nr] = spec[ii]
	nr = nr +1
        
rate_mat = np.zeros((len(pars['pg_rate_exc']),len(pars['g']),len(pars['inh2_ratio'])))
ff_mat = np.zeros((len(pars['pg_rate_exc']),len(pars['g']),len(pars['inh2_ratio'])))
spec_mat = np.zeros((len(pars['pg_rate_exc']),len(pars['g']),len(pars['inh2_ratio'])))

# CAUTION: g and eta are reversed in dictionary
nr = 0
#nr2 = 0

exc_range = pars['pg_rate_exc']
inh_range = pars['extra_inh']
g_range = pars['g']
ratio_range = pars['inh2_ratio']
spb_range = pars['num_spikes_burst']
nr2 = 0
for id1,exc_val in enumerate(exc_range):
      for id2,g_val in enumerate(g_range):
	  for id3,inh_ratio in enumerate(ratio_range):
	    sim_id = get_sim_3d({'pg_rate_exc':exc_val,'extra_inh':9000.,'inh2_ratio':inh_ratio,'g':g_val,'num_spikes_burst':5.})
	    rate_mat[id1,id2,id3] = rate2[sim_id]
	    ff_mat[id1,id2,id3] = ff2[sim_id]
	    spec_mat[id1,id2,id3] = spec2[sim_id]
	    nr = nr+1

masked_rate = np.ma.array(rate_mat, mask=np.isnan(rate_mat))
masked_ff = np.ma.array(ff_mat, mask=np.isnan(ff_mat))
masked_spec = np.ma.array(spec_mat, mask=np.isnan(spec_mat))

exc_range_l = exc_range.tolist()
g_range_l = g_range.tolist()

exc_range_l.append(351000)
g_range_l.append(11.0)
#try:	
ratio_range = ratio_range[[0,5,10]]
fig = pl.figure(5, figsize = (8,12))
for id2,ratio_val in enumerate(ratio_range):
    print ratio_val
    #fig = pl.figure(4, figsize = (8,12))
    #pl.subplot(1,3,id2+1) 
    #ax1 = pl.pcolor(np.array(g_range_l),np.array(exc_range_l),masked_spec[:,:,id2*5])
    #pl.clim(0.,1.0)
    #pl.gca().set_xticks(g_range_l[::2])
    #pl.gca().set_yticks(exc_range_l[::2])
    #pl.xlabel('g',fontsize = 20)
    #pl.ylabel('External input', fontsize = 20)       
    
    #cax = fig.add_axes([0.92,0.25,0.015,0.45])
    #cbar = pl.colorbar(ax1,cax = cax)
    #cbar.set_label('Spectral entropy',position= (2.9,0.5),fontsize = 20, rotation = 'vertical')
##      pl.savefig('corr_coef.pdf', format = 'pdf')

    
    ax1 = fig.add_axes([0.125 + id2*0.25,0.8,0.2,0.133]) 
    pl.pcolor(np.array(g_range_l),np.array(exc_range_l),masked_rate[:,:,id2*5],cmap = cm.jet)
    pl.clim(0.,10.0)
    pl.gca().set_xticks(g_range[::2])
    for tl in ax1.get_xticklabels():
	tl.set_fontsize(ticksize)
    ax1.set_xlim(g_range_l[0],g_range_l[-1])
    ax1.set_ylim(exc_range_l[0],exc_range_l[-1])
    if id2 == 0:
      pl.gca().set_yticks(exc_range_l[::2])
      for tl in ax1.get_yticklabels():
	tl.set_fontsize(ticksize)
	ax1.yaxis.get_major_formatter().set_powerlimits((0,1))
	ax1.set_ylabel('External input rate',fontsize = labelsize)    
    else:
      pl.gca().set_yticks([])
    pl.xlabel('g',fontsize = labelsize)      

    cr = np.arange(0.,10.1,2.5)
    cax = fig.add_axes([0.875,0.825,0.02,0.1])
    cbar = pl.colorbar(orientation = 'vertical',cax = cax, cmap = cm.jet)
    cbar.set_label('Firing rate (Hz)',position= (1.3,0.5),fontsize = labelsize, rotation = 'vertical')
    cbar.set_ticks(cr[::2])
    cticks = pl.getp(cbar.ax.axes,'yticklabels')
    pl.setp(cticks,fontsize = ticksize)
    
    #pl.savefig('firing_rate.pdf', format = 'pdf')
    ax2 = fig.add_axes([0.125 + id2*0.25,0.6,0.2,0.133])
    pl.pcolor(np.array(g_range_l),np.array(exc_range_l),masked_ff[:,:,id2*5],cmap = cm.jet)
    pl.clim(0.,2.5)
    pl.gca().set_xticks(g_range[::2])
    for tl in ax2.get_xticklabels():
	tl.set_fontsize(ticksize)
    ax2.set_xlim(g_range_l[0],g_range_l[-1])
    ax2.set_ylim(exc_range_l[0],exc_range_l[-1])
    if id2 == 0:
      pl.gca().set_yticks(exc_range_l[::2])
      for tl in ax1.get_yticklabels():
	tl.set_fontsize(ticksize)
	ax2.yaxis.get_major_formatter().set_powerlimits((0,1))
	ax2.set_ylabel('External input rate',fontsize = labelsize)    
      #ax1.set_yticks(fontsize = 45)
    else:
      pl.gca().set_yticks([])
    pl.xlabel('g',fontsize = labelsize)      

    cr = np.arange(0.,2.6,1.)
    cax = fig.add_axes([0.875,0.625,0.02,0.1])
    cbar = pl.colorbar(orientation = 'vertical',cax = cax, cmap = cm.jet)
    cbar.set_label('Fano Factor(Hz)',position= (1.3,0.5),fontsize = labelsize, rotation = 'vertical')
    cbar.set_ticks(cr[::2])
    cticks = pl.getp(cbar.ax.axes,'yticklabels')
    pl.setp(cticks,fontsize = ticksize)
    #pl.savefig('fano_factor.pdf', format = 'pdf')
pl.show()

