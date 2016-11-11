
# script to analyze batch data : firing rates, synchronicity and regularity

# NOTE: run analyze_batch.py mps_izh_rs_fs_ch-1 rs_1 0 1539 regular_spiking none


import numpy as np
import pylab as pl
import scipy.io as sio
import os
import re
import socket
import common.analyze_data_psd as adata
reload(adata)
import sys
from numpy import array
import cPickle as cp

# simulation id
#prefix0 = 'mps_izh_fs_4'
#prefix0 = 'mps_izh_rs_4'
#prefix0 = 'mps_izh_ch_4'
#prefix0 = 'mps_izh_ib_4'
#prefix0 = 'mps_izh_res_4'
#prefix0 = 'mps_izh_tc_4'



prefix0 = sys.argv[1]
prefix1 = sys.argv[2]
min_nr = int(sys.argv[3])
max_nr = int(sys.argv[4])
#NOTE : important to check if clause in order to exclude some names; use ' ' for no exclusion
pop_name = sys.argv[5]
exclude_files = sys.argv[6]	


# path to load results
host = socket.getfqdn()
#if 'lusitania' in host:
  #path1 = '/home/vlachos/bcfgrid2/results/levcomp/'
  #if not os.path.exists('/home/vlachos/bcfgrid2/results/levcomp/res_analysis/'):
    #os.makedirs('/home/vlachos/bcfgrid2/results/levcomp/res_analysis/')
  #path2 = '/home/vlachos/bcfgrid2/results/levcomp/res_analysis/'
  #path3 = '/home/vlachos/bcfgrid_root/bcfgrid/data/vlachos/scripts/levcomp/batch/'
  #print_time = True
#elif 'bcfgrid' in host:
  #if 'vlachos' in os.getcwd():
    #path1 = '/bcfgrid/data/vlachos/levcomp/results/'
    #if not os.path.exists('/bcfgrid/data/vlachos/levcomp/results/res_analysis/'):
      #os.makedirs('/bcfgrid/data/vlachos/levcomp/results/res_analysis/')
    #path2 = '/bcfgrid/data/vlachos/levcomp/results/res_analysis/'	
  #elif 'padmanabhan' in os.getcwd():
    #path1 = '/bcfgrid/data/padmanabhan/results/levcomp/' 
    #path2 = '/bcfgrid/data/padmanabhan/results/levcomp/res_analysis/'
  #print_time = False
#elif 'lew' in host:
  #path1 = '/bcfgrid/data/padmanabhan/results/levcomp/' 
  #path2 = '/bcfgrid/data/padmanabhan/results/levcomp/res_analysis/'
  #print_time = True

#path1 = '/media/Transcend/ffcon_chcon-10/'
#path1 = '/media/Transcend/simdata_3/' 
path1 = '/media/DISK_C/results/levcomp/'
path2 = '/media/DISK_B/results/levcomp/res_analysis/'
path3 = '/media/DISK_B/scripts_backup/levcomp/batch/'

## spikes and vm path
#path1 = '/home/vlachos/bcfgrid_root/bcfgrid/data/padmanabhan/results/levcomp/' 
## rate,cv,ff, path
#path2 = '/home/vlachos/bcfgrid_root/bcfgrid/data/vlachos/results/levcomp/res_analysis/'
## batch-scripts path
#path3 = '/home/vlachos/bcfgrid_root/bcfgrid/data/padmanabhan/scripts/levcomp/batch/'


contents = os.listdir(path1)
prefix = []
#for fname in contents:
    #if 'spikes.npy' in fname:	
	#idx = fname.find('_spikes')
	#prefix.append(fname[:idx])

for fname in contents:
    #if ('spikes.npy' in fname) and (('bursting' in fname) or ('spiking' in fname) or ('adaptation' in fname)):	
    if ('spikes.npy' in fname) and (prefix0 in fname) and (exclude_files not in fname):
	idx = fname.find('_spikes')
	prefix.append(fname[:idx])
	print idx,fname[:idx]

file_nr = np.zeros((len(prefix),))
for ii,string in enumerate(prefix):
    file_nr[ii] = re.findall('[0-9]{1,4}',string)[-1]
order = np.argsort(file_nr)
prefix = np.array(prefix)[order]

# NOTE: get the ids of missing files, e.g. killed due to walltime
pars = eval(open(path3+prefix0+'.log').read())
g_range = pars['g']
exc_range = pars['pg_rate_exc']
inh2_ratio = pars['inh2_ratio']
inh_range = pars['extra_inh']
spb_range = pars['num_spikes_burst']
#spa_pg_rate_range = pars['spa_pg_rate']
#spa_epsilon_range = pars['spa_epsilon']
#inh2_ratio_range = pars['inh2_ratio']

#xc2_ratio = pars['xc2_ratio']
tot_files = len(exc_range) * len(g_range) * len(inh2_ratio)* len(spb_range) * len(inh_range)
#tot_files = len(spa_pg_rate_range) * len(spa_epsilon_range) * len(inh2_ratio_range)
files0 = np.zeros((tot_files,))
files0[array(file_nr,int)] =1
missing_ids = pl.find(files0==0)
#missing_ids = [1032,1033]

mrate_tot = []
ff_tot = []
cv_tot = []
cv_kl_tot = []
cc_tot = []
osc_ind_tot = []
burst_ind_tot = []
inh_ratio_tot = []
comp_fil_ff_tot = []
max_freq_tot = []
power_tot = []
inh_rate_tot = []
pg_rate_tot = []
freq1_tot = []
freq2_tot = []
time1_tot = []
time2_tot = []
spec_tot = []
peak_val_tot = []


prefix0 = prefix0 + '_' + prefix1


for nr,pfx in enumerate(prefix[min_nr:max_nr+1]):
    print nr,pfx
    res = adata.analyze_data(pfx,path1)
    #pl.figure()
    #res.plot_rate_hist()
    #res.plot_raster()
    #pl.subplot(211)
    #pl.title(pfx)
    #res.load_vm()
    #vm = res.events['vm']
    #pl.figure()
    #pl.plot(vm[:,2])
    nmax = 1000
    #bin_size = 1 # (ms) for cv_kl
    #nr_peaks = 5
    mrate = res.comp_mean_rate(pop_id=pop_name,nmax=nmax,time_range = [res.pars['T_wup'],res.pars['T_wup'] + res.pars['T_sim']])	# FIRING RATES
    #inh_rate = res.inh_mean_rate(time_range = [res.pars['T_wup'],res.pars['T_wup'] + res.pars['T_sim']])	# FIRING RATES
    #pg_rate = res.coll_pg_rate()
    ff = res.comp_ff(kernel_w=2.,pop_id=pop_name,nmax = nmax,time_range = [res.pars['T_wup'],res.pars['T_wup'] + res.pars['T_sim']]) 	# SYNCHRONY; comp_cc maybe more accurate but takes more time
    #cv,dummy = res.comp_cv(nmax=nmax,pop_id=pop_name,time_range = [res.pars['T_wup'],res.pars['T_wup'] + res.pars['T_sim']])# REGULARITY  \
    #dummy,dummy,dummy,freq = res.psd(nmax=nmax,pop_id=pop_name,time_range = [res.pars['T_wup'],res.pars['T_wup'] + res.pars['T_sim']], bin_w = 10.)    
    #freq1,dummy,dummy,dummy = res.comp_osc_ind2(pop_id=pop_name,time_range = [res.pars['T_wup'],res.pars['T_wup'] + res.pars['T_sim']], nr_peaks = 2)
    #dummy,freq2,dummy,dummy = res.comp_osc_ind2(pop_id=pop_name,time_range = [res.pars['T_wup'],res.pars['T_wup'] + res.pars['T_sim']], nr_peaks = 2)
    #dummy,dummy,time1,dummy = res.comp_osc_ind2(pop_id=pop_name,time_range = [res.pars['T_wup'],res.pars['T_wup'] + res.pars['T_sim']], nr_peaks = 2)
    #dummy,dummy,dummy,time2 = res.comp_osc_ind2(pop_id=pop_name,time_range = [res.pars['T_wup'],res.pars['T_wup'] + res.pars['T_sim']], nr_peaks = 2)
    #power,dummy,dummy,dummy = res.psd(nmax=nmax,pop_id=pop_name,time_range = [res.pars['T_wup'],res.pars['T_wup'] + res.pars['T_sim']])
    #dummy,dummy,fil_ff = res.comp_fil_ff(nmax=nmax,pop_id=pop_name,time_range = [res.pars['T_wup'],res.pars['T_wup']+res.pars['T_sim']])
    #cv_kl = res.comp_cv_kl(pop_id=pop_name,bin_size=bin_size,nmax = nmax,time_range = [res.pars['T_wup'],(res.pars['T_sim']+res.pars['T_wup'])])	# REGULARITY KL
    #dummy,cc,dummy = res.comp_cc(nmax=nmax,pop_id=pop_name,time_range = [res.pars['T_wup'],res.pars['T_wup']+res.pars['T_sim']])
    #osc_ind,dummy,dummy,dummy,dummy = res.comp_osc_ind2(pop_id=pop_name)		# OSCILLATION INDEX
    #burst_ind,dummy,dummy,dummy = res.compute_burst_index(pop_id=pop_name,nmax=nmax)	#only for inh2
    spec = res.spec_entropy(time_range = [res.pars['T_wup']+200.,res.pars['T_wup']+res.pars['T_sim']],freq_range = [0.,100.])
#    peak_val,dummy = res.find_peaks(thr = 15.,time_range = [res.pars['st_val'],res.pars['st_val']+50.])
    #inh_ratio = res.pars['inh2_ratio']
    mrate_tot.append(mrate)
    #inh_rate_tot.append(inh_rate)
    ff_tot.append(ff)
    #pg_rate_tot.append(pg_rate)
    #cv_tot.append(cv)
    #max_freq_tot.append(freq)
    #freq1_tot.append(freq1)
    #freq2_tot.append(freq2)
    #time1_tot.append(time1)
    #time2_tot.append(time2)
    spec_tot.append(spec)
    #power_tot.append(power)
    #comp_fil_ff_tot.append(fil_ff)
    #inh_ratio_tot.append(inh_ratio)
    #cv_kl_tot.append(cv_kl)
    #cc_tot.append(cc)
    #osc_ind_tot.append(osc_ind)
    #burst_ind_tot.append(burst_ind)
#    peak_val_tot.append(peak_val)        
    
    np.save(path2+prefix0+'_rate_%d_%d'%(min_nr,max_nr),mrate_tot)
    #np.save(path2+prefix0+'_inh_rate_%d_%d'%(min_nr,max_nr),inh_rate_tot)
    np.save(path2+prefix0+'_ff_%d_%d'%(min_nr,max_nr),ff_tot)
    #np.save(path2+prefix0+'_cv_%d_%d'%(min_nr,max_nr),cv_tot)
    #np.save(path2+prefix0+'_freq_%d_%d'%(min_nr,max_nr),max_freq_tot)
    #np.save(path2+prefix0+'_freq1_%d_%d'%(min_nr,max_nr),freq1_tot)
    #np.save(path2+prefix0+'_freq2_%d_%d'%(min_nr,max_nr),freq2_tot)
    #np.save(path2+prefix0+'_time1_%d_%d'%(min_nr,max_nr),time1_tot)
    #np.save(path2+prefix0+'_time2_%d_%d'%(min_nr,max_nr),time2_tot)
    np.save(path2+prefix0+'_spec_%d_%d'%(min_nr,max_nr),spec_tot)
    #np.save(path2+prefix0+'_power_%d_%d'%(min_nr,max_nr),power_tot)
    #np.save(path2+prefix0+'_pg_rate_%d_%d'%(min_nr,max_nr),pg_rate_tot)
    #np.save(path2+prefix0+'_fil_ff_%d_%d'%(min_nr,max_nr),comp_fil_ff_tot)
    #np.save(path2+prefix0+'_inh_ratio_%d_%d'%(min_nr,max_nr),inh_ratio_tot)
    #np.save(path2+prefix0+'_cv_kl_%d_%d'%(min_nr,max_nr),cv_kl_tot)
    #np.save(path2+prefix0+'_cc_%d_%d'%(min_nr,max_nr),cc_tot)
    #np.save(path2+prefix0+'_osc_ind_%d_%d'%(min_nr,max_nr),osc_ind_tot)
    #np.save(path2+prefix0+'_burst_ind_%d_%d'%(min_nr,max_nr),burst_ind_tot)
#    np.save(path2+prefix0+'_peak_val_%d_%d'%(min_nr,max_nr),peak_val_tot)

# NOTE: fill the missing files with nan
for ii in missing_ids:
    mrate_tot.insert(ii,np.nan)
    #inh_rate_tot.insert(ii,np.nan)
    ff_tot.insert(ii,np.nan)
    #cv_tot.insert(ii,np.nan)
    #max_freq_tot.insert(ii,np.nan)
    #power_tot.insert(ii,np.nan)
    spec_tot.insert(ii,np.nan)
    #pg_rate_tot.insert(ii,np.nan)
    #cv_kl_tot.insert(ii,np.nan)
    #osc_ind_tot.insert(ii,np.nan)
    #inh_ratio_tot.insert(ii,np.nan)
##pl.show()

np.save(path2+prefix0+'_rate_%d_%d'%(min_nr,max_nr),mrate_tot)
#np.save(path2+prefix0+'_inh_rate_%d_%d'%(min_nr,max_nr),inh_rate_tot)
np.save(path2+prefix0+'_ff_%d_%d'%(min_nr,max_nr),ff_tot)
#np.save(path2+prefix0+'_cv_%d_%d'%(min_nr,max_nr),cv_tot)
#np.save(path2+prefix0+'_freq_%d_%d'%(min_nr,max_nr),max_freq_tot)
#np.save(path2+prefix0+'_freq1_%d_%d'%(min_nr,max_nr),freq1_tot)
#np.save(path2+prefix0+'_freq2_%d_%d'%(min_nr,max_nr),freq2_tot)
#np.save(path2+prefix0+'_time1_%d_%d'%(min_nr,max_nr),time1_tot)
#np.save(path2+prefix0+'_time2_%d_%d'%(min_nr,max_nr),time2_tot)
np.save(path2+prefix0+'_spec_%d_%d'%(min_nr,max_nr),spec_tot)
#np.save(path2+prefix0+'_power_%d_%d'%(min_nr,max_nr),power_tot)
#np.save(path2+prefix0+'_pg_rate_%d_%d'%(min_nr,max_nr),pg_rate_tot)
#np.save(path2+prefix0+'_fil_ff_%d_%d'%(min_nr,max_nr),comp_fil_ff_tot)
#np.save(path2+prefix0+'_inh_ratio_%d_%d'%(min_nr,max_nr),inh_ratio_tot)
#np.save(path2+prefix0+'_cv_kl_%d_%d'%(min_nr,max_nr),cv_kl_tot)
#np.save(path2+prefix0+'_cc_%d_%d'%(min_nr,max_nr),cc_tot)
#np.save(path2+prefix0+'_osc_ind_%d_%d'%(min_nr,max_nr),osc_ind_tot)
#np.save(path2+prefix0+'_burst_ind_%d_%d'%(min_nr,max_nr),burst_ind_tot)
#np.save(path2+prefix0+'_peak_val_%d_%d'%(min_nr,max_nr),peak_val_tot)
