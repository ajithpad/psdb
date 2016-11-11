# Code to plot the changes in STN and GPe bursting ratio against change in population activity

import numpy as np
import numpy as np
import itertools
import sys
sys.path.append('/users/padmanabhan/simulations/basal/scripts/common/')
import analyze_data_psd as adata
#import common.ff_analyze_data as adata
reload(adata)
import sys
from numpy import *
reload(adata)
import pdb
import pylab as pl
import matplotlib.cm as cm



p1 = []
lab_names = []
legend_dict = dict()
def get_sim_3d(match_pars):
    ''' get simulation object corresponding to (g,eta,inh2_ratio)
    match_pars is a dictionary, e.g. {'g':5.6,'eta':7.6, 'inh2_ratio':0.5}
    CAUTION: check order of keys'''  
    batch_pars ={'pg_rate_exc':ext_rate,'extra_inh':add_inh_rate,'exc2_ratio':exc_range,'inh2_ratio':inh_range }
    #batch_pars ={'pg_rate_exc':ext_rate,'pg_rate_inh':add_inh_rate,'exc2_ratio':exc_range,'inh2_ratio':inh_range }
    comb_pars = [[value for (key, value) in zip(batch_pars, values)] for values in itertools.product(*batch_pars.values())] 
    sim_id = -1
    keys_sort = batch_pars.keys() 
    match_pars_lis = [match_pars[k] for k in keys_sort]
    for ii,pars in enumerate(comb_pars):
	if match_pars_lis==np.array(pars).round(2).tolist():
	    sim_id = ii	  	
	    #print pars
    if sim_id>-1:
	res = adata.analyze_data(prefix0+'_%d'%sim_id,path1)
	#print 'sim_id is: ',sim_id
	return res
    else:
	print 'No %s found'%match_pars.keys()
	return comb_pars,batch_pars.keys()
	
	
prefix1 = 'small_extern_input-2'
path1 = '/data/padmanabhan/results/basal/'
prefix0 = prefix1
log_path = '/data/padmanabhan/scripts/basal/batch/'
log_file = open(log_path + prefix1 +'.log')
log_values = eval(log_file.read())
exc_range = log_values['exc2_ratio']
inh_range = log_values['inh2_ratio']
ext_rate = log_values['pg_rate_exc']
#add_inh_rate = log_values['extra_inh']
add_inh_rate = log_values['extra_inh']


#####Go through the different values of busrting ratio of GPe for the given external input values and calculate the spectral entropy and firing rate
inh_exc_stn_fr = []
inh_exc_stn_ff = []
inh_exc_gpe_fr = []
inh_exc_gpe_ff = []
inh_exc_gpe_spec = []
inh_exc_stn_spec = []
for inh_val in inh_range:
  exc_stn_fr = []
  exc_stn_ff = []
  exc_gpe_fr = []
  exc_gpe_ff = []
  exc_gpe_spec = []
  exc_stn_spec = []
  for exc_val in exc_range:
    ad1 = get_sim_3d({'pg_rate_exc':ext_rate[2],'extra_inh':add_inh_rate[2],'inh2_ratio':inh_val,'exc2_ratio':exc_val})
    stn_fr_val = ad1.comp_mean_rate(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_exc')
    stn_ff_val = ad1.comp_ff(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_exc', kernel_w= 2.)
    exc_stn_fr.append(stn_fr_val)
    exc_stn_ff.append(stn_ff_val)
    gpe_fr_val = ad1.comp_mean_rate(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh')
    gpe_ff_val = ad1.comp_ff(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh', kernel_w= 2.)
    exc_gpe_fr.append(gpe_fr_val)
    exc_gpe_ff.append(gpe_ff_val)
    xx,xx,xx,peak_freq_gpe, peak_pow_gpe = ad1.psd(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh')
    gpe_spec_val =  ad1.spec_entropy(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh', freq_range = [0.,100])
    exc_gpe_spec.append(gpe_spec_val)
    xx,xx,xx,peak_freq_stn, peak_pow_stn = ad1.psd(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_exc')
    stn_spec_val = ad1.spec_entropy(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_exc', freq_range = [0.,100])
    exc_stn_spec.append(stn_spec_val)
  inh_exc_stn_fr.append(exc_stn_fr)
  inh_exc_stn_ff.append(exc_stn_ff)
  inh_exc_gpe_fr.append(exc_gpe_fr)
  inh_exc_gpe_ff.append(exc_gpe_ff)
  inh_exc_gpe_spec.append(exc_gpe_spec)
  inh_exc_stn_spec.append(exc_stn_spec)

fig = pl.figure(1,(10,15))
ax1 = fig.add_subplot(321)
ax2 = fig.add_subplot(322)
ax3 = fig.add_subplot(323)
ax4 = fig.add_subplot(324)
ax5 = fig.add_subplot(325)
ax6 = fig.add_subplot(326)
for inh_count, inh_val in enumerate(inh_range):
  ax1.plot(exc_range, inh_exc_stn_fr[inh_count], color = cm.Reds_r(inh_count*30),lw = 3.,label = str(inh_val))
  ax3.plot(exc_range, inh_exc_stn_ff[inh_count], color = cm.Greens_r(inh_count*30),lw = 3.,label = str(inh_val))
  ax2.plot(exc_range, inh_exc_gpe_fr[inh_count], color = cm.Reds_r(inh_count*30),lw = 3.)
  ax4.plot(exc_range, inh_exc_gpe_ff[inh_count], color = cm.Greens_r(inh_count*30),lw = 3.)
  ax5.plot(exc_range, inh_exc_stn_spec[inh_count], color = cm.Blues_r(inh_count*30),lw = 3.,label = str(inh_val))
  ax6.plot(exc_range, inh_exc_gpe_spec[inh_count], color = cm.Blues_r(inh_count*30),lw = 3.)
ax1.legend(prop = {'size':7.}, title = 'BS frac. GPe')
ax1.set_ylabel('Firing rate (Bq)')

ax3.legend(prop = {'size':7.}, title = 'BS frac. GPe')
ax3.set_ylabel('Fano Factor')

ax5.legend(prop = {'size':7.}, title = 'BS frac. GPe')
ax5.set_ylabel('spectral entropy')
ax6.set_xlabel('BS ratio STN')
ax5.set_xlabel('BS ratio STN')



#inh_exc_stn_fr = []
#inh_exc_stn_ff = []
#inh_exc_gpe_fr = []
#inh_exc_gpe_ff = []
#inh_exc_gpe_spec = []
#inh_exc_stn_spec = []
#for inh_val in add_inh_rate:
  #exc_stn_fr = []
  #exc_stn_ff = []
  #exc_gpe_fr = []
  #exc_gpe_ff = []
  #exc_gpe_spec = []
  #exc_stn_spec = []
  #for exc_val in ext_rate:
    #ad1 = get_sim_3d({'pg_rate_exc':exc_val,'extra_inh':inh_val,'inh2_ratio':0.,'exc2_ratio':0.})
    #stn_fr_val = ad1.comp_mean_rate(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_exc')
    #stn_ff_val = ad1.comp_ff(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_exc', kernel_w= 2.)
    #exc_stn_fr.append(stn_fr_val)
    #exc_stn_ff.append(stn_ff_val)
    #gpe_fr_val = ad1.comp_mean_rate(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh')
    #gpe_ff_val = ad1.comp_ff(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh', kernel_w= 2.)
    #exc_gpe_fr.append(gpe_fr_val)
    #exc_gpe_ff.append(gpe_ff_val)
    #xx,xx,xx,peak_freq_gpe, peak_pow_gpe = ad1.psd(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh')
    #gpe_spec_val =  ad1.spec_entropy(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh', freq_range = [0.,100])
    #exc_gpe_spec.append(gpe_spec_val)
    #xx,xx,xx,peak_freq_stn, peak_pow_stn = ad1.psd(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_exc')
    #stn_spec_val = ad1.spec_entropy(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_exc', freq_range = [0.,100])
    #exc_stn_spec.append(stn_spec_val)
  #inh_exc_stn_fr.append(exc_stn_fr)
  #inh_exc_stn_ff.append(exc_stn_ff)
  #inh_exc_gpe_fr.append(exc_gpe_fr)
  #inh_exc_gpe_ff.append(exc_gpe_ff)
  #inh_exc_gpe_spec.append(exc_gpe_spec)
  #inh_exc_stn_spec.append(exc_stn_spec)

#fig = pl.figure(2,(10,15))
#ax1 = fig.add_subplot(321)
#ax2 = fig.add_subplot(322)
#ax3 = fig.add_subplot(323)
#ax4 = fig.add_subplot(324)
#ax5 = fig.add_subplot(325)
#ax6 = fig.add_subplot(326)
#for inh_count, inh_val in enumerate(add_inh_rate):
  #ax1.plot(ext_rate, inh_exc_stn_fr[inh_count], color = cm.Reds_r(inh_count*30),lw = 3.,label = str(inh_val))
  #ax3.plot(ext_rate, inh_exc_stn_ff[inh_count], color = cm.Greens_r(inh_count*30),lw = 3.,label = str(inh_val))
  #ax2.plot(ext_rate, inh_exc_gpe_fr[inh_count], color = cm.Reds_r(inh_count*30),lw = 3.)
  #ax4.plot(ext_rate, inh_exc_gpe_ff[inh_count], color = cm.Greens_r(inh_count*30),lw = 3.)
  #ax5.plot(ext_rate, inh_exc_stn_spec[inh_count], color = cm.Blues_r(inh_count*30),lw = 3.,label = str(inh_val))
  #ax6.plot(ext_rate, inh_exc_gpe_spec[inh_count], color = cm.Blues_r(inh_count*30),lw = 3.)
#ax1.legend(prop = {'size':7.}, title = 'Input to GPe')
#ax1.set_ylabel('Firing rate (Bq)')

#ax3.legend(prop = {'size':7.}, title = 'Input to GPe')
#ax3.set_ylabel('Fano Factor')

#ax5.legend(prop = {'size':7.}, title = 'Input to GPe')
#ax5.set_ylabel('spectral entropy')
#ax6.set_xlabel('Input to STN')
#ax5.set_xlabel('Input to STN')
pl.show()
