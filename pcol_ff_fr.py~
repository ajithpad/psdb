#Code to make pcolor maps of Fano Factor, firing rate of the STN and GPe populations

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
    batch_pars ={'exc2_ratio':exc_range,'inh2_ratio':inh_range }
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
	
	
prefix1 = 'start_asyn-1'
path1 = '/data/padmanabhan/results/basal/'
prefix0 = prefix1
log_path = '/data/padmanabhan/scripts/basal/batch/'
log_file = open(log_path + prefix1 +'.log')
log_values = eval(log_file.read())
exc_range = log_values['exc2_ratio']
inh_range = log_values['inh2_ratio']

big_stn_ff_lis = []
big_stn_fr_lis = []
big_gpe_ff_lis = []
big_gpe_fr_lis = []


for exc_count,exc_val in enumerate(exc_range):
  stn_ff_lis = []
  stn_fr_lis = []
  gpe_ff_lis = []
  gpe_fr_lis = []
  for inh_count,inh_val in enumerate(inh_range):
    ad1 = get_sim_3d({'inh2_ratio':inh_val,'exc2_ratio':exc_val})
    stn_ff_val = ad1.comp_ff(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_exc', kernel_w= 2.)
    gpe_ff_val = ad1.comp_ff(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh', kernel_w= 2.)
    stn_fr_val = ad1.comp_mean_rate(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_exc')
    gpe_fr_val = ad1.comp_mean_rate(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh')
    stn_ff_lis.append(stn_ff_val)
    stn_fr_lis.append(stn_fr_val)
    gpe_ff_lis.append(gpe_ff_val)
    gpe_fr_lis.append(gpe_fr_val)
  big_gpe_ff_lis.append(gpe_ff_lis)
  big_gpe_fr_lis.append(gpe_fr_lis)
  big_stn_ff_lis.append(stn_ff_lis)
  big_stn_fr_lis.append(stn_fr_lis)

fig = pl.figure(1)
ax1 = fig.add_subplot(221)
pl.pcolor(exc_range,inh_range, array(big_stn_ff_lis).T)
cb = pl.colorbar()
cb.set_label('FF STN')
#ax1.set_xlabel('Bursting ratio STN')
ax1.set_ylabel('Bursting ratio GPe')
ax2 = fig.add_subplot(222)
pl.pcolor(exc_range,inh_range, array(big_gpe_ff_lis).T)
cb = pl.colorbar()
cb.set_label('FF Gpe')
#ax2.set_xlabel('Bursting ratio STN')
#ax2.set_ylabel('Bursting ratio GPe')
ax3 = fig.add_subplot(223)
pl.pcolor(exc_range,inh_range, array(big_stn_fr_lis).T)
cb = pl.colorbar()
cb.set_label('FR STN')
ax3.set_xlabel('Bursting ratio STN')
ax3.set_ylabel('Bursting ratio GPe')
ax3 = fig.add_subplot(224)
pl.pcolor(exc_range,inh_range, array(big_gpe_fr_lis).T)
cb = pl.colorbar()
cb.set_label('FR Gpe')
ax3.set_xlabel('Bursting ratio STN')
#ax3.set_ylabel('Bursting ratio GPe')

pl.show()