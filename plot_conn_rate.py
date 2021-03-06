#Code to plot the connectivity against changes in the population activity

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
    #batch_pars ={'pg_rate_exc':ext_rate,'extra_inh':add_inh_rate,'exc2_ratio':exc_range,'inh2_ratio':inh_range }
    batch_pars ={'exc2_ratio':exc_range,'inh2_ratio':inh_range,'epsilon_gpe_gpe':conn_gpe_gpe,'epsilon_stn_gpe':conn_stn_gpe,'epsilon_gpe_stn':conn_gpe_stn}
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
	
	
prefix1 = 'conn_change-2'
path1 = '/data/padmanabhan/results/basal/'
prefix0 = prefix1
log_path = '/data/padmanabhan/scripts/basal/batch/'
log_file = open(log_path + prefix1 +'.log')
log_values = eval(log_file.read())
exc_range = log_values['exc2_ratio']
inh_range = log_values['inh2_ratio']
conn_gpe_gpe = log_values['epsilon_gpe_gpe']
conn_gpe_stn = log_values['epsilon_gpe_stn']
conn_stn_gpe = log_values['epsilon_stn_gpe']
  
gs_gg_sg_stn_fr = []
gs_gg_sg_stn_ff = []
gs_gg_sg_gpe_fr = []
gs_gg_sg_gpe_ff = []
gs_gg_sg_gpe_spec = []
gs_gg_sg_stn_spec = []


for gs_val in conn_gpe_stn:
  gg_sg_stn_fr = []
  gg_sg_stn_ff = []
  gg_sg_gpe_fr = []
  gg_sg_gpe_ff = []
  gg_sg_gpe_spec = []
  gg_sg_stn_spec = []
  for gg_val in conn_gpe_gpe:
    sg_stn_fr = []
    sg_stn_ff = []
    sg_gpe_fr = []
    sg_gpe_ff = []
    sg_gpe_spec = []
    sg_stn_spec = []
    for sg_val in conn_stn_gpe:
      ad1 = get_sim_3d({'epsilon_gpe_gpe':gg_val,'epsilon_stn_gpe':sg_val,'epsilon_gpe_stn':gs_val,'inh2_ratio':0.4,'exc2_ratio':0.4})
      stn_fr_val = ad1.comp_mean_rate(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_exc')
      stn_ff_val = ad1.comp_ff(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_exc', kernel_w= 2.)
      sg_stn_fr.append(stn_fr_val)
      sg_stn_ff.append(stn_ff_val)
      gpe_fr_val = ad1.comp_mean_rate(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh')
      gpe_ff_val = ad1.comp_ff(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh', kernel_w= 2.)
      sg_gpe_fr.append(gpe_fr_val)
      sg_gpe_ff.append(gpe_ff_val)
      xx,xx,xx,peak_freq_gpe, peak_pow_gpe = ad1.psd(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh')
      gpe_spec_val =  ad1.spec_entropy(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh', freq_range = [0.,100])
      sg_gpe_spec.append(gpe_spec_val)
      xx,xx,xx,peak_freq_stn, peak_pow_stn = ad1.psd(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_exc')
      stn_spec_val = ad1.spec_entropy(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_exc', freq_range = [0.,100])
      sg_stn_spec.append(stn_spec_val)
    gg_sg_stn_fr.append(exc_stn_fr)
    gg_sg_stn_ff.append(exc_stn_ff)
    gg_sg_gpe_fr.append(exc_gpe_fr)
    gg_sg_gpe_ff.append(exc_gpe_ff)
    gg_sg_gpe_spec.append(exc_gpe_spec)
    gg_sg_stn_spec.append(exc_stn_spec)
    
  gs_gg_sg_stn_fr.append(gg_sg_stn_fr)
  gs_gg_sg_stn_ff.append(gg_sg_stn_ff)
  gs_gg_sg_gpe_fr.append(gg_sg_gpe_fr)
  gs_gg_sg_gpe_ff.append(gg_sg_gpe_ff)
  gs_gg_sg_gpe_spec.append(gg_sg_gpe_spec)
  gs_gg_sg_stn_spec.append(gg_sg_stn_spec)
  