import numpy as np
import itertools
import common.analyze_data_psd as adata
#import common.ff_analyze_data as adata
reload(adata)
import sys
from numpy import *
reload(adata)
import pdb

prefix0 = sys.argv[1]
#prefix0 = 'mps_izh_rs_fs_ch-4'
log_path = '/media/DISK_B/scripts_backup/levcomp/batch/'
log_file = open(log_path + prefix0 +'.log')
log_values = eval(log_file.read())
exc_range = log_values['pg_rate_exc']
g_range = log_values['g']
ratio_range = log_values['inh2_ratio']
spb_range = log_values['num_spikes_burst']
inh_range = log_values['extra_inh']
#phi_range = np.array([np.round(ii,1) for ii in log_values['ff_phi']])
path1 = '/media/DISK_C/results/levcomp/'

def get_sim_3d(match_pars):
    ''' get simulation object corresponding to (g,eta,inh2_ratio)
    match_pars is a dictionary, e.g. {'g':5.6,'eta':7.6, 'inh2_ratio':0.5}
    CAUTION: check order of keys'''  
    batch_pars ={'pg_rate_exc':exc_range,'num_spikes_burst':spb_range,'inh2_ratio':ratio_range,'g':g_range,'extra_inh':inh_range}
    comb_pars = [[value for (key, value) in zip(batch_pars, values)] for values in itertools.product(*batch_pars.values())] 
    sim_id = -1
    keys_sort = batch_pars.keys() 
    match_pars_lis = [match_pars[k] for k in keys_sort]
    for ii,pars in enumerate(comb_pars):
	if match_pars_lis==np.array(pars).round(2).tolist():
	    sim_id = ii	  	
	    print pars
    if sim_id>-1:
	res = adata.analyze_data(prefix0+'_%d'%sim_id,path1)
	print 'sim_id is: ',sim_id
	return res
    else:
	print 'No %s found'%match_pars.keys()
	return comb_pars,batch_pars.keys()
	
