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
    #batch_pars ={'pg_rate_exc':ext_rate,'extra_inh':add_inh_rate,'exc2_ratio':exc_range,'inh2_ratio':inh_range }
    batch_pars ={'pg_rate_exc':ext_rate,'pg_rate_inh':add_inh_rate,'exc2_ratio':exc_range,'inh2_ratio':inh_range }
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
	
	
prefix1 = 'beta_oscil-3'
path1 = '/export/data-padmanabhan/results/basal/'
prefix0 = prefix1
log_path = '/export/data-padmanabhan/scripts/basal/batch/'
log_file = open(log_path + prefix1 +'.log')
log_values = eval(log_file.read())
exc_range = log_values['exc2_ratio']
inh_range = log_values['inh2_ratio']
ext_rate = log_values['pg_rate_exc']
#add_inh_rate = log_values['extra_inh']
add_inh_rate = log_values['pg_rate_inh']

pre_str = 'inp_change'
#pre_str = 'big'
post_str = 'diff_burs_len'
#post_str = '2'
#post_str = '3'

#big_stn_ff_lis = []
#big_stn_fr_lis = []
#big_gpe_ff_lis = []
#big_gpe_fr_lis = []
#big_gpe_spec_lis = []
#big_stn_spec_lis = []
#big_gpe_freq_lis = []
#big_stn_freq_lis = []

#for ext_count, ext_val in enumerate(ext_rate):
  #stn_ff_lis = []
  #stn_fr_lis = []
  #gpe_ff_lis = []
  #gpe_fr_lis = []
  #gpe_spec_lis = []
  #stn_spec_lis = []
  #gpe_freq_lis = []
  #stn_freq_lis = []
  #for add_inh_count, add_inh_val in enumerate(add_inh_rate):
    #stn_ff_val = np.zeros((len(exc_range), len(inh_range)))
    #stn_fr_val = np.zeros((len(exc_range), len(inh_range)))
    #gpe_fr_val = np.zeros((len(exc_range), len(inh_range)))
    #gpe_ff_val = np.zeros((len(exc_range), len(inh_range)))
    #gpe_spec_val = np.zeros((len(exc_range), len(inh_range)))
    #stn_spec_val = np.zeros((len(exc_range), len(inh_range)))
    #gpe_freq_val = np.zeros((len(exc_range), len(inh_range)))
    #stn_freq_val = np.zeros((len(exc_range), len(inh_range)))
    #for exc_count,exc_val in enumerate(exc_range):
      #for inh_count,inh_val in enumerate(inh_range):
	##ad1 = get_sim_3d({'pg_rate_exc':ext_val,'extra_inh':add_inh_val,'inh2_ratio':inh_val,'exc2_ratio':exc_val})
	#ad1 = get_sim_3d({'pg_rate_exc':ext_val,'pg_rate_inh':add_inh_val,'inh2_ratio':inh_val,'exc2_ratio':exc_val})
	#stn_ff_val[exc_count,inh_count] = ad1.comp_ff(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_exc', kernel_w= 2.)
	#gpe_ff_val[exc_count,inh_count] = ad1.comp_ff(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh', kernel_w= 2.)
	#stn_fr_val[exc_count,inh_count] = ad1.comp_mean_rate(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_exc')
	#gpe_fr_val[exc_count,inh_count] = ad1.comp_mean_rate(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh')
	#xx,xx,xx,peak_freq_gpe, peak_pow_gpe = ad1.psd(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh')
	#gpe_freq_val[exc_count, inh_count] = peak_freq_gpe
	#gpe_spec_val[exc_count, inh_count] = ad1.spec_entropy(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh', freq_range = [0.,1.8 * peak_freq_gpe])
	#xx,xx,xx,peak_freq_stn, peak_pow_stn = ad1.psd(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh')
	#stn_freq_val[exc_count, inh_count] = peak_freq_stn
	#stn_spec_val[exc_count, inh_count] = ad1.spec_entropy(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh', freq_range = [0.,1.8 * peak_freq_stn])
	
    #stn_ff_lis.append(stn_ff_val)
    #stn_fr_lis.append(stn_fr_val)
    #gpe_ff_lis.append(gpe_ff_val)
    #gpe_fr_lis.append(gpe_fr_val)
    #gpe_spec_lis.append(gpe_spec_val)
    #gpe_freq_lis.append(gpe_freq_val)
    #stn_spec_lis.append(stn_spec_val)
    #stn_freq_lis.append(stn_freq_val)
  #big_gpe_ff_lis.append(gpe_ff_lis)
  #big_gpe_fr_lis.append(gpe_fr_lis)
  #big_stn_ff_lis.append(stn_ff_lis)
  #big_stn_fr_lis.append(stn_fr_lis)
  #big_stn_freq_lis.append(stn_freq_lis)
  #big_stn_spec_lis.append(stn_spec_lis)
  #big_gpe_freq_lis.append(gpe_freq_lis)
  #big_gpe_spec_lis.append(gpe_spec_lis)

  
#np.save('analysis/'+pre_str +'_stn_ff_val-'+post_str+'.npy',big_stn_ff_lis)
#np.save('analysis/'+pre_str + '_gpe_ff_val-'+post_str+'.npy',big_gpe_ff_lis)
#np.save('analysis/'+pre_str + '_stn_fr_val-'+post_str+'.npy',big_stn_fr_lis)
#np.save('analysis/'+pre_str + '_gpe_fr_val-'+post_str+'.npy',big_gpe_fr_lis)
#np.save('analysis/'+pre_str + '_stn_spec_val-'+post_str+'.npy',big_stn_spec_lis)
#np.save('analysis/'+pre_str + '_gpe_spec_val-'+post_str+'.npy',big_gpe_spec_lis)
#np.save('analysis/'+pre_str + '_stn_freq_val-'+post_str+'.npy',big_stn_freq_lis)
#np.save('analysis/'+pre_str + '_gpe_freq_val-'+post_str+'.npy',big_gpe_freq_lis)



## Plotting the pcolors###
big_stn_ff_val = np.load('analysis/'+pre_str + '_stn_ff_val-'+post_str+'.npy')
big_stn_fr_val = np.load('analysis/'+pre_str + '_stn_fr_val-'+post_str+'.npy')
big_gpe_ff_val = np.load('analysis/'+pre_str + '_gpe_ff_val-'+post_str+'.npy')
big_gpe_fr_val = np.load('analysis/'+pre_str + '_gpe_fr_val-'+post_str+'.npy')
big_stn_spec_val = np.load('analysis/'+pre_str + '_stn_spec_val-'+post_str+'.npy')
big_stn_freq_val = np.load('analysis/'+pre_str + '_stn_freq_val-'+post_str+'.npy')
big_gpe_spec_val = np.load('analysis/'+pre_str + '_gpe_spec_val-'+post_str+'.npy')
big_gpe_freq_val = np.load('analysis/'+pre_str + '_gpe_freq_val-'+post_str+'.npy')

exc_range = exc_range.tolist()
exc_range.append(1.2)

exc_range = np.array(exc_range)


inh_range = inh_range.tolist()
inh_range.append(1.2)
inh_range = np.array(inh_range)

nr = 1
fig = pl.figure(10,(20,20)) 
for ext_count, ext_val in enumerate(ext_rate):
  for add_inh_count, add_inh_val in enumerate(add_inh_rate):
    ax1 = fig.add_subplot(len(ext_rate), len(add_inh_rate), nr)
    pl.pcolor(exc_range,inh_range,(big_stn_ff_val[ext_count][add_inh_count]).T)
    pl.clim(0,10 )
    if add_inh_count == len(add_inh_rate) -1: 

      pl.colorbar()
    if ext_count == len(ext_rate)-1:
      pl.gca().set_xticks(exc_range[::2])
      ax1.set_xlabel('Bursting ratio STN')
    else:
      pl.gca().set_xticks([])
    if add_inh_count == 0:
      pl.gca().set_yticks(inh_range[::2])
      ax1.set_ylabel('Burs. ratio GPe')
    else:
      pl.gca().set_yticks([])
    nr += 1
fig.text(0.05,0.5,'Excitatatory input', rotation = 'vertical')
fig.text(0.5,0.05,'Inhibitory input')
fig.suptitle('STN FF vals.')
fig.savefig('figures/'+pre_str+'_stn_ff-'+post_str+'.pdf',format = 'pdf')
nr = 1
fig = pl.figure(11,(20,20)) 
for ext_count, ext_val in enumerate(ext_rate):
  for add_inh_count, add_inh_val in enumerate(add_inh_rate):
    ax1 = fig.add_subplot(len(ext_rate), len(add_inh_rate), nr)
    pl.pcolor(exc_range,inh_range,(big_stn_fr_val[ext_count][add_inh_count]).T)
    pl.clim(0,40)
    if add_inh_count == len(add_inh_rate) -1: 

      pl.colorbar()
    if ext_count == len(ext_rate)-1:
      pl.gca().set_xticks(exc_range[::2])
      ax1.set_xlabel('Bursting ratio STN')
    else:
      pl.gca().set_xticks([])
    if add_inh_count == 0:
      pl.gca().set_yticks(inh_range[::2])
      ax1.set_ylabel('Burs. ratio GPe')
    else:
      pl.gca().set_yticks([])
    nr += 1
fig.text(0.05,0.5,'Excitatatory input', rotation = 'vertical')
fig.text(0.5,0.05,'Inhibitory input')
fig.suptitle('STN FR vals.')
fig.savefig('figures/'+pre_str+'_stn_fr-'+post_str+'.pdf',format = 'pdf')
nr = 1
fig = pl.figure(12,(20,20)) 
for ext_count, ext_val in enumerate(ext_rate):
  for add_inh_count, add_inh_val in enumerate(add_inh_rate):
    ax1 = fig.add_subplot(len(ext_rate), len(add_inh_rate), nr)
    pl.pcolor(exc_range,inh_range,(big_gpe_ff_val[ext_count][add_inh_count]).T)
    pl.clim(0,10)
    if add_inh_count == len(add_inh_rate) -1: 

      pl.colorbar()
    if ext_count == len(ext_rate)-1:
      pl.gca().set_xticks(exc_range[::2])
      ax1.set_xlabel('Bursting ratio STN')
    else:
      pl.gca().set_xticks([])
    if add_inh_count == 0:
      pl.gca().set_yticks(inh_range[::2])
      ax1.set_ylabel('Burs. ratio GPe')
    else:
      pl.gca().set_yticks([])
    nr += 1
fig.text(0.05,0.5,'Excitatatory input', rotation = 'vertical')
fig.text(0.5,0.05,'Inhibitory input')
fig.suptitle('GPe FF vals.')
fig.savefig('figures/'+pre_str+'_gpe_ff-'+post_str+'.pdf',format = 'pdf')
nr = 1
fig = pl.figure(13,(20,20)) 
for ext_count, ext_val in enumerate(ext_rate):
  for add_inh_count, add_inh_val in enumerate(add_inh_rate):
    ax1 = fig.add_subplot(len(ext_rate), len(add_inh_rate), nr)
    pl.pcolor(exc_range,inh_range,(big_gpe_fr_val[ext_count][add_inh_count]).T)
    pl.clim(10,80)
    if add_inh_count == len(add_inh_rate) -1: 

      pl.colorbar()
    if ext_count == len(ext_rate)-1:
      pl.gca().set_xticks(exc_range[::2])
      ax1.set_xlabel('Bursting ratio STN')
    else:
      pl.gca().set_xticks([])
    if add_inh_count == 0:
      pl.gca().set_yticks(inh_range[::2])
      ax1.set_ylabel('Burs. ratio GPe')
    else:
      pl.gca().set_yticks([])
    nr += 1
fig.text(0.05,0.5,'Excitatatory input', rotation = 'vertical')
fig.text(0.5,0.05,'Inhibitory input')
fig.suptitle('GPe FR vals.')
fig.savefig('figures/'+pre_str+'_gpe_fr-'+post_str+'.pdf',format = 'pdf')

nr = 1
fig = pl.figure(14,(20,20)) 
for ext_count, ext_val in enumerate(ext_rate):
  for add_inh_count, add_inh_val in enumerate(add_inh_rate):
    ax1 = fig.add_subplot(len(ext_rate), len(add_inh_rate), nr)
    pl.pcolor(exc_range,inh_range,(big_stn_spec_val[ext_count][add_inh_count]).T)
    pl.clim(0,1)
    if add_inh_count == len(add_inh_rate) -1: 

      pl.colorbar()
    if ext_count == len(ext_rate)-1:
      pl.gca().set_xticks(exc_range[::2])
      ax1.set_xlabel('Bursting ratio STN')
    else:
      pl.gca().set_xticks([])
    if add_inh_count == 0:
      pl.gca().set_yticks(inh_range[::2])
      ax1.set_ylabel('Burs. ratio GPe')
    else:
      pl.gca().set_yticks([])
    nr += 1
fig.text(0.05,0.5,'Excitatatory input', rotation = 'vertical')
fig.text(0.5,0.05,'Inhibitory input')
fig.suptitle('STN spec vals.')
fig.savefig('figures/'+pre_str+'_stn_spec-'+post_str+'.pdf',format = 'pdf')


nr = 1
fig = pl.figure(15,(20,20)) 
for ext_count, ext_val in enumerate(ext_rate):
  for add_inh_count, add_inh_val in enumerate(add_inh_rate):
    ax1 = fig.add_subplot(len(ext_rate), len(add_inh_rate), nr)
    pl.pcolor(exc_range,inh_range,(big_stn_freq_val[ext_count][add_inh_count]).T)
    pl.clim(0,100)
    if add_inh_count == len(add_inh_rate) -1: 

      pl.colorbar()
    if ext_count == len(ext_rate)-1:
      pl.gca().set_xticks(exc_range[::2])
      ax1.set_xlabel('Bursting ratio STN')
    else:
      pl.gca().set_xticks([])
    if add_inh_count == 0:
      pl.gca().set_yticks(inh_range[::2])
      ax1.set_ylabel('Burs. ratio GPe')
    else:
      pl.gca().set_yticks([])
    nr += 1
fig.text(0.05,0.5,'Excitatatory input', rotation = 'vertical')
fig.text(0.5,0.05,'Inhibitory input')
fig.suptitle('STN freq vals.')
fig.savefig('figures/'+pre_str+'_stn_freq-'+post_str+'.pdf',format = 'pdf')


nr = 1
fig = pl.figure(16,(20,20)) 
for ext_count, ext_val in enumerate(ext_rate):
  for add_inh_count, add_inh_val in enumerate(add_inh_rate):
    ax1 = fig.add_subplot(len(ext_rate), len(add_inh_rate), nr)
    pl.pcolor(exc_range,inh_range,(big_gpe_spec_val[ext_count][add_inh_count]).T)
    pl.clim(0,1)
    if add_inh_count == len(add_inh_rate) -1: 

      pl.colorbar()
    if ext_count == len(ext_rate)-1:
      pl.gca().set_xticks(exc_range[::2])
      ax1.set_xlabel('Bursting ratio STN')
    else:
      pl.gca().set_xticks([])
    if add_inh_count == 0:
      pl.gca().set_yticks(inh_range[::2])
      ax1.set_ylabel('Burs. ratio GPe')
    else:
      pl.gca().set_yticks([])
    nr += 1
fig.text(0.05,0.5,'Excitatatory input', rotation = 'vertical')
fig.text(0.5,0.05,'Inhibitory input')
fig.suptitle('GPe spec vals.')
fig.savefig('figures/'+pre_str+'_gpe_spec-'+post_str+'.pdf',format = 'pdf')


nr = 1
fig = pl.figure(17,(20,20)) 
for ext_count, ext_val in enumerate(ext_rate):
  for add_inh_count, add_inh_val in enumerate(add_inh_rate):
    ax1 = fig.add_subplot(len(ext_rate), len(add_inh_rate), nr)
    pl.pcolor(exc_range,inh_range,(big_gpe_freq_val[ext_count][add_inh_count]).T)
    pl.clim(0,100)
    if add_inh_count == len(add_inh_rate) -1: 

      pl.colorbar()
    if ext_count == len(ext_rate)-1:
      pl.gca().set_xticks(exc_range[::2])
      ax1.set_xlabel('Bursting ratio STN')
    else:
      pl.gca().set_xticks([])
    if add_inh_count == 0:
      pl.gca().set_yticks(inh_range[::2])
      ax1.set_ylabel('Burs. ratio GPe')
    else:
      pl.gca().set_yticks([])
    nr += 1
fig.text(0.05,0.5,'Excitatatory input', rotation = 'vertical')
fig.text(0.5,0.05,'GPe Freq vals.')
fig.savefig('figures/'+pre_str+'_gpe_freq-'+post_str+'.pdf',format = 'pdf')
pl.show()