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
	
	
prefix1 = 'conn_change-1'
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
#add_inh_rate = log_values['extra_inh']
#add_inh_rate = log_values['pg_rate_inh']



#big_big_stn_ff_lis = []
#big_big_stn_fr_lis = []
#big_big_gpe_ff_lis = []
#big_big_gpe_fr_lis = []
#big_big_gpe_spec_lis = []
#big_big_stn_spec_lis = []
#big_big_gpe_pow_lis = []
#big_big_stn_pow_lis = []


#for gpe_gpe_count, gpe_gpe_val in enumerate(conn_gpe_gpe):

  #big_stn_ff_lis = []
  #big_stn_fr_lis = []
  #big_gpe_ff_lis = []
  #big_gpe_fr_lis = []
  #big_gpe_spec_lis = []
  #big_stn_spec_lis = []
  #big_gpe_pow_lis = []
  #big_stn_pow_lis = []

  #for gpe_stn_count, gpe_stn_val in enumerate(conn_gpe_stn):
    #stn_ff_lis = []
    #stn_fr_lis = []
    #gpe_ff_lis = []
    #gpe_fr_lis = []
    #gpe_spec_lis = []
    #stn_spec_lis = []
    #gpe_pow_lis = []
    #stn_pow_lis = []
    #for stn_gpe_count, stn_gpe_val in enumerate(conn_stn_gpe):
      #stn_ff_val = np.zeros((len(exc_range), len(inh_range)))
      #stn_fr_val = np.zeros((len(exc_range), len(inh_range)))
      #gpe_fr_val = np.zeros((len(exc_range), len(inh_range)))
      #gpe_ff_val = np.zeros((len(exc_range), len(inh_range)))
      #gpe_spec_val = np.zeros((len(exc_range), len(inh_range)))
      #stn_spec_val = np.zeros((len(exc_range), len(inh_range)))
      #gpe_pow_val = np.zeros((len(exc_range), len(inh_range)))
      #stn_pow_val = np.zeros((len(exc_range), len(inh_range)))
      #for exc_count,exc_val in enumerate(exc_range):
	#for inh_count,inh_val in enumerate(inh_range):
	  ##ad1 = get_sim_3d({'pg_rate_exc':ext_val,'extra_inh':add_inh_val,'inh2_ratio':inh_val,'exc2_ratio':exc_val})
	  #ad1 = get_sim_3d({'epsilon_gpe_gpe':gpe_gpe_val,'epsilon_stn_gpe':stn_gpe_val,'epsilon_gpe_stn':gpe_stn_val,'inh2_ratio':inh_val,'exc2_ratio':exc_val})
	  #stn_ff_val[exc_count,inh_count] = ad1.comp_ff(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_exc', kernel_w= 2.)
	  #gpe_ff_val[exc_count,inh_count] = ad1.comp_ff(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh', kernel_w= 2.)
	  #stn_fr_val[exc_count,inh_count] = ad1.comp_mean_rate(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_exc')
	  #gpe_fr_val[exc_count,inh_count] = ad1.comp_mean_rate(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh')
	  #xx,xx,xx,peak_freq_gpe, peak_pow_gpe = ad1.psd(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh')
	  #gpe_pow_val[exc_count, inh_count] = peak_pow_gpe
	  #gpe_spec_val[exc_count, inh_count] = ad1.spec_entropy(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh', freq_range = [0.,1.8 * peak_freq_gpe])
	  #xx,xx,xx,peak_freq_stn, peak_pow_stn = ad1.psd(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh')
	  #stn_pow_val[exc_count, inh_count] = peak_pow_stn
	  #stn_spec_val[exc_count, inh_count] = ad1.spec_entropy(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh', freq_range = [0.,1.8 * peak_freq_stn])
	  
      #stn_ff_lis.append(stn_ff_val)
      #stn_fr_lis.append(stn_fr_val)
      #gpe_ff_lis.append(gpe_ff_val)
      #gpe_fr_lis.append(gpe_fr_val)
      #gpe_spec_lis.append(gpe_spec_val)
      #gpe_pow_lis.append(gpe_pow_val)
      #stn_spec_lis.append(stn_spec_val)
      #stn_pow_lis.append(stn_pow_val)
    #big_gpe_ff_lis.append(gpe_ff_lis)
    #big_gpe_fr_lis.append(gpe_fr_lis)
    #big_stn_ff_lis.append(stn_ff_lis)
    #big_stn_fr_lis.append(stn_fr_lis)
    #big_stn_pow_lis.append(stn_pow_lis)
    #big_stn_spec_lis.append(stn_spec_lis)
    #big_gpe_pow_lis.append(gpe_pow_lis)
    #big_gpe_spec_lis.append(gpe_spec_lis)
  #big_big_gpe_ff_lis.append(big_gpe_ff_lis)
  #big_big_gpe_fr_lis.append(big_gpe_fr_lis)
  #big_big_stn_ff_lis.append(big_stn_ff_lis)
  #big_big_stn_fr_lis.append(big_stn_fr_lis)
  #big_big_gpe_spec_lis.append(big_gpe_spec_lis)
  #big_big_gpe_pow_lis.append(big_gpe_pow_lis)
  #big_big_stn_spec_lis.append(big_stn_spec_lis)
  #big_big_stn_pow_lis.append(big_stn_pow_lis)

#np.save('analysis/big_big_stn_ff_val-2.npy',big_big_stn_ff_lis)
#np.save('analysis/big_big_gpe_ff_val-2.npy',big_big_gpe_ff_lis)
#np.save('analysis/big_big_stn_fr_val-2.npy',big_big_stn_fr_lis)
#np.save('analysis/big_big_gpe_fr_val-2.npy',big_big_gpe_fr_lis)
#np.save('analysis/big_big_stn_spec_val-2.npy',big_big_stn_spec_lis)
#np.save('analysis/big_big_gpe_spec_val-2.npy',big_big_gpe_spec_lis)
#np.save('analysis/big_big_stn_pow_val-2.npy',big_big_stn_pow_lis)
#np.save('analysis/big_big_gpe_pow_val-2.npy',big_big_gpe_pow_lis)



## Plotting the pcolors###
big_big_stn_ff_val = np.load('analysis/big_big_stn_ff_val-2.npy')
big_big_stn_fr_val = np.load('analysis/big_big_stn_fr_val-2.npy')
big_big_gpe_ff_val = np.load('analysis/big_big_gpe_ff_val-2.npy')
big_big_gpe_fr_val = np.load('analysis/big_big_gpe_fr_val-2.npy')
big_big_stn_spec_val = np.load('analysis/big_big_stn_spec_val-2.npy')
big_big_stn_pow_val = np.load('analysis/big_big_stn_pow_val-2.npy')
big_big_gpe_spec_val = np.load('analysis/big_big_gpe_spec_val-2.npy')
big_big_gpe_pow_val = np.load('analysis/big_big_gpe_pow_val-2.npy')


for gpe_gpe_count, gpe_gpe_val in enumerate(conn_gpe_gpe):
  nr = 1
  fig = pl.figure(10+gpe_gpe_count,(20,20)) 
  for gpe_stn_count, gpe_stn_val in enumerate(conn_gpe_stn):
    for stn_gpe_count, stn_gpe_val in enumerate(conn_stn_gpe):
      ax1 = fig.add_subplot(len(conn_gpe_stn), len(conn_stn_gpe), nr)
      big_stn_ff_val = big_big_stn_ff_val[gpe_gpe_count]
      #mask_big_stn_ff_val = np.ma.array(big_stn_ff_val,mask = np.isnan(big_stn_ff_val))
      pl.pcolor(exc_range,inh_range,(big_stn_ff_val[gpe_stn_count][stn_gpe_count]).T)
      pl.clim(0,10)
      if stn_gpe_count == len(conn_stn_gpe) -1: 

	pl.colorbar()
      if gpe_stn_count == len(conn_gpe_stn)-1:
	pl.gca().set_xticks(exc_range[::2])
	ax1.set_xlabel('Bursting ratio STN')
      else:
	pl.gca().set_xticks([])
      if stn_gpe_count == 0:
	pl.gca().set_yticks(inh_range[::2])
	ax1.set_ylabel('Burs. ratio GPe')
      else:
	pl.gca().set_yticks([])
      nr += 1
  fig.text(0.05,0.5,'GPe- STN connectivity', rotation = 'vertical')
  fig.text(0.5,0.05,'STN- GPe connectivity')
  fig.suptitle('STN FF vals. for gpe-gpe'+ str(gpe_gpe_val))
  fig.savefig('figures/stn_ff_'+str(gpe_gpe_count)+'.pdf',format = 'pdf')
for gpe_gpe_count, gpe_gpe_val in enumerate(conn_gpe_gpe):
  nr = 1
  fig = pl.figure(20+gpe_gpe_count,(20,20)) 
  for gpe_stn_count, gpe_stn_val in enumerate(conn_gpe_stn):
    for stn_gpe_count, stn_gpe_val in enumerate(conn_stn_gpe):
      ax1 = fig.add_subplot(len(conn_gpe_stn), len(conn_stn_gpe), nr)
      big_stn_ff_val = big_big_stn_fr_val[gpe_gpe_count]
      #mask_big_stn_ff_val = np.ma.array(big_stn_ff_val,mask = np.isnan(big_stn_ff_val))
      pl.pcolor(exc_range,inh_range,(big_stn_ff_val[gpe_stn_count][stn_gpe_count]).T)
      pl.clim(0,20)
      if stn_gpe_count == len(conn_stn_gpe) -1: 

	pl.colorbar()
      if gpe_stn_count == len(conn_gpe_stn)-1:
	pl.gca().set_xticks(exc_range[::2])
	ax1.set_xlabel('Bursting ratio STN')
      else:
	pl.gca().set_xticks([])
      if stn_gpe_count == 0:
	pl.gca().set_yticks(inh_range[::2])
	ax1.set_ylabel('Burs. ratio GPe')
      else:
	pl.gca().set_yticks([])
      nr += 1
  fig.text(0.05,0.5,'GPe- STN connectivity', rotation = 'vertical')
  fig.text(0.5,0.05,'STN- GPe connectivity')
  fig.suptitle('STN FR vals. for gpe-gpe'+ str(gpe_gpe_val))
  fig.savefig('figures/stn_fr_'+str(gpe_gpe_count)+'.pdf',format = 'pdf')
for gpe_gpe_count, gpe_gpe_val in enumerate(conn_gpe_gpe):
  nr = 1
  fig = pl.figure(30+gpe_gpe_count,(20,20)) 
  for gpe_stn_count, gpe_stn_val in enumerate(conn_gpe_stn):
    for stn_gpe_count, stn_gpe_val in enumerate(conn_stn_gpe):
      ax1 = fig.add_subplot(len(conn_gpe_stn), len(conn_stn_gpe), nr)
      big_stn_ff_val = big_big_gpe_ff_val[gpe_gpe_count]
      #mask_big_stn_ff_val = np.ma.array(big_stn_ff_val,mask = np.isnan(big_stn_ff_val))
      pl.pcolor(exc_range,inh_range,(big_stn_ff_val[gpe_stn_count][stn_gpe_count]).T)
      pl.clim(0,10)
      if stn_gpe_count == len(conn_stn_gpe) -1: 

	pl.colorbar()
      if gpe_stn_count == len(conn_gpe_stn)-1:
	pl.gca().set_xticks(exc_range[::2])
	ax1.set_xlabel('Bursting ratio STN')
      else:
	pl.gca().set_xticks([])
      if stn_gpe_count == 0:
	pl.gca().set_yticks(inh_range[::2])
	ax1.set_ylabel('Burs. ratio GPe')
      else:
	pl.gca().set_yticks([])
      nr += 1
  fig.text(0.05,0.5,'GPe- STN connectivity', rotation = 'vertical')
  fig.text(0.5,0.05,'STN- GPe connectivity')
  fig.suptitle('GPe FF vals. for gpe-gpe'+ str(gpe_gpe_val))
  fig.savefig('figures/gpe_ff_'+str(gpe_gpe_count)+'.pdf',format = 'pdf')

for gpe_gpe_count, gpe_gpe_val in enumerate(conn_gpe_gpe):
  nr = 1
  fig = pl.figure(40+gpe_gpe_count,(20,20)) 
  for gpe_stn_count, gpe_stn_val in enumerate(conn_gpe_stn):
    for stn_gpe_count, stn_gpe_val in enumerate(conn_stn_gpe):
      ax1 = fig.add_subplot(len(conn_gpe_stn), len(conn_stn_gpe), nr)
      big_stn_ff_val = big_big_gpe_fr_val[gpe_gpe_count]
      #mask_big_stn_ff_val = np.ma.array(big_stn_ff_val,mask = np.isnan(big_stn_ff_val))
      pl.pcolor(exc_range,inh_range,(big_stn_ff_val[gpe_stn_count][stn_gpe_count]).T)
      pl.clim(0,40)
      if stn_gpe_count == len(conn_stn_gpe) -1: 

	pl.colorbar()
      if gpe_stn_count == len(conn_gpe_stn)-1:
	pl.gca().set_xticks(exc_range[::2])
	ax1.set_xlabel('Bursting ratio STN')
      else:
	pl.gca().set_xticks([])
      if stn_gpe_count == 0:
	pl.gca().set_yticks(inh_range[::2])
	ax1.set_ylabel('Burs. ratio GPe')
      else:
	pl.gca().set_yticks([])
      nr += 1
  fig.text(0.05,0.5,'GPe- STN connectivity', rotation = 'vertical')
  fig.text(0.5,0.05,'STN- GPe connectivity')
  fig.suptitle('GPe FR vals. for gpe-gpe'+ str(gpe_gpe_val))
  fig.savefig('figures/gpe_fr_'+str(gpe_gpe_count)+'.pdf',format = 'pdf')
#nr = 1
#fig = pl.figure(14,(20,20)) 
#for ext_count, ext_val in enumerate(ext_rate):
  #for add_inh_count, add_inh_val in enumerate(add_inh_rate):
    #ax1 = fig.add_subplot(len(ext_rate), len(add_inh_rate), nr)
    #pl.pcolor(exc_range,inh_range,(big_stn_spec_val[ext_count][add_inh_count]).T)
    #pl.clim(0,1)
    #if add_inh_count == len(add_inh_rate) -1: 

      #pl.colorbar()
    #if ext_count == len(ext_rate)-1:
      #pl.gca().set_xticks(exc_range[::2])
      #ax1.set_xlabel('Bursting ratio STN')
    #else:
      #pl.gca().set_xticks([])
    #if add_inh_count == 0:
      #pl.gca().set_yticks(inh_range[::2])
      #ax1.set_ylabel('Burs. ratio GPe')
    #else:
      #pl.gca().set_yticks([])
    #nr += 1
#fig.text(0.05,0.5,'Excitatatory input', rotation = 'vertical')
#fig.text(0.5,0.05,'Inhibitory input')
#fig.suptitle('STN spec vals.')
#fig.savefig('stn_spec-2.pdf',format = 'pdf')


#nr = 1
#fig = pl.figure(15,(20,20)) 
#for ext_count, ext_val in enumerate(ext_rate):
  #for add_inh_count, add_inh_val in enumerate(add_inh_rate):
    #ax1 = fig.add_subplot(len(ext_rate), len(add_inh_rate), nr)
    #pl.pcolor(exc_range,inh_range,(big_stn_pow_val[ext_count][add_inh_count]).T)
    #pl.clim(0,5e10)
    #if add_inh_count == len(add_inh_rate) -1: 

      #pl.colorbar()
    #if ext_count == len(ext_rate)-1:
      #pl.gca().set_xticks(exc_range[::2])
      #ax1.set_xlabel('Bursting ratio STN')
    #else:
      #pl.gca().set_xticks([])
    #if add_inh_count == 0:
      #pl.gca().set_yticks(inh_range[::2])
      #ax1.set_ylabel('Burs. ratio GPe')
    #else:
      #pl.gca().set_yticks([])
    #nr += 1
#fig.text(0.05,0.5,'Excitatatory input', rotation = 'vertical')
#fig.text(0.5,0.05,'Inhibitory input')
#fig.suptitle('STN pow vals.')
#fig.savefig('stn_pow-2.pdf',format = 'pdf')


#nr = 1
#fig = pl.figure(16,(20,20)) 
#for ext_count, ext_val in enumerate(ext_rate):
  #for add_inh_count, add_inh_val in enumerate(add_inh_rate):
    #ax1 = fig.add_subplot(len(ext_rate), len(add_inh_rate), nr)
    #pl.pcolor(exc_range,inh_range,(big_gpe_spec_val[ext_count][add_inh_count]).T)
    #pl.clim(0,1)
    #if add_inh_count == len(add_inh_rate) -1: 

      #pl.colorbar()
    #if ext_count == len(ext_rate)-1:
      #pl.gca().set_xticks(exc_range[::2])
      #ax1.set_xlabel('Bursting ratio STN')
    #else:
      #pl.gca().set_xticks([])
    #if add_inh_count == 0:
      #pl.gca().set_yticks(inh_range[::2])
      #ax1.set_ylabel('Burs. ratio GPe')
    #else:
      #pl.gca().set_yticks([])
    #nr += 1
#fig.text(0.05,0.5,'Excitatatory input', rotation = 'vertical')
#fig.text(0.5,0.05,'Inhibitory input')
#fig.suptitle('GPe spec vals.')
#fig.savefig('gpe_spec-2.pdf',format = 'pdf')


#nr = 1
#fig = pl.figure(17,(20,20)) 
#for ext_count, ext_val in enumerate(ext_rate):
  #for add_inh_count, add_inh_val in enumerate(add_inh_rate):
    #ax1 = fig.add_subplot(len(ext_rate), len(add_inh_rate), nr)
    #pl.pcolor(exc_range,inh_range,(big_gpe_pow_val[ext_count][add_inh_count]).T)
    #pl.clim(0,5e10)
    #if add_inh_count == len(add_inh_rate) -1: 

      #pl.colorbar()
    #if ext_count == len(ext_rate)-1:
      #pl.gca().set_xticks(exc_range[::2])
      #ax1.set_xlabel('Bursting ratio STN')
    #else:
      #pl.gca().set_xticks([])
    #if add_inh_count == 0:
      #pl.gca().set_yticks(inh_range[::2])
      #ax1.set_ylabel('Burs. ratio GPe')
    #else:
      #pl.gca().set_yticks([])
    #nr += 1
#fig.text(0.05,0.5,'Excitatatory input', rotation = 'vertical')
#fig.text(0.5,0.05,'Inhibitory input')
#fig.suptitle('GPe Power vals.')
#fig.savefig('gpe_pow-2.pdf',format = 'pdf')
#pl.show()