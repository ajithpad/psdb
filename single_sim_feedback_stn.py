import sys
sys.path.append('/users/padmanabhan/simulations/basal/scripts/common/')
import misc2
import analyze_data_psd as adata
reload(adata)
import simulation_hysterisis_stn as simd
import numpy as np
import pdb

reload(simd)
import pylab as pl

prefix = 'test_osc'

pars = dict()
#pars['model_type'] = 'iaf_cond_alpha'
pars['model_type_exc'] = 'ssbn'
pars['model_type_inh'] = 'ssbn'

ff_lis = []
inh_lis = []
ff_inc_lis = []
ff_min = 2.0
ff_max = 9.56
max_spb = 4.
state_lis_stn = []
state_lis_gpe = []
pars['add_spikes'] = 0
pars['T_sim'] = 1000.
#pars['delay'] = 5.0
#pars['g'] = 9.0
pars['inh2_ratio'] = 1.
pars['exc2_ratio'] = 0.
pars['pg_rate_exc'] = 5200.#3.6
pars['pg_rate_inh'] = 2100.
pars['num_spikes_burst'] = 4.
pars['rnd_seeds'] = 1
ext_rate_steps = 100.
nr = 0

max_state_var_stn = 4.

state_var_stn = 0.
state_var_gpe = 0.5
ff = ff_max
chk_val = 0
#while ff <= ff_max:
while state_var_gpe <= max_spb+1:
#while state_var >= 0:
  if nr == 0:
    #pars['add_input'] = -200.
    net1 = simd.Simulation(prefix,pars)
    net1.sim()
    ad1 = adata.analyze_data(prefix)
    ff = ad1.comp_ff(pop_id = 'pops_exc', time_range = [200.,1000.])
    ff_lis.append(ff)
    print ff
  else:
    #pars['add_input'] = -200.
    if ff >= ff_max-3:
      pars['exc2_ratio'] = 0.8
      pars['state_var_stn'] = max_state_var_stn
      print 'critical point'
      chk_val += 1 
      state_var_stn = max_state_var_stn
    pars['state_var_gpe'] = state_var_gpe
    net1 = simd.Simulation(prefix,pars)
    net1.sim()
    ad1 = adata.analyze_data(prefix)
    pars['inh2_ratio'] = 1.
    ff = ad1.comp_ff(pop_id = 'pops_exc', time_range = [200.,1000.])
    ff_lis.append(ff)
    #pdb.set_trace()
    state_var_gpe = ((ff -ff_min)/ff_max)*max_spb*2.5
    print ff,state_var_stn, state_var_gpe
    state_lis_gpe.append(state_var_gpe)
    np.save('analysis/add_inp_state_lis_stn_gpe_stn_new.npy',state_lis_stn)
    np.save('analysis/add_inp_state_lis_stn_gpe_gpe_new.npy',state_lis_gpe)
    np.save('analysis/add_inp_lis_stn_gpe_new.npy',ff_lis) 

    
  nr+=1    
if chk_val == 0:
  pars['exc2_ratio'] = 0.8
  pars['state_var_stn'] = max_state_var_stn
  state_var_stn = max_state_var_stn
  pars['state_var_gpe'] = state_var_gpe
  net1 = simd.Simulation(prefix,pars)
  net1.sim()
  ad1 = adata.analyze_data(prefix)   
  ff = ad1.comp_ff(pop_id = 'pops_exc', time_range = [200.,1000.])
  ff_lis.append(ff)    
  np.save('analysis/add_inp_state_lis_stn_gpe_stn_new.npy',state_lis_stn)
  np.save('analysis/add_inp_state_lis_stn_gpe_gpe_new.npy',state_lis_gpe)
  np.save('analysis/add_inp_lis_stn_gpe_new.npy',ff_lis) 

#fig = pl.figure(1,(6,5))
#labelsize = 14.
#ff_lis = np.load('analysis/add_inp_ff_lis_stn_gpe_new.npy')
#ax1 = fig.add_subplot(111)
#ax1.plot(ff_lis,'-o', lw = 3., color = pl.cm.Greens(200))
#ax1.plot(5*np.ones(len(np.arange(2,12))),np.arange(2,12),'--', lw = 2.)
#ax1.set_ylabel('Fano Factor',size = labelsize)
#ax1.set_xlabel('Time', size = labelsize)
#ax1.set_yticks(np.arange(2.,11.,2.))
#ax1.set_yticklabels(np.arange(2.,11.,2.),size = labelsize)
#ax1.set_xticks([])

