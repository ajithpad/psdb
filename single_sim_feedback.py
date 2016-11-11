import sys
sys.path.append('/users/padmanabhan/simulations/basal/scripts/common/')
import misc2
import analyze_data_psd as adata
reload(adata)
import simulations_hysterisis as simd
import numpy as np
import pdb

reload(simd)
import pylab as pl

prefix = 'test_osc'

pars = dict()
#pars['model_type'] = 'iaf_cond_alpha'
pars['model_type'] = 'ssbn'


ff_lis = []
inh_lis = []
ff_inc_lis = []
ff_min = 2.66
ff_max = 10.03
max_spb = 4.
state_lis = []
pars['add_spikes'] = 0
pars['T_sim'] = 1000.
#pars['delay'] = 5.0
#pars['g'] = 9.0
pars['inh2_ratio'] = 1.
pars['exc2_ratio'] = 0.
pars['pg_rate_exc'] = 5100.#3.6
pars['pg_rate_inh'] = 2100.
pars['num_spikes_burst'] = 4.
pars['rnd_seeds'] = 1
ext_rate_steps = 100.
nr = 0



state_var = 0.5
while state_var <= max_spb+1:
#while state_var >= 0:
  if nr == 0:
    #pars['add_input'] = -200.
    net1 = simd.Simulation(prefix,pars)
    net1.sim()
    ad1 = adata.analyze_data(prefix)
    ff = ad1.comp_ff(pop_id = 'pops_inh', time_range = [200.,1000.])
    print ff
  else:
    #pars['add_input'] = -200.
    pars['state_var'] = state_var
    net1 = simd.Simulation(prefix,pars)
    net1.sim()
    ad1 = adata.analyze_data(prefix)
    pars['inh2_ratio'] = 1.
    ff = ad1.comp_ff(pop_id = 'pops_inh', time_range = [200.,1000.])
    ff_lis.append(ff)
    #pdb.set_trace()
    print ff
    state_var = ((ff -ff_min)/ff_max)*max_spb*4.5
    print state_var
    state_lis.append(state_var)

    
  nr+=1    

pars['state_var'] = state_var
net1 = simd.Simulation(prefix,pars)
net1.sim()
ad1 = adata.analyze_data(prefix)   
ff = ad1.comp_ff(pop_id = 'pops_inh', time_range = [200.,1000.])
ff_lis.append(ff)    

#fig = pl.figure(1,(6,5))
#labelsize = 14.
#ax1 = fig.add_subplot(111)
#ax1.plot(ff_lis[:-1], lw = 3., color = pl.cm.Greens(200))
#ax1.set_ylabel('Fano Factor',size = labelsize)
#ax1.set_xlabel('Time', size = labelsize)
#ax1.set_yticks(np.arange(2.,11.,2.))
#ax1.set_yticklabels(np.arange(2.,11.,2.),size = labelsize)
#ax1.set_xticks([])

