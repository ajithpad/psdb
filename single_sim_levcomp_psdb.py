import sys

sys.path.append('/users/padmanabhan/simulations/basal/scripts/common/')
import misc2
import analyze_data_psd as adata
reload(adata)
import simulation_nest_psdb as simd
reload(simd)
import pylab as pl

prefix = 'test'

pars = dict()
pars['model_type'] = 'psdb'
pars['add_spikes'] = 0
pars['T_sim'] = 1000.
#pars['delay'] = 5.0
#pars['g'] = 9.0
pars['inh2_ratio'] = 0.7
pars['exc2_ratio'] = 0.
pars['pg_rate_exc'] = 5000.#3.6
pars['pg_rate_inh'] = 2100.
pars['num_spikes_burst_stn'] = 4.
pars['num_spikes_burst_gpe'] = 4.
pars['add_spikes']= 0
pars['extra_spk'] = 'inh'
pars['st_val']= 700.
pars['num_gen']= 300
pars['ex_width']= 8.
pars['rnd_seeds'] = 1
pars['change_type'] = 1
pars['num_chg'] = 500


net1 = simd.Simulation(prefix,pars)
net1.sim()

ad2 = adata.analyze_data(prefix)
#pl.figure()
ad2.plot_raster(4000)

