# script to map phase space (mps) for different Izhikevich neurons profiles
# creates copies of main file to be executed
# NOTE: run script from grid


import sys
import numpy as np
#import params_d
import params_hubs 
reload(params_hubs)
import itertools
import shutil
import os


sim_name = 'mps_psdb_ing-5'
#sim_name = 'test'

#path1 = '/home/vlachos/bcfgrid2/scripts/levcomp/batch/'
pars_file = 'params_d_psdb'+sim_name
#pars = params_d.Parameters()	
pars = params_hubs.Parameters()
data_path = pars.data_path

path1 = '/storage/users/ap1012/levcomp/scripts/batch/' # Path to store the results
path2 = '/storage/users/ap1012/levcomp/jdf_files/' # Path to jdf files
#path3 = '/bcfgrid/data/padmanabhan/scripts/levcomp/params_d.py' # Path of generic params_d file
path3 = '/storage/users/ap1012/levcomp/scripts/params_d_psdb.py'
path4 = '/storage/users/ap1012/levcomp/scripts/batch/'+pars_file+'.py'# Path of specific params_d file
path5 = '/storage/users/ap1012/output/' # Path for qsub output 


#batch_pars = {
    #'loss_ratio_exc':np.arange(0.1,0.21,0.1),
#}


batch_pars = {
    'g': np.arange(1.0,10.1,2.0),
    'pg_rate_exc': np.arange(345000.,350001,1000.),
    'inh2_ratio': np.arange(0.,1.01,0.1),
    'num_spikes_burst':np.arange(3.,5.1,2.),
    'extra_inh': np.arange(7000.,9001,2000.)
}

#}
# get all possible combinations and store them in a list
# caution: pars.values() *sorts* the key,value pairs; therefore the sorted keys
# have to be recovered

comb_pars = [[value for (key, value) in zip(batch_pars, values)]
  for values in itertools.product(*batch_pars.values())]
keys_sort = batch_pars.keys()

nr = 0
for ii in comb_pars:
    #fh = open(path1 + 'mps_%d_%d'%(g,eta)+'.py','w')
    fh = open(path1 + '%s_%d.py'%(sim_name,nr),'w')
    fh.write('import sys\n')
    fh.write("sys.path.insert(0,'/storage/users/ap1012/levcomp/scripts/')\n")
    fh.write('import simulation_nest_psdb as simd\n')
    fh.write('prefix = "%s_%d"\n'%(sim_name,nr))
    fh.write('pars_file="%s"\n'%pars_file)
    #fh.write('prefix = "sim_1_%d_%d"\n'%(g,eta))
    fh.write('pars = dict()\n')
    for jj in np.arange(np.shape(comb_pars)[1]):
	fh.write('pars["%s"] = %.2f\n'%(keys_sort[jj],ii[jj]))	
    fh.write('net1 = simd.Simulation(prefix,pars,pars_file)\n')
    fh.write('net1.sim()\n')
    fh.close()
    print sim_name,'%d/%d'%(nr,len(comb_pars))
    nr = nr + 1



    

# store pars dictionary for later use
fname = path1 + '%s.log'%sim_name
file(fname,'w').write(str(batch_pars))
# copy generic params_d file and store with modified name;
# this will be later used by the batch scripts
shutil.copy(path3,path4)


#create jdf file
for jj,ii in enumerate(comb_pars):
  fh = open(path2 + '%s_%s.jdf'%(sim_name,str(jj)),'w')
  #fh = open(path2 + '%s_%s.jdf'%(sim_name,str(nr)),'w')
  content = [
      '#!/bin/bash\n',
      #'#PBS -t 0-%d\n'%(len(comb_pars)-1),
      '#PBS -d %s\n'%path5,
      '#PBS -j oe\n',
      '#PBS -N %s\n'%sim_name,
      '#PBS -l mem=3500mb\n',
      '#PBS -l walltime=0:30:00\n',
      'export PYTHONPATH=/clustersw/cns/nest-2.2.2-padmanabhan/lib/python2.7/dist-packages/:$PYTHONPATH\n',
      'python /storage/users/ap1012/levcomp/scripts/batch/%s_%s.py\n'%(sim_name,str(jj))
      #'python /bcfgrid/data/padmanabhan/scripts/levcomp/batch/%s_%s.py'%(sim_name,str(nr))
      ]
  fh.writelines(content)
  fh.close()
  filename = path2 + '%s_%s.jdf'%(sim_name,str(jj))	
  os.system('qsub  %s'%filename )


