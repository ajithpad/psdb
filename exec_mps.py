# Code to execute files in a loop
import sys
import subprocess
import numpy as np
import multiprocessing
path = '/space/padmanabhan/scripts/basal/batch/'
prefix0 = 'beta_oscil-3'

num_core = multiprocessing.cpu_count()

num_file = 4356

procs = []
for ii in range(2080,num_file):
  fil_name = path + prefix0 + '_' +str(ii) + '.py'
  print fil_name
  proc = subprocess.Popen([sys.executable, fil_name])
  procs.append(proc)
  if np.mod(ii+1,num_core) == 0: 
    for proc in procs:
      proc.wait()
    procs = []
    
  

  
  
