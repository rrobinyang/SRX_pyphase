#This file is part of the PyPhase software.
#
#Copyright (c) Max Langer (2019) 
#
#max.langer@creatis.insa-lyon.fr
#
#This software is a computer program whose purpose is to allow development,
#implementation, and deployment of phase retrieval algorihtms.
#
#This software is governed by the CeCILL  license under French law and
#abiding by the rules of distribution of free software.  You can  use, 
#modify and/ or redistribute the software under the terms of the CeCILL
#license as circulated by CEA, CNRS and INRIA at the following URL
#"http://www.cecill.info". 
#
#As a counterpart to the access to the source code and  rights to copy,
#modify and redistribute granted by the license, users are provided only
#with a limited warranty  and the software's author,  the holder of the
#economic rights,  and the successive licensors  have only  limited
#liability. 
#
#In this respect, the user's attention is drawn to the risks associated
#with loading,  using,  modifying and/or developing or reproducing the
#software by the user in light of its specific status of free software,
#that may mean  that it is complicated to manipulate,  and  that  also
#therefore means  that it is reserved for developers  and  experienced
#professionals having in-depth computer knowledge. Users are therefore
#encouraged to load and test the software's suitability as regards their
#requirements in conditions enabling the security of their systems and/or 
#data to be ensured and,  more generally, to use and operate it in the 
#same conditions as regards security. 
#
#The fact that you are presently reading this means that you have had
#knowledge of the CeCILL license and that you accept its terms.

import os
from math import *
import glob
#import utilities
import time
from pyphase.config import *
import inspect
from functools import wraps
import pathlib
import pickle
import numpy as np

def Serial(func):
    """Decorator for serial processing. Parallelisation decorators work on 
    functions with signature (self, *, dataset, projections)"""
    @wraps(func)
    def wrapper(*argv, **kwargs):
        func(*argv, **kwargs)
    return wrapper
    
def SLURM(func):
    """Decorator for parallellisation on SLURM"""
    @wraps(func)
    def wrapper(*argv, **kwargs):
        # Create a "worker" script that just calls align_projections
        print('SLURM')
        job_number = str(np.random.randint(0, 1e9))
        projections = kwargs.pop('projections') #Necessary argument
        parent_object = argv[0]
        current_dir = os.getcwd()
     
        tmpdir = '/.pyphase/SLURM/'
        serialisation_filename = 'pyphase_'+ job_number+'.pickle'

        pathlib.Path(current_dir+tmpdir).mkdir(parents=True, exist_ok=True)

        with open(current_dir+tmpdir+serialisation_filename, 'wb') as f:
            pickle.dump([parent_object, kwargs], f, protocol=pickle.HIGHEST_PROTOCOL)
        
        pathlib.Path(current_dir+tmpdir).mkdir(parents=True, exist_ok=True)    

        worker_file_name = current_dir + tmpdir + 'pyphase_worker_'+job_number+'.py'
        
        with open(worker_file_name, 'w') as worker_file:
            print("import sys", file=worker_file)
            print("import dataset", file=worker_file)
            print("import phaseretrieval", file=worker_file)
            print("import pickle", file=worker_file)
            print("with open('{}', 'rb') as f:".format(current_dir+tmpdir+serialisation_filename), file=worker_file)
            print("    stream = pickle.load(f)", file=worker_file)
            print("projections = [int(sys.argv[1]), int(sys.argv[2])]", file=worker_file)
            print("print(str(sys.argv))", file=worker_file)
            print("stream[0].{}.__wrapped__(stream[0], projections=projections, **stream[1])".format(func.__name__), file=worker_file)

        # Create a SLURM batch file
        number_of_CPU_per_node = min(number_of_nodes*number_of_cores, int(floor(memory_per_node/min_memory_per_core))) # Use the number of CPU corresponding to memory needs
        slurm_batch_file_name = current_dir + tmpdir + 'pyphase_slurm_'+job_number+'.sh'
        tasks = int(projections[1]) - int(projections[0]) + 1
        number_of_CPU = min(number_of_CPU_per_node * number_of_nodes, tasks)
        step = tasks // number_of_CPU - 1
        remainder = tasks % (number_of_CPU)
        with open(slurm_batch_file_name, 'w') as slurm_batch_file:
            print("#!/bin/sh", file=slurm_batch_file)
            print("#SBATCH -N {}".format(number_of_nodes), file=slurm_batch_file)
            #print("#SBATCH --exclusive", file=slurm_batch_file)
            print("#SBATCH --mem-per-cpu {}".format(min_memory_per_core), file=slurm_batch_file)
            print("#SBATCH --tasks-per-node {}".format(number_of_CPU_per_node), file=slurm_batch_file)            
            print("#SBATCH -t 20:00:00", file=slurm_batch_file)
            print("#SBATCH -J pyphase", file=slurm_batch_file)
            print("#SBATCH -o {}pyphase_%j.out".format(current_dir+tmpdir), file=slurm_batch_file)
            print("#SBATCH -e {}pyphase_%j.out".format(current_dir+tmpdir), file=slurm_batch_file)
            print("cat $0", file=slurm_batch_file)
            print("export tasks={}".format(tasks), file=slurm_batch_file)
            print("export start={}".format(projections[0]), file=slurm_batch_file)
            print("export step={}".format(step), file=slurm_batch_file)
            print("export remainder={}".format(remainder), file=slurm_batch_file)
            print("for (( i=0; i<{}; i++ ))".format(number_of_CPU), file=slurm_batch_file)
            print("do", file=slurm_batch_file)
            print("export end=$(( start+step+(i<remainder) ))", file=slurm_batch_file)
            print("srun -Q --exclusive -n 1 -N 1 python {} $start $end &> {}pyphase_worker_${{SLURM_JOB_ID}}_${{i}} &".format(worker_file_name, current_dir+tmpdir), file=slurm_batch_file)           
            print("echo \"python {} $start $end\"".format(worker_file_name), file=slurm_batch_file)
            print("export start=$(( end+1 ))", file=slurm_batch_file)
            print("sleep 1", file=slurm_batch_file)
            print("done", file=slurm_batch_file)        
            print("wait", file=slurm_batch_file)
            print("touch {}pyphase_{}.sig".format(current_dir+tmpdir, job_number), file= slurm_batch_file)

        # Submit
        cmd="sbatch ."+tmpdir+"pyphase_slurm_"+job_number+".sh"
        print("Executing: " + cmd)
        os.system(cmd)        
        timer=0        
        sigfile = 'pyphase_'+job_number+'.sig'
        while not os.path.exists(current_dir+tmpdir+sigfile):
            #if not(timer%30):
            #    print("Executing: " + func.__name__ + ", {} s".format(timer), end='\r')                
            time.sleep(1)
            timer += 1
    
        os.remove(current_dir+tmpdir+sigfile)
        os.remove(current_dir+tmpdir+serialisation_filename)
    
    return wrapper
    # pass
     
class OAR:
    """Legacy class for parallellisation on OAR. Will be reimplemented as 
    decorator"""
    def __init__(self):
        self.cores = 100
        self.executable_name = 'pyphase'
        self.path = os.path.split(os.path.realpath(__file__))[0]
        self.executable = os.path.join(self.path, self.executable_name)
        pass
    
    def WriteOarFiles(self, DS):
        pass
    
    #TODO: Don likne this implementation at all. Should probably be a decorator?
    def Launch(self, dataset, operator, **kwargs):
        # algorithm='', parameter='', distance=''
        #TODO: Don't like this solution at all, needs reimplementation... this is not elegant
        #TODO: standardised file names...oars
        #fname = dataset.path+'/'+dataset.name+'_/'+dataset.name+'_'+dataset.version+'_'+operator
        
        if operator=='retrieve':
            command_prefix = '{} {}'.format(operator, dataset.name)
            command_postfix = ''
            fname = dataset.phase_prefix
            if 'algorithm' in kwargs:
                command_postfix+=' --algorithm {}'.format(kwargs.get('algorithm'))
        elif operator=='retrieve_difference':
            fname = dataset.update_prefix
            if 'algorithm' in kwargs:
                command+=' --algorithm {}'.format(kwargs.get('algorithm'))
        elif operator=='propagate':
            fname = dataset.propagated_prefix # TODO:surely needs distance somehow?
            if 'algorithm' in kwargs:
                command+=' --algorithm {}'.format(kwargs.get('algorithm'))  
        elif operator=='difference':
            fname = dataset.difference_prefix
            
        if 'distance' in kwargs:
            command+=' --distance {}'.format(kwargs.get('distance'))
 
        #TODO: remove edf files (?)
        if os.path.exists(fname+"_.oar"):
            os.remove(fname+"_.oar")
        with open(fname+'_.oar', 'w') as oar_file:
            print("#!/bin/bash", file=oar_file)                 
            print("#OAR -l {mem_core_mb>=8000}/core=1,walltime=24", file=oar_file)
            print("#OAR -l {mem_core_mb>=4000 and mem_core_mb <8000}/nodes=1/core=2,walltime=48", file=oar_file)
            print("#OAR --array-param-file {}".format(fname+'_.params'), file=oar_file)
            print("#OAR --name {}".format(dataset.name), file=oar_file)
            print("#OAR --type besteffort", file=oar_file)
            print("#OAR --type idempotent", file=oar_file)
            print("{} $@".format(self.executable), file=oar_file)
        if os.path.exists(fname+"_.params"):
            os.remove(fname+"_.params")
        with open(fname+'_.params', 'w') as oar_file:
            interval = ceil(dataset.nbprojections/self.cores)           
            for x in range(self.cores):
                print(command_prefix+' -p {} {} '.format(x*interval, (x+1)*interval)+command_postfix, file=oar_file)
        #with open(fname+'.params', 'w') as oar_file:
        #    print("/mntdirect/_users/mlanger/Python/pyPhase/main.py 0 1", file=oar_file)
        
        # remove edf files before launching
        if operator in ('difference', 'propagate'):
            for file in glob.glob(fname+'_'+str(kwargs.get('distance'))+"_????.edf"):
                os.remove(file)
        else:
            for file in glob.glob(fname+"_????.edf"):
                os.remove(file)

        os.chmod(fname+'_.oar', 511)
        os.chmod(fname+'_.params', 511)
        # submit oar job that runs a script to signal end of job somehow 
        #program_path=os.path.dirname(__file__)
        cmd="oarsub -S "+fname+"_.oar > /dev/null"
        print("Executing: " + cmd)
        os.system(cmd)
        # wait for job to finish (track number of edf's finally, I guess) 
        # progress bar for niceness 
        
        #TODO: this won't work now. Needs proper file name handling
        
        edfCounter = 0

        if operator in ('difference', 'propagate'):
            nbprojections=dataset.nbprojections
            while edfCounter < nbprojections:
                utilities.update(operator.capitalize() + ' Distance ' + str(kwargs.get('distance')), edfCounter, nbprojections)
                edfCounter = len(glob.glob(fname+'_'+str(kwargs.get('distance'))+"_????.edf"))
                time.sleep(1)
            utilities.update(operator.capitalize() + ' Distance ' + str(kwargs.get('distance')), edfCounter, nbprojections)
        else:
            while edfCounter < dataset.nbprojections:
                utilities.update(operator.capitalize(), edfCounter, dataset.nbprojections)
                edfCounter = len(glob.glob(fname+"_????.edf"))
                time.sleep(1)
            utilities.update(operator.capitalize(), edfCounter, dataset.nbprojections)
#        pass
