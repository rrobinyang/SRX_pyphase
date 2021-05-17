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

import dataset as Dataset
import phaseretrieval as Phaseretrieval
import argparse
import os
import tomography as Tomography
import propagator as Propagator
import utilities as Utilities
from config import *
import matplotlib.pyplot as plt
import yaml
from pathlib import Path
import sys

#path = '/data/id19/bones01/bones1/max/Holocalibration/' #TODO: HOW TF DO I AUTOMATE THIS FOR TESTING

#TODO: Hey this shouldn't be at import? How to solve this, well, enough with a check if the file exists
#TODO: assume there is a parameter file for the moment
this = sys.modules[__name__]

this.parameter_filename = pyphase_path / 'parameters.yml'

def ReadPyphaseState():
    with open(parameter_filename, 'r') as f:
        this.pyphase_parameters=yaml.safe_load(f)
    this.version = pyphase_parameters['version']
    this.dataset_name = pyphase_parameters['dataset']

    version_filename = pyphase_path / dataset_name / version / 'parameters.yml'
    #read version file (?)
    with open(version_filename, 'r') as f:
        this.data_parameters=yaml.safe_load(f)

    this.data=Dataset.ESRF(this.dataset_name, './')


def Initialise(name, path='.', force=1):
    """Tries to find the dataset parameters in the accompanying *.info* and *xml* files.
    It is called when a parameter file ``pyphase_parameter.yml`` is not found.

    Returns
    -------
    self
        Returns *self* with the updated attributes.
    """
    
    #TODO: Implement (pseudocode)
    # If not exist, create ./.pyphase, ./.pyphase
    #TODO: OK, algo & prior should go in version parameter file. Dataset also represents reconstruction
    pyphase_path.mkdir(exist_ok=True)
    
    parameters = dict(version = 'master',
                      algorithm = default_method,
                      dataset = name,
                      prior = 'homogeneous')
    with open(parameter_filename, 'w') as f:
        yaml.dump(parameters, f)
    
    data = Dataset.ESRF(name)
    
    # Find data sets? Create .pyphase (./.pyphase/dataset_name? ./.pyphase/dataset_name/master ?)
    
    # Populate parameter files (.pyphase/parameters.yml, .pyphase/dataset_name/master/parameters.yml)
    pass

def Version():
    #List/Create/set version (like branch a bit?)
    pass
    
def Status():
    #Show current state (data set, version, )

    #Read pyphase parameter file
    ReadPyphaseState()
    print(yaml.dump(pyphase_parameters))
    print(yaml.dump(data_parameters))
    
    pass

#def Dataset(args): #Watch out, clash with module
#    pass

def CalculatetMotion():
    ReadPyphaseState()
    data.CalculateMotion()

def Align():
    ReadPyphaseState()
    data.Align()
    pass

def DisplayShifts():
    ReadPyphaseState()
    data.DisplayShifts()

def Preprocess():
    ReadPyphaseState()
    data.Preprocess()

def Reconstruct():
    ReadPyphaseState()
    data.Preprocess() #TODO: I guess these "compound" operations should be here, not in dataset?

def DisplayImage(projection=0):
    ReadPyphaseState()
    viewer.Display(data.GetImage(projection), 'Phase map ' + str(projection))
   
def DisplayProjection(projection=0, distance=1):   
    ReadPyphaseState()
    viewer.Display(data.GetProjection(projection, distance), 'Projection ' + str(projection) + ' Distance ' + str(distance))

def Retrieve(projection=0):
    """Phase retrieval with the current state"""
    #TODO: Decide on the interface...
    ReadPyphaseState()
    #Use directly the string in the parameter file as name of the class...:
    retriever = eval('Phaseretrieval.'+pyphase_parameters['algorithm']+'(data)')
    retriever.ReconstructProjection(data, projection)

def Process():
    pass

def Lcurve():
    pass

if __name__ == '__main__':
    # CLI
    
    path =  os.getcwd()

    def parse_status(args):
        Status()

    def parse_initialise(args):
        Initialise(args.name)

    def parse_motion(args):
        CalculatetMotion()
    
    def parse_align(args):
        if args.display:
            DisplayShifts()
        else:
            Align()
    
    # def parse_preprocess(args):
    #     #TODO: choice/automation of different data sources
    #     dataset = Dataset.ESRF(path, args.name)
    
    #     if args.initialise:
    #         dataset.Initialise()    
    #     elif args.motion:
    #         dataset.CalculateMotion()
    #     elif args.align:
    #         if args.projections:
    #             dataset.AlignProjections(args.projections)
    #         else:
    #             dataset.Align()
    #             dataset.DisplayShifts()
    #     elif args.lcurve:
    #         print('L-curve')
    #         retriever = Phaseretrieval.Mixed(dataset)
    #         # TODO: make possible to give projection for L-curve, possibility for interpolation of parameter etc
    #         retriever.Lcurve(dataset, 0)
    #     elif args.axis_position:
    #         dataset.CalculateAxisPositions()
    #     elif args.display_lcurve:
    #         retriever = Phaseretrieval.Mixed(dataset)
    #         retriever.DisplayLcurve(dataset)        
    #     else:
    #         dataset.Preprocess()
    
    
    def parse_retrieve(args):
        if type(args.proj) == int:
            Retrieve(args.proj)
        elif len(args.proj) == 2:
            Retrieve(args.proj[0], args.proj[1])
        else:
            print("Wrong number of arguments")

    def parse_display(args):
        #TODO: choice/automation of different data sources
        dataset = Dataset.ESRF(path, args.name)
    
        #TODO: How should the dataset malarky be handled
        if args.distance:
            if args.no_ff:
                preprocess = dataset.preprocess
                dataset.preprocess=0
                viewer.DisplayImage(dataset.GetProjection(args.projection[0]), args.distance, 'Projection '+str(args.projection[0])+', Distance '+str(args.distance)+' (uncorrected)')    
                dataset.preprocess = preprocess
            else:
                viewer.Display(dataset.GetProjection(args.projection[0]), args.distance , 'Projection '+str(args.projection[0])+', Distance '+str(args.distance))
        elif args.attenuation:
            viewer.Display(dataset.GetImage(args.projection[0], 'attenuation'), 'Attenuation map '+str(args.projection[0]))  
        else:
            viewer.Display(dataset.GetImage(args.projection[0]), 'Phase map ' + str(args.projection[0]))

            
    # def parse_reconstruct(args):
    #     #TODO: choice/automation of different data sources
    #     dataset = Dataset.ESRF(path, args.name)
    
    #     if args.beta:
    #         tomography.Reconstruct(dataset, 'retrieved_attenuation')
    #     elif args.delta_update:
    #         tomography.Reconstruct(dataset, 'update')
    #     elif args.beta_update:
    #         tomography.Reconstruct(dataset, 'attenuation_update')
    #     elif args.update_delta:
    #         dataset.UpdatePhase() #TODO: uh, phase? delta?
    #     elif args.update_beta:
    #         dataset.UpdateAttenuation()
    #     else: #args.delta
    #         tomography.Reconstruct(dataset, 'phase')


        # #TODO: choice/automation of different data sources
        # dataset = Dataset.ESRF(path, args.name)
    
        # if args.parameter:
        #     dataset.alpha = args.parameter
    
        # if args.algorithm[0].lower() == 'mixed':
        #     print("Phase retrieval: Mixed approach")
        #     retriever = Phaseretrieval.Mixed(dataset)
        #     if not args.projection:
        #         print("Phase retrieval: {}".format(dataset.name))
        #         retriever.Reconstruct(dataset)
        #     elif len(args.projection) == 1:
        #         print("Phase retrieval: projection {}".format(args.projection[0]))
        #         retriever.ReconstructProjection(dataset,args.projection[0])
        #     elif len(args.projection) == 2:
        #         print("Phase retrieval: projections {} {}".format(args.projection[0], args.projection[1]))
        #         retriever.ReconstructProjections(dataset, args.projection[0], args.projection[1])
        #     else:
        #         print('Error, incorrect number of projections') #TODO: Correct error handling
    
        # elif args.algorithm[0].lower() == 'iterative':
        #     print("Phase retrieval: 3D CTF-SIRT")
        #     retriever = Phaseretrieval.Iterative()
        #     retriever.Reconstruct(dataset, args.iterations[0], args.parameter[0])
        
        # elif args.algorithm[0].lower() == 'ctf' and args.difference:
        #     print('Phase retrieval from contrast difference: CTF')
        #     retriever = Phaseretrieval.CTF()
        #     if not args.projection:
        #         retriever.ReconstructDifference(dataset)
        #     elif len(args.projection) == 1:
        #         retriever.ReconstructDifferenceProjection(dataset,args.projection[0])
        #     elif len(args.projection) == 2:
        #         retriever.ReconstructDifferenceProjections(dataset, args.projection[0], args.projection[1])
        #     else:
        #         print('Error, incorrect number of projections') #TODO: Correct error handling
        # else:
        #     print('Phase retrieval: CTF (default)')
        #     retriever = Phaseretrieval.CTF(dataset)
        #     if not args.projection:
        #         retriever.Reconstruct(dataset)
        #     elif len(args.projection) == 1:
        #         retriever.ReconstructProjection(dataset,args.projection[0])
        #     elif len(args.projection) == 2:
        #         retriever.ReconstructProjections(dataset, args.projection[0], args.projection[1])
        #     else:
        #         print('Error, incorrect number of projections') #TODO: Correct error handling
    
    # def parse_project(args):
    #     #TODO: choice/automation of different data sources
    #     dataset = Dataset.ESRF(path, args.name)
    
    #     if delta:
    #         tomography.ForwardProject(dataset, 'phase')
    #     elif beta:
    #         tomography.ForwardProject(dataset, 'retrieved_attenuation')
    
    # def parse_propagate(args):
    #     #TODO: choice/automation of different data sources
    #     dataset = Dataset.ESRF(path, args.name)
    
    #     if args.projections and args.distance:
    #         propagator.PropagateProjections(dataset, args.projections[0], args.distance[0])
    #     else:
    #         propagator.Propagate(dataset)
    
    # def parse_process(args):
    #     #TODO: choice/automation of different data sources
    #     dataset = Dataset.ESRF(path, args.name)
    
    #     if args.difference:
    #         if args.projections and args.distance:
    #             dataset.DifferenceProjections(args.projections, args.distance[0])
    #         else:
    #             dataset.Difference()
    
        
    parser = argparse.ArgumentParser(description='PyPhase phase retrieval', prog='pyphase')
    subparsers = parser.add_subparsers(title='actions')
    
    parser_status = subparsers.add_parser('status', description='Current state', help='Display current state')
    parser_status.set_defaults(func=parse_status)

    parser_initialise = subparsers.add_parser('initialise', description='Initialise pyphase', help='Initialises pyphase in this directory')
    parser_initialise.add_argument('name', metavar='name', help='name of the data set')
    parser_initialise.set_defaults(func=parse_initialise)

    parser_motion = subparsers.add_parser('motion', description='Correct motion', help='Calculates drift correction of tomographic scans')
    parser_motion.set_defaults(func=parse_motion)

    parser_align = subparsers.add_parser('align', description='Align projections', help='Aligns multidistance projections')
    parser_align_mutex = parser_align.add_mutually_exclusive_group()
    parser_align_mutex.add_argument('--display', action='store_true', help='Display alignment curves and an aligned projection')
    parser_align.set_defaults(func=parse_align)
   
    #parser_preprocess = subparsers.add_parser('preprocess')
    #parser_preprocess.add_argument('projections', metavar='proj', type=int, nargs='*', action='store', help='projection to propagate')
    #parser_preprocess_mutex = parser_preprocess.add_mutually_exclusive_group()
    #parser_preprocess_mutex.add_argument('--initialise', action='store_true', help='force initialisation')
    #parser_preprocess_mutex.add_argument('--align', action='store_true', help='align projections')
    #parser_preprocess_mutex.add_argument('--motion', action='store_true', help='calculate motion')
    #parser_preprocess_mutex.add_argument('--lcurve', action='store_true', help='find regularisation parameter with L-curve')
    #parser_preprocess_mutex.add_argument('--display_lcurve', action='store_true', help='display L-curve')
    #parser_preprocess_mutex.add_argument('--axis_position', action='store_true', help='Caclulate rotation axes')
    #parser_preprocess.set_defaults(func=parse_preprocess)
      
    #parser_retrieve = subparsers.add_parser('retrieve', description='phase retrieval', help='Phase retrieval')
    #parser_retrieve.add_argument('proj', metavar='proj', nargs='*', default=0, action='store', help='Can be one value for one projection or two values for a range') #TODO: Doesn't look so good on help string
    #parser_retrieve_mutex = parser_retrieve.add_mutually_exclusive_group()
    #parser_retrieve_mutex.add_argument('--difference', action='store_true', help='retrieve phase from difference images')
    #parser_retrieve.set_defaults(func=parse_retrieve)
    
    # parser_reconstruct = subparsers.add_parser('reconstruct', help='tomo help')
    # parser_reconstruct.add_argument('name', metavar='name', help='name of the data set')
    # parser_reconstruct_mutex = parser_reconstruct.add_mutually_exclusive_group()
    # parser_reconstruct_mutex.add_argument('--delta', action='store_true', help='tomographic reconstruction of delta from phase maps')
    # parser_reconstruct_mutex.add_argument('--beta', action='store_true', help='tomographic reconstruction of beta from attenuation maps')
    # parser_reconstruct_mutex.add_argument('--delta_update', action='store_true', help='reconstruct delta update from retrieved difference image')
    # parser_reconstruct_mutex.add_argument('--beta_update', action='store_true', help='reconstruct beta update from retrieved difference image')
    # parser_reconstruct_mutex.add_argument('--update_delta', action='store_true', help='update current solution for delta')
    # parser_reconstruct_mutex.add_argument('--update_beta', action='store_true', help='update current solution for beta')
    # parser_reconstruct.set_defaults(func=parse_reconstruct)
    
    # parser_propagate = subparsers.add_parser('propagate')
    # parser_propagate.add_argument('name', metavar='name', help='name of the data set')
    # parser_propagate.add_argument('projection', metavar='proj', type=int, nargs='*', action='store', help='projection to propagate')
    # parser_propagate.add_argument('distance', metavar='dist', type=int, nargs='*', action='store', help='distance to propagate to')
    # parser_propagate_mutex = parser_propagate.add_mutually_exclusive_group()
    # parser_propagate.set_defaults(func=parse_propagate)
    
    # parser_project = subparsers.add_parser('project')
    # parser_project.add_argument('name', metavar='name', help='name of the data set')
    # parser_project_mutex = parser_project.add_mutually_exclusive_group()
    # parser_project_mutex.add_argument('--delta', action='store_true', help='create projections of delta')
    # parser_project_mutex.add_argument('--beta', action='store_true', help='create projections of beta')
    # parser_project.set_defaults(func=parse_project)
    
    # parser_display = subparsers.add_parser('display', help='display an image')
    # parser_display.add_argument('name', metavar='name', help='name of the data set')
    # parser_display.add_argument('projection', metavar='proj', type=int, nargs=1, action='store', help='projection to display')
    # parser_display.add_argument('distance', metavar='dist', type=int, nargs='?', action='store', help='if given, displays a recorded intensity')
    # parser_display_mutex = parser_display.add_mutually_exclusive_group()
    # parser_display_mutex.add_argument('--no_ff', action='store_true', help='shows uncorrected projections')
    # parser_display_mutex.add_argument('-a', '--attenuation', action='store_true', help='shows retrieved attenuation')
    # parser_display.set_defaults(func=parse_display)
    
    # parser_process = subparsers.add_parser('process', help='Other processing')
    # parser_process.add_argument('name', metavar='name', help='Name of the data set')
    # parser_process.add_argument('projection', metavar='proj', type=int, nargs='?', action='store', help='projection to display')
    # parser_process.add_argument('distance', metavar='dist', type=int, nargs='?', action='store', help='if given, displays a recorded intensity')
    # parser_process_mutex = parser_preprocess.add_mutually_exclusive_group()
    # parser_process_mutex.add_argument('--difference', action='store_true', help='Calculates difference of measured and generated intensity')
    
    args = parser.parse_args()
    
    args.func(args)
    
     # How should this be handled anyway? It should try to figure out the data source automatically, right? How can this be done anyway?
    #dataset.preprocess=0
    #dataset.correct_shifts=0
    #dataset.correct_motion=0
    
    # TODO: CONFIGURATION
    #tomography = Tomography.PyHST()
    #propagator = Propagator.CTF(dataset)
    #retriever = Phaseretrieval.CTF(dataset)


