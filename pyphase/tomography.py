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
import sys 
import time

#TODO: Parallelisation and reconstruction needs to be split

#TODO: Import of tomopy can only be done if tomopy is to be used!
#import tomopy


class Tomography:
    """Class for tomographic operations, for use with 3D iterative methods
    
    
    Note
    ----
    Legacy code to be aligned with current API
    """
    pass

# class TomoPy:
#     def __init__(self):
#         pass
    
#     def Reconstruct(self, data, volume='phase'):
#         # Since tomopy takes volumes I baically have to create them? 
        
#         # How to parallellise ? 
#         # How to split with paralleliser?
        
#         sinogram = data.GetSinogram() #TODO: is this wise?
        
#         pass

class PyHST:
    """Tomographic operations with PyHST"""
    def __init__(self):
        pass

    def ForwardProject(self, DS, volume='phase'):
        """Generate projections from a reconstructed dataset"""
        # Attenuation
        # retrieved attenuation
        # Phase
        
        #self._create_forward_frojection_parameter_file(DS, volume)
        self._create_parameter_file(DS, volume, direction='forward')
        self._create_reconstruction_files(DS, volume, direction='forward')
        
        if volume == 'phase':
            prefix = DS.phase_forward_prefix
            vol_file_name = DS.phase_vol_filename
        elif volume=='attenuation':
            prefix=DS.attenuation_forward_prefix
            vol_file_name=DS.attenuation_vol_filename
        elif volume=='retrieved_attenuation':
            prefix=DS.retrieved_attenuation_forward_prefix
            vol_file_name=DS.retrieved_attenuation_vol_filename
        elif volume=='prior':
            prefix=DS.prior_forward_prefix
            vol_file_name=DS.prior_vol_filename
            
            
        if os.path.exists(prefix+"_.sig"):
            os.remove(prefix+"_.sig")
        
        cmd="oarsub -S "+prefix+"_.oar"
        os.system(cmd)
        edfCounter = 0
        timer=0        
        while not os.path.exists(prefix+"_.sig"):
            if not(timer%1):
                print("Executing: " + cmd + ", {} s".format(timer), end='\r')                
            time.sleep(1)
            timer += 1

        print("Executing: " + cmd + ", {} s".format(timer))                

    pass

    def _create_reconstruction_files(self, DS, volume='phase', direction='back'):
        """Generate supporting files necessary for reconstruction with PyHST"""
        if direction == 'back': 
            if volume == 'phase':
                prefix = DS.phase_prefix
                vol_file_name = DS.phase_vol_filename
                name=DS.phase_name
            elif volume=='attenuation':
                prefix=DS.attenuation_prefix
                vol_file_name=DS.attenuation_vol_filename
                name=DS.attenuation_name
            elif volume=='retrieved_attenuation':
                prefix=DS.retrieved_attenuation_prefix
                vol_file_name=DS.retrieved_attenuation_vol_filename
                name=DS.retrieved_attenuation_name
            elif volume=='update':
                prefix=DS.update_prefix
                vol_file_name=DS.update_vol_filename
                name=DS.update_name
            elif volume=='attenuation_update':
                prefix=DS.attenuation_update_prefix
                vol_file_name=DS.attenuation_update_vol_filename
                name=DS.attenuation_update_name
            pass
        
        elif direction == 'forward':
            if volume == 'phase':
                prefix = DS.phase_forward_prefix
                vol_file_name = DS.phase_vol_filename
                name=DS.phase_forward_name
            elif volume=='attenuation':
                prefix=DS.attenuation_forward_prefix
                vol_file_name=DS.attenuation_vol_filename
                name=DS.attenuation_forward_name
            elif volume=='retrieved_attenuation':
                prefix=DS.retrieved_attenuation_forward_prefix
                vol_file_name=DS.retrieved_attenuation_vol_filename
                name=DS.retrieved_attenuation_forward_name
            elif volume=='prior':
                prefix=DS.prior_forward_prefix
                vol_file_name=DS.prior_vol_filename
                name=DS.prior_forward_name                
            pass

        if os.path.exists(prefix+"_.oar"):
            os.remove(prefix+"_.oar")

        with open(prefix+'_.oar', 'w') as oar_file:
            print("#!/bin/bash", file=oar_file)
            print("# This file has been generated by PyPhase", file=oar_file)
            print("#", file=oar_file)
            print("# OAR directives start with one hash followed by OAR (#OAR ....)", file=oar_file)
            print("# commented (optional) OAR directives start with two hashes followed by OAR (##OAR ....)", file=oar_file)
            print("# comments start with one hash and blanc (# here is my comment) ", file=oar_file)
            print("#", file=oar_file)
            print("#OAR -l {gpu='YES'}/nodes=1,walltime=24", file=oar_file)
            print("#", file=oar_file)
            print("# optional start calculation on CPUs", file=oar_file)
            print("##OAR -l {gpu='NO'}/nodes=1,walltime=48", file=oar_file)
            print("#", file=oar_file)
            print("#OAR --name {}".format(name+'_'), file=oar_file)
            print("#OAR --project pyhst2", file=oar_file)
            print("echo Executing /usr/bin/pyhst2 {}".format(prefix+'_.par'), file=oar_file)
            print("/usr/bin/pyhst2 {}".format(prefix+'_.par'), file=oar_file)
            print("touch {}".format(prefix+'_.sig'), file=oar_file)
            os.chmod(prefix+'_.oar', 511)
        pass

    # def _create_forward_frojection_parameter_file(self, DS, volume='phase'):

    #     if volume == 'phase':
    #         parname=DS.phase_forward_prefix
    #         file_prefix=DS.phase_forward_prefix
    #         output_file=DS.phase_vol_filename
    #         proj_output_file=DS.phase_forward_prefix
    #         axis_position=DS.axis_position
    #     elif volume == 'attenuation':
    #         parname = DS.attenuation_forward_prefix
    #         file_prefix=DS.attenuation_forward_prefix
    #         output_file=DS.attenuation_vol_filename
    #         proj_output_file=DS.attenuation_forward_prefix
    #         axis_position=DS.attenuation_axis_position
    #     elif volume == 'retrieved_attenuation': 
    #         parname=DS.retrieved_attenuation_forward_prefix
    #         file_prefix=DS.retrieved_attenuation_forward_prefix
    #         output_file=DS.retrieved_attenuation_vol_filename
    #         proj_output_file=DS.retrieved_attenuation_forward_prefix
    #         axis_position=DS.attenuation_axis_position
    #     elif volume=='prior':
    #         prefix=DS.prior_forward_prefix
    #         vol_file_name=DS.prior_vol_filename
    #         subtract_background='NO'
    #         correct_flatfield='NO'
    #         take_logarithm = 'NO' # Take log of projection values
    #         axis_position=DS.axis_position[DS.reference_distance-1] #TODO: Verify (should be like this, all proj to reference plane axis?)
    #         proj_output_file=DS.prior_forward_prefix

    #     else:
    #         print("ForwardProject: no such volume {}".format(volume))
    #         return

    #     proj_output_file+='_%04d.edf'

    #     with open(parname+'.par', 'w') as oar_file:
    #         print("## PYHST PARAMETER FILE", file=oar_file)                 
    #         print("## Parameters defining the projection file series", file=oar_file)                 
    #         print("MULTIFRAME=0", file=oar_file)                 
    #         print("FILE_PREFIX = {}".format(file_prefix), file=oar_file)
    #         print("NUM_FIRST_IMAGE = 0 # No. of first projection file", file=oar_file)
    #         print("NUM_LAST_IMAGE = {} # No. of last projection file".format(DS.nbprojections-1), file=oar_file)
    #         print("TRYEDFCONSTANTHEADER = 0 # assume constant edf header size", file=oar_file)
    #         print("NUMBER_LENGTH_VARIES = NO", file=oar_file)
    #         print("LENGTH_OF_NUMERICAL_PART = 4 # No. of characters", file=oar_file)
    #         print("FILE_POSTFIX = .edf", file=oar_file)
    #         print("FILE_INTERVAL = 1 # Interval between input files", file=oar_file)
    #         print("", file=oar_file)            
    #         print("## Parameters defining the projection file format", file=oar_file)
    #         print("NUM_IMAGE_1 = {} # Number of pixels horizontally".format(DS.nx), file=oar_file)
    #         print("NUM_IMAGE_2 = {} # Number of pixels vertically".format(DS.ny), file=oar_file)
    #         print("IMAGE_PIXEL_SIZE_1 = {} # Pixel size horizontally (microns)".format(DS.pixel_size), file=oar_file)
    #         print("IMAGE_PIXEL_SIZE_2 = {} # Pixel size vertically".format(DS.pixel_size), file=oar_file)
    #         print("", file=oar_file)
    #         print("## Parameters defining background treatment", file=oar_file)
    #         print("SUBTRACT_BACKGROUND = NO # No background subtraction", file=oar_file)
    #         print("BACKGROUND_FILE = N.A.", file=oar_file)
    #         print("", file=oar_file)
    #         print("## Parameters defining flat-field treatment", file=oar_file)
    #         print("CORRECT_FLATFIELD = NO # No flat-field correction", file=oar_file)
    #         print("FLATFIELD_CHANGING = N.A.", file=oar_file)
    #         print("FLATFIELD_FILE = N.A.", file=oar_file)
    #         print("FF_PREFIX = N.A.", file=oar_file)
    #         print("FF_NUM_FIRST_IMAGE = N.A.", file=oar_file)
    #         print("FF_NUM_LAST_IMAGE = N.A.", file=oar_file)
    #         print("FF_NUMBER_LENGTH_VARIES = N.A.", file=oar_file)
    #         print("FF_LENGTH_OF_NUMERICAL_PART = N.A.", file=oar_file)
    #         print("FF_POSTFIX = N.A.", file=oar_file)
    #         print("TAKE_LOGARITHM = YES # Take log of projection values", file=oar_file)
    #         print("FF_FILE_INTERVAL = N.A.", file=oar_file)
    #         print("", file=oar_file)
    #         print("## Parameters defining experiment", file=oar_file)
    #         print("ANGLE_BETWEEN_PROJECTIONS = {} # Increment angle in degrees".format(DS.angle_increment), file=oar_file)
    #         print("ROTATION_VERTICAL = YES", file=oar_file)
    #         #print("ROTATION_AXIS_POSITION = {:f} # Position in pixels".format((DS.nx-1)/2), file=oar_file)
    #         #print("ROTATION_AXIS_POSITION = {:f} # Position in pixels".format(1060), file=oar_file)
    #         #TODO: Axis position is not necessarily the same for the attenuation distance and the phase/reconstructed att
    #         print("ROTATION_AXIS_POSITION = {:f} # Position in pixels".format(DS.axis_position), file=oar_file)
    #         print("", file=oar_file)
    #         print("## Parameters defining reconstruction", file=oar_file)
    #         print("OUTPUT_SINOGRAMS = NO # Output sinograms to files or not", file=oar_file)
    #         print("OUTPUT_RECONSTRUCTION = YES # Reconstruct and save or not", file=oar_file)
    #         print("START_VOXEL_1 =      1 # X-start of reconstruction volume", file=oar_file)
    #         print("START_VOXEL_2 =      1 # Y-start of reconstruction volume", file=oar_file)
    #         print("START_VOXEL_3 =      1 # Z-start of reconstruction volume", file=oar_file)
    #         print("END_VOXEL_1 =   {} # X-end of reconstruction volume".format(DS.nx), file=oar_file)
    #         print("END_VOXEL_2 =   {} # Y-end of reconstruction volume".format(DS.nx), file=oar_file)
    #         print("END_VOXEL_3 =   {} # Z-end of reconstruction volume".format(DS.ny), file=oar_file)
    #         print("OVERSAMPLING_FACTOR = 4 # 0 = Linear, 1 = Nearest pixel", file=oar_file)
    #         print("ANGLE_OFFSET = 0.000000 # Reconstruction rotation offset angle in degrees", file=oar_file)
    #         print("CACHE_KILOBYTES = 256 # Size of processor cache (L2) per processor (KBytes)", file=oar_file)
    #         print("SINOGRAM_MEGABYTES = 1000 # Maximum size of sinogram storage (megabytes)", file=oar_file)
    #         print("", file=oar_file)
    #         print("## Parameters extra features PyHST", file=oar_file)
    #         print("DO_CCD_FILTER = YES # CCD filter (spikes)", file=oar_file)
    #         print("""CCD_FILTER = "CCD_Filter" """, file=oar_file)
    #         print("""CCD_FILTER_PARA = {"threshold": 0.040000 }""", file=oar_file)
    #         print("DO_SINO_FILTER = NO # Sinogram filter (rings)", file=oar_file)
    #         print("""SINO_FILTER = "SINO_Filter" """, file=oar_file)
    #         print("ar = Numeric.ones(3216,'f')", file=oar_file)
    #         print("ar[0]=0.0", file=oar_file)
    #         print("ar[2:18]=0.0", file=oar_file)
    #         print("""SINO_FILTER_PARA = {"FILTER": ar }""", file=oar_file)
    #         print("DO_AXIS_CORRECTION = NO # Axis correction", file=oar_file)
    #         print("AXIS_CORRECTION_FILE = correct.txt", file=oar_file)
    #         print("OPTIONS= { 'padding':'E' , 'axis_to_the_center':'Y' , 'avoidhalftomo':'Y'} # Padding and position axis", file=oar_file)
    #         print("TRYEDFCONSTANTHEADER=1", file=oar_file)
    #         print("ZEROOFFMASK = 1 # Mask to zero region that is not covered by field of view", file=oar_file)
    #         print("#IGNORE_FILE = ignore_angles.txt", file=oar_file)
    #         print("ZEROCLIPVALUE = 0.001 # Minimum value of radiographs after flat / before log", file=oar_file)
    #         print("", file=oar_file)
    #         print("## Parameters defining output file / format", file=oar_file)
    #         print("OUTPUT_FILE = {}".format(output_file), file=oar_file)
    #         print("", file=oar_file)
    #         print("## Reconstruction program options", file=oar_file)
    #         print("DISPLAY_GRAPHICS = NO # No images", file=oar_file)
    #         print("PROJ_OUTPUT_FILE = {}".format(proj_output_file), file=oar_file)
    #         print("STEAM_INVERTER = 1", file=oar_file)
    #         print("NSLICESATONCE = 50", file=oar_file)
    #         pass
    #     pass
    # pass
    
    def reconstruct(self, DS, volume='phase'):
        """Tomographic reconstruction"""
        self._create_parameter_file(DS, volume)
        self._create_reconstruction_files(DS, volume)
        
        if volume == 'phase':
            prefix = DS.phase_prefix
            vol_file_name = DS.phase_vol_filename
        elif volume=='attenuation':
            prefix=DS.attenuation_prefix
            vol_file_name=DS.attenuation_vol_filename
        elif volume=='retrieved_attenuation':
            prefix=DS.retrieved_attenuation_prefix
            vol_file_name=DS.retrieved_attenuation_vol_filename
        elif volume=='update':
            prefix=DS.update_prefix
            vol_file_name=DS.update_vol_filename
        elif volume=='attenuation_update':
            prefix=DS.attenuation_update_prefix
            vol_file_name=DS.attenuation_update_vol_filename
        elif volume=='prior':
            prefix=DS.prior_forward_prefix
            vol_file_name=DS.prior_vol_filename


        if os.path.exists(vol_file_name):        
            os.remove(vol_file_name)
        if os.path.exists(prefix+"_.sig"):
            os.remove(prefix+"_.sig")

        cmd="oarsub -S "+prefix+"_.oar"
        os.system(cmd)
        timer=0
        while not os.path.exists(prefix+"_.sig"):
            if not(timer%1):
                text = "\r"+"Reconstructing " + DS.name + " " + volume + ", {} s".format(timer)      
                sys.stdout.write(text)
                sys.stdout.flush()
            time.sleep(1)
            timer += 1

        text = "\r"+"Reconstructing " + DS.name + " " + volume + ", {} s".format(timer)    
        sys.stdout.write(text)
        sys.stdout.flush()
        pass
         
    def _create_parameter_file(self, DS, volume='phase', direction='back'):
        #TODO: Combine with forward projection files
        background_file=DS.path+'/'+DS.attenuation_name+'_/dark.edf'
        flatfield_file='refHST.edf' 
        ff_prefix = '/mntdirect/_data_id19_bones01/bones1/max/Holocalibration/phantom_0.7um_20.5_1_/refHST'
        
        if direction=='back':    
            if volume == 'phase':
                prefix=DS.phase_prefix
                vol_file_name = DS.phase_vol_filename
                subtract_background='NO'
                correct_flatfield='NO'
                axis_position=DS.axis_position[DS.reference_distance-1]
                take_logarithm='NO'
            elif volume=='attenuation':
                prefix=DS.attenuation_prefix
                vol_file_name=DS.attenuation_vol_filename[0]
                subtract_background='YES'
                correct_flatfield='YES'
                take_logarithm = 'YES' # Take log of projection values
                axis_position=DS.attenuation_axis_position
            elif volume=='retrieved_attenuation':
                prefix=DS.retrieved_attenuation_prefix
                vol_file_name=DS.retrieved_attenuation_vol_filename
                subtract_background='NO'
                correct_flatfield='NO'
                take_logarithm='NO'
                axis_position=DS.axis_position[DS.reference_distance-1]
            elif volume == 'update':
                prefix=DS.update_prefix
                vol_file_name = DS.update_vol_filename
                subtract_background='NO'
                correct_flatfield='NO'
                axis_position=DS.axis_position[DS.reference_distance-1]
                take_logarithm='NO'
            elif volume=='attenuation_update':
                prefix=DS.attenuation_update_prefix
                vol_file_name=DS.attenuation_update_vol_filename
                subtract_background='NO'
                correct_flatfield='NO'
                take_logarithm = 'NO' # Take log of projection values
                axis_position=DS.axis_position[DS.reference_distance-1]

        elif direction=='forward':
            if volume == 'phase': 
                prefix=DS.phase_forward_prefix
                vol_file_name=DS.phase_vol_filename
                subtract_background='NO'
                correct_flatfield='NO'
                axis_position=DS.axis_position[DS.reference_distance-1]
                take_logarithm='NO'
                proj_output_file=DS.phase_forward_prefix
               
            elif volume=='attenuation':
                prefix=DS.attenuation_forward_prefix
                vol_file_name=DS.attenuation_vol_filename
                subtract_background='NO'
                correct_flatfield='NO'
                take_logarithm = 'YES' # Take log of projection values
                axis_position=DS.axis_position[DS.reference_distance-1] #TOOO: Verify
                proj_output_file=DS.attenuation_forward_prefix
                pass
            elif volume=='retrieved_attenuation':
                prefix=DS.retrieved_attenuation_forward_prefix
                vol_file_name=DS.retrieved_attenuation_vol_filename
                subtract_background='NO'
                correct_flatfield='NO'
                take_logarithm = 'NO' # Take log of projection values
                axis_position=DS.axis_position[DS.reference_distance-1] #TODO: Verify (should be like this, all proj to reference plane axis?)
                proj_output_file=DS.retrieved_attenuation_forward_prefix
                pass
            elif volume=='prior':
                prefix=DS.prior_forward_prefix
                vol_file_name=DS.prior_vol_filename
                subtract_background='NO'
                correct_flatfield='NO'
                take_logarithm = 'NO' # Take log of projection values
                axis_position=DS.axis_position[DS.reference_distance-1] #TODO: Verify (should be like this, all proj to reference plane axis?)
                proj_output_file=DS.prior_forward_prefix
                
                
            pass
            proj_output_file+='_%04d.edf'
        
        with open(prefix+'_.par', 'w') as oar_file:
            print("## PYHST PARAMETER FILE", file=oar_file)                 
            print("## Parameters defining the projection file series", file=oar_file)                 
            print("", file=oar_file)                 
            print("FILE_PREFIX = {}".format(prefix+'_'), file=oar_file)
            print("NUM_FIRST_IMAGE = 0 # No. of first projection file", file=oar_file)
            print("NUM_LAST_IMAGE = {} # No. of last projection file".format(DS.nbprojections-1), file=oar_file)
            print("TRYEDFCONSTANTHEADER = 0 # assume constant edf header size", file=oar_file)
            print("NUMBER_LENGTH_VARIES = NO", file=oar_file)
            print("LENGTH_OF_NUMERICAL_PART = 4 # No. of characters", file=oar_file)
            print("FILE_POSTFIX = .edf", file=oar_file)
            print("FILE_INTERVAL = 1 # Interval between input files", file=oar_file)
            print("", file=oar_file)            
            print("## Parameters defining the projection file format", file=oar_file)
            print("NUM_IMAGE_1 = {} # Number of pixels horizontally".format(DS.nx), file=oar_file)
            print("NUM_IMAGE_2 = {} # Number of pixels vertically".format(DS.ny), file=oar_file)
            print("IMAGE_PIXEL_SIZE_1 = {} # Pixel size horizontally (microns)".format(DS.pixel_size), file=oar_file)
            print("IMAGE_PIXEL_SIZE_2 = {} # Pixel size vertically".format(DS.pixel_size), file=oar_file)
            print("IMAGE_PIXEL_SIZE_2 = {} # Pixel size vertically".format(DS.pixel_size), file=oar_file)
            print("", file=oar_file)
            print("## Parameters defining background treatment", file=oar_file)
            print("SUBTRACT_BACKGROUND = {} # No background subtraction".format(subtract_background), file=oar_file)
            print("BACKGROUND_FILE = {}".format(background_file), file=oar_file)
            print("", file=oar_file)
            print("## Parameters defining flat-field treatment", file=oar_file)
            print("CORRECT_FLATFIELD = {} # No flat-field correction".format(correct_flatfield), file=oar_file)
            print("FLATFIELD_CHANGING = N.A.", file=oar_file)
            print("FLATFIELD_FILE = {}".format(flatfield_file), file=oar_file)
            print("FF_PREFIX = {}".format(ff_prefix), file=oar_file)
            print("FF_NUM_FIRST_IMAGE = 0", file=oar_file)
            print("FF_NUM_LAST_IMAGE = {}".format(DS.nbprojections), file=oar_file)
            print("FF_NUMBER_LENGTH_VARIES = NO", file=oar_file)
            print("FF_LENGTH_OF_NUMERICAL_PART = 4", file=oar_file)
            print("FF_POSTFIX = .edf", file=oar_file)
            print("TAKE_LOGARITHM = {} # Take log of projection values".format(take_logarithm), file=oar_file)
            print("FF_FILE_INTERVAL = {}".format(DS.reference_interval), file=oar_file)
            print("", file=oar_file)
            print("## Parameters defining experiment", file=oar_file)
            print("ANGLE_BETWEEN_PROJECTIONS = {} # Increment angle in degrees".format(DS.angle_increment), file=oar_file)
            print("ROTATION_VERTICAL = YES", file=oar_file)
            print("ROTATION_AXIS_POSITION = {:f} # Position in pixels".format(axis_position), file=oar_file)
            print("", file=oar_file)
            print("## Parameters defining reconstruction", file=oar_file)
            print("OUTPUT_SINOGRAMS = NO # Output sinograms to files or not", file=oar_file)
            print("OUTPUT_RECONSTRUCTION = YES # Reconstruct and save or not", file=oar_file)
            print("START_VOXEL_1 =      1 # X-start of reconstruction volume", file=oar_file)
            print("START_VOXEL_2 =      1 # Y-start of reconstruction volume", file=oar_file)
            print("START_VOXEL_3 =      1 # Z-start of reconstruction volume", file=oar_file)
            print("END_VOXEL_1 =   {} # X-end of reconstruction volume".format(DS.nx), file=oar_file)
            print("END_VOXEL_2 =   {} # Y-end of reconstruction volume".format(DS.nx), file=oar_file)
            print("END_VOXEL_3 =   {} # Z-end of reconstruction volume".format(DS.ny), file=oar_file)
            print("OVERSAMPLING_FACTOR = 4 # 0 = Linear, 1 = Nearest pixel", file=oar_file)
            print("ANGLE_OFFSET = 0.000000 # Reconstruction rotation offset angle in degrees", file=oar_file)
            print("", file=oar_file)
            print("## Parameters extra features PyHST", file=oar_file)
            print("DO_CCD_FILTER = NO # CCD filter (spikes)", file=oar_file)
            print("""CCD_FILTER = "CCD_Filter" """, file=oar_file)
            print("""CCD_FILTER_PARA = {"threshold": 0.040000 }""", file=oar_file)
            print("DO_SINO_FILTER = NO # Sinogram filter (rings)", file=oar_file)
            print("""SINO_FILTER = "SINO_Filter" """, file=oar_file)
            print("ar = Numeric.ones(3216,'f')", file=oar_file)
            print("ar[0]=0.0", file=oar_file)
            print("ar[2:18]=0.0", file=oar_file)
            print("""SINO_FILTER_PARA = {"FILTER": ar }""", file=oar_file)
            print("DO_AXIS_CORRECTION = NO # Axis correction", file=oar_file)
            print("AXIS_CORRECTION_FILE = correct.txt", file=oar_file)
            print("OPTIONS= { 'padding':'E' , 'axis_to_the_center':'Y' , 'avoidhalftomo':'N'} # Padding and position axis", file=oar_file)
            print("ZEROOFFMASK = 1 # Mask to zero region that is not covered by field of view", file=oar_file)
            print("#IGNORE_FILE = ignore_angles.txt", file=oar_file)
            print("FBFILTER = 0 # ", file=oar_file)
            print("", file=oar_file)
            print("## Parameters defining output file / format", file=oar_file)
            print("OUTPUT_FILE = {}".format(vol_file_name), file=oar_file)
            print("", file=oar_file)
            print("## Reconstruction program options", file=oar_file)
            print("DISPLAY_GRAPHICS = NO # No images", file=oar_file)
            
            if direction=='forward':            
                print("PROJ_OUTPUT_FILE = {}".format(proj_output_file), file=oar_file)
                print("STEAM_INVERTER = 1", file=oar_file)
                print("NSLICESATONCE = 50", file=oar_file)
            
            pass
        pass
    pass
