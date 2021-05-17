#This file is part of the PyPhase software.alpha
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
import os.path
import sys
import numpy as np
import yaml
import re
#from sortedcontainers import SortedList
from vendor.EdfFile import EdfFile
import pickle
from math import *
import time
from xml.dom import minidom
import matplotlib.pyplot as pyplot
import scipy.ndimage as ndimage
from pathlib import Path
import h5py
if os.name == 'posix':
    import fcntl #TODO:from importlib import reload # python 2.7 does not require this
from pyphase.config import *
import vendor.pyelastix as pyelastix


class Dataset:
    """
    Abstract class representing a set of recorded images, acquisition 
    parameters, reconstructed images and any intermediary images. 
    
    Attributes
    ----------
    
    """

    def __init__(self, name: str, path: str='.', version: str='master'):
        """
        Template for initialisation of a Dataset.
        
        Parameters
        ----------
        name : str
            Name of the dataset. Corresponds to the name of the file(s) or 
            directory(ies) of the recorded data.
        path : str, optional
            Path to the file(s) or directory(ies) of the recorded data
        version : str, optional.
            Version of dataset file. To permit tests of different alignment
            and reconstruction parameters.
    
        """
        pass

    def _initialise(self):
        self.motion_corrected = 0
        self.aligned = 0
        # TODO: Should erase files also
        self.populate()
        self._write_parameter_file()

    def align_projection(self, *, projection=0, save_projection=False):
        """
        Aligns images at different positions at one projection angle 
        
        
        Aligns (registers) images taken at different positions using the
        default registration algorithm (pyphase.registrator).
        
        Parameters
        ----------
        
        projection : int, optional
            Number of projection to register.
        save_projection : bool, optional
            Save aligned images to file (creates new data files)
        """
        
        #%% Images read should not be aligned if alignent exist. 
        # TODO: necessary?
        correct_alignment = self.correct_alignment 
        self.correct_alignment = 0
        
        #%% Read intensities 
        I_D = np.zeros((len(self.position), self.ny, self.nx))
        for N in self.position:
            I_D[N]= self.get_projection(projection=projection, position=N, aligned=False, pad=False)
        
        #%% register images
        alignment = np.zeros((len(self.position),registrator.number_of_parameters)) #TODO: check existance

        for N in self.position:
            if N != self.reference_position:
                print('Aligning position '+ str(N))
                field, im1reg, transform_parameters = registrator.register(I_D[N], I_D[self.reference_position])
                alignment[N] = transform_parameters
            else:
                im1reg = I_D[self.reference_position] #TODO: Maybe not the safest to take im1reg but rather transform image?
            if save_projection:         #TODO: This should be generalised somehow I suspect
                aligned_filename = self.path + '/' + self.name + '_aligned' + '_dist' + str(N) + '.h5'
                with h5py.File(aligned_filename, 'a') as f:
                    if not 'data' in f.keys():
                        f.create_group(self.cpr_prefix)
                    if str(projection).zfill(6) in f[self.cpr_prefix].keys():
                        f[self.cpr_prefix+str(projection).zfill(6)][:] = im1reg
                    else:
                        f.create_dataset(self.cpr_prefix+str(projection).zfill(6), data=im1reg)
                    
        #%% Store alignment
        #TODO: Check for length of trasnform parameterrs
        #TODO: Quick and dirty implementation wait for write access parallel
        while True:
            try:
                with h5py.File(self.dataset_filename, 'a') as f:
                    if 'alignment' in f.keys():
                        f['alignment'][projection] = alignment
                    else:
                        f.create_dataset('alignment', data=np.zeros((self.number_of_projections, len(self.position), registrator.number_of_parameters)))
                        f['alignment'][projection] = alignment
                break  # Success!
            except OSError:
                time.sleep(1)  # Wait a bit

        self.correct_alignment = correct_alignment

    def get_alignment(self, *, projection, position):
        """
        Get transform parameters for alignment.
        
        Returns
        -------
        transform_parameters : numpy array
            Array of transform parameters. Size depends on transform used.
        
        """
        with h5py.File(self.dataset_filename, 'r') as f:    
            transform_parameters = f['alignment'][projection, position, :]
        return transform_parameters

    def write_image(self, *, image, projection, projection_type='phase'):
        # TODO: raise error where needed! number of args, args type...
        ''' Saves images into file.

        Parameters
        ----------
        image : ndarray
           An ndarray of the image to save.
        projection : int
           The number of the projection to be saved.
        proj_type : str
           A string containing the prefix for the name of the file to be saved.
        position: int
           Number of the position of the projection (for propagation).
        Returns
        -------
        None
           Saves the images into the file ``[prefix]_[projection_number].edf`` or
           ``[prefix]_[distance]_[projection_number].edf`` into the ``name_`` directory.
        '''
        if projection_type.lower() == 'phase':
            prefix = self.phase_prefix          
        elif projection_type.lower() == 'attenuation':
            prefix =  self.attenuation_prefix
        # elif proj_type.lower() == 'phase update':
        #     pref = self.update_prefix
        # elif proj_type.lower() == 'attenuation update':
        #     pref = self.attenuation_update_prefix
        # elif proj_type.lower() == 'intensity':
        #     pref = self.propagated_prefix
        else:
            raise TypeError("proj_type must be one of the following: 'phase', 'attenuation', 'phase update', 'attenuation update', 'intensity'" )
        
        projection_name = str(projection).zfill(6)

        # TODO: Initial file locking implementation for parallelisation on LINUX/SLURM
        lock_filename = self.dataset_filename.with_suffix(self.dataset_filename.suffix+'.lock')
        try:
            lock_file = open(lock_filename, 'a') # Create a lock file next to dataset file (TODO: Consider making invisible at least)
            if os.name == 'posix': #TODO: quick fix for lock file. Shouldn be necessary on windows due to no such parallelism? WHat about future multithread?
                fcntl.lockf(lock_file, fcntl.LOCK_EX) # Acquire a lockf lock (holds inter-process but not intra?)
            with h5py.File(self.dataset_filename, 'a') as f:    
                data = f.require_dataset(prefix+projection_name, shape=image.shape, dtype=image.dtype)
                data[:] = image   
            lock_file.close() # Release lock
        except OSError:
            time.sleep(1)  # Wait a bit

    @Parallelize # TODO:How to get to reload from config?
    def align_projections(self, *, projections):
        """
        Align a series of projections.
        
        Parameters
        ----------
        projections : tuple of ints
            In the form [start, end]
        """
        for projection in range(projections[0], projections[1]+1):
            self.align_projection(projection=projection)

    def get_image(self, *, projection, image_type='phase', Fourier=False, pad=False):
        """Get reconstructed images (phase by default)

        Parameters
        ----------
        projection : int
           Number of the projection.
        image_type : str, optional
            Can be phase or attenuation.
        Fourier : bool
            Return Fourier transform of image.
        pad : bool
            Return padded image.

        Returns
        -------
        image : ndarray
        """
        
        #TODO: handling majuscule minuscule
        # if 'distance' in kwargs:
        #     distance = kwargs.get('distance') #TODO: Needs error handling
        if image_type == 'attenuation':
            prefix = self.attenuation_prefix
        # elif 'propagated' in argv:        
        #     fname = self.propagated_prefix.as_posix()+'_'+str(distance)+'_'+str(projection).zfill(4)+'.edf'
        # elif 'difference' in argv:
        #     fname = self.difference_prefix.as_posix()+'_'+str(distance)+'_'+str(projection).zfill(4)+'.edf'
        # elif 'prior' in argv:
        #     fname = self.prior_forward_prefix.as_posix()+('_'+str(projection).zfill(4)+'.edf')
        else:
            prefix = self.phase_prefix #TODO: separate file type...

        with h5py.File(self.dataset_filename, 'r') as f:    
            image=f[prefix+str(projection).zfill(6)][:]
                        
        if Fourier:
            #TODO: Refactor into an FT function in utilities (f ex)
            return np.fft.fft2(Utilities.resize(image, (self.nfx,self.nfy)))
        elif pad:
            image = Utilities.resize(image, (self.nfx,self.nfy))
        return image

class Nanomax(Dataset):
    """Class for raw Nanomax data sets.
    
    Attributes
    ----------
    ztot : float
        Focus to detector distance in m.
    path : str
        Path to dataset.
    version : str
        Version of Dataset backend structure, e.g. 'TIEHOM'
    name : str
        Name of dataset
    frames_per_position : int
        Number of frames per position.
    initial_frame : int
        Number of the initial frame.
    projection_prefix : str
        HDF5 path to projection data.
    reference_position : int
        Position number as reference for alignment (usually highest resolution)
    diode_name : str
        Name of diode (HDF path) for normalisation
    darkfield_prefix : str
        HDF5 path prefix to darkfields
    flatfield_prefix : str
        HDF5 path prefix to flatfields
    cpr_prefix : str
        HDF5 path prefix to flatfield corrected images
    phase_prefix : str
        HDF5 path prefix to phase maps
    attenuation_prefix : str
        HDF5 path prefix to attenuation images
    aligned : int
        Flag if dataset is aligned or not.
    """
    
    def __init__(self, name: str,path: str='.', version: str= 'ḿaster'):        
        self.ztot = 1.123 # Focus to detector distance in m
        self.path = path
        self.version = version
        self.name = name
                         
        self.data_basename = self.path + '/' + self.name + '/' #TODO: Attribute?
        self.dataset_filename = pyphase_path / (self.name + '_' + self.version + '.h5') #TODO: Attribute?

        self.frames_per_position = 4
        self.initial_frame = 0

        self.projection_prefix = 'entry/measurement/andor/frames/'

        self.reference_position = 0
        self.correct_alignment = 1
        
        self.measurement_prefix = 'entry/snapshot/'
        self.pixel_size_filename = path + '/' 'pixelsizes.txt'
        self.diode_name = 'entry/measurement/alba2/1'
        self.darkfield_prefix = 'reference_images/darkfield/'
        self.flatfield_prefix = 'reference_images/flatfield/'
        self.cpr_prefix = 'data/cpr/'
        self.phase_prefix = '/phase/'
        self.attenuation_prefix = '/attenuation/'
        
        self.aligned=0 # TODO: Should be bool
        
        super().__init__(name, path, version)

        #TODO: Quick and dirty concession to parallel requirements
        while True:
            try:
                with h5py.File(self.dataset_filename, 'a') as dataset_file:    
                    if 'scan_parameters' in dataset_file.keys():
                        #TODO: initialise must override reading of parameter file
                        #self.ReadParameters()
                        data = dataset_file['scan_parameters']
                        self.energy = data['energy'][()]
                        self.effective_distance = data['effective_distance'][()]
                        self.pixel_size_x = data['pixel_size_x'][()]
                        self.pixel_size_y = data['pixel_size_y'][()]
                        self.number_of_projections = data['number_of_projections'][()]
                        self.nx = data['nx'][()]
                        self.ny = data['ny'][()]
                    else:
                        #self.Initialise()
                        #TODO: STub. Should take values from the data h5 file
                        with open(self.pixel_size_filename) as pixel_size_file:
                            lines = pixel_size_file.readlines()
                
                        del lines[0] # remove header line
                        self.pixel_size_x = np.zeros(len(lines))
                        self.pixel_size_y = np.zeros_like(self.pixel_size_x)
                        z1 = np.zeros_like(self.pixel_size_x) # Focus to sample distance (given in um in text files)
                
                        for n, line in enumerate(lines):
                            line = line.split()
                            z1[n] = float(line[0])/1e6 # Focus to sample distance (given in um in text files)
                            self.pixel_size_x[n] = float(line[1])
                            self.pixel_size_y[n] = float(line[2])
                            pass
        
                        z2 = self.ztot - z1 # Sample to detector distance
                        self.effective_distance =  z1*z2 / (z1+z2)
         
                        # Read parameters from data file
                        projection_data_filename = self.data_basename + str(self.initial_frame+2).zfill(6) + '.h5'
                        with h5py.File(projection_data_filename, 'r') as data_file:
                            self.energy = float(data_file[self.measurement_prefix+'energy'][:]) / 1e3 # eV -> keV
                            self.number_of_projections = len(data_file[self.projection_prefix])
                            self.nx = data_file[self.projection_prefix][:].shape[2]
                            self.ny = data_file[self.projection_prefix][:].shape[1]
                            
                        # Calculate darks & flats
                        print('Calculating darkfield and flatfield images')
        
                        for position in range(len(self.effective_distance)):
                            file_base_number = self.initial_frame + position*self.frames_per_position
                            
                            # Calculate darks
                            dark_data_filename = self.data_basename + str(file_base_number).zfill(6) + '.h5'
                            with h5py.File(dark_data_filename, 'r') as data_file:
                                darkfields = data_file[self.projection_prefix][:]
                                weights = data_file[self.diode_name][:] # Read diode measurements
                                darkfield = np.average(darkfields,axis=0,weights=weights)
                                dataset_dark = dataset_file.require_dataset(self.darkfield_prefix + str(position).zfill(6), darkfield.shape, darkfield.dtype)
                                dataset_dark[:] = darkfield
        
                            # Calculate flats
                            # TODO: not cleanly implemented wrt intermediate flats (but YAGNI)
                            
                            flat_beginning_data_filename = self.data_basename + str(file_base_number+1).zfill(6) + '.h5' 
                            flatfield = np.zeros((2, self.ny, self.nx)) # TODO: Adapt to flexible number of flats (but YAGNI)
                            with h5py.File(flat_beginning_data_filename, 'r') as data_file:
                                flatfields = data_file[self.projection_prefix][:]
                                weights = data_file[self.diode_name][:]
                                flatfield[0] = np.average(flatfields, axis=0, weights=weights)
                                
                            flat_end_data_filename = self.data_basename + str(file_base_number+3).zfill(6) + '.h5' 
                            with h5py.File(flat_end_data_filename, 'r') as data_file:
                                flatfields = data_file[self.projection_prefix][:]
                                weights = data_file[self.diode_name][:]
                                flatfield[1] = np.average(flatfields, axis=0, weights=weights)
                                dataset_flat = dataset_file.require_dataset(self.flatfield_prefix + str(position).zfill(6), flatfield.shape, flatfield.dtype)
                                dataset_flat[:] = flatfield
        
                        g = dataset_file.create_group('scan_parameters')
                        g.create_dataset('energy', data=self.energy)
                        g.create_dataset('effective_distance', data=self.effective_distance)
                        g.create_dataset('pixel_size_x', data=self.pixel_size_x)
                        g.create_dataset('pixel_size_y', data=self.pixel_size_y)
                        g.create_dataset('number_of_projections', data=self.number_of_projections)
                        g.create_dataset('nx', data=self.nx)
                        g.create_dataset('ny', data=self.ny)
                        if 'alignment' in dataset_file.keys():
                            self.aligned = 1
                        pass #With
                break  # Success!
            except OSError:
                time.sleep(1)  # Wait a bit

        self.padding = 2

        self.position=np.arange(len(self.effective_distance))

        self.magnification_x = self.pixel_size_x[0] / self.pixel_size_x
        self.magnification_y = self.pixel_size_x[0] / self.pixel_size_x
        self.pixel_size = np.array([self.pixel_size_x[0], self.pixel_size_y[0]])

        self.reference_interval = self.number_of_projections


    @property
    def nfx(self):
        return self.padding * self.nx
    
    @property
    def nfy(self):
        return self.padding * self.ny
        
    def get_projection(self, *, projection, position, pad=True, magnification=True, aligned=True, Fourier=False):
        """
        Read one recorded image.
        
        Parameters
        ----------
        projection : int
            Number of projection to read.
        position : int
            Number of positiion ("distance") to read.
        pad : bool, optional
            Pads the image.
        magnification : bool, optional
            Brings the image to the magnification of reference_position.
        aligned : bool, optional
            Corrects alignment (requires alignment of projections).
        Fourier : bool
            Returns the Fourier transform of the image.
            
        """
        
        # TODO: consider renaming just projection() to get data.projection(1,23)
        
        # Nanomax data: frame inital_frame=dark, initial_frame+1=flat_beginning, initial_frame+2=projections, initial_frame+3=flat_end
        
        projection_name = str(projection).zfill(6)
        position_name = str(position).zfill(6)
        
        #TODO: Quick and nasty concession to parallel computing
        
        while True:
            try:
                with h5py.File(self.dataset_filename, 'r') as f:    
                    dark_image = f[self.darkfield_prefix+position_name][:] # Read projection
                    flatfield_image = f[self.flatfield_prefix+position_name][:] # Read projection
                break  # Success!
            except OSError:
                time.sleep(1)  # Wait a bit

        file_base_number = self.initial_frame + position*self.frames_per_position + 2
        #TODO: Quick and nasty concession to parallel computing
        while True:
            try:
                with h5py.File(self.data_basename + str(file_base_number).zfill(6) + '.h5', 'r') as f:
                    projection_image = f[self.projection_prefix][projection] # Read projection
                break  # Success!
            except OSError:
                time.sleep(1)  # Wait a bit

        reference_index_1 = (self.reference_interval * (projection // self.reference_interval))
        reference_index_2 = min(reference_index_1+self.reference_interval, self.number_of_projections)
        
        # Weights flat fields
        weight_1 = (reference_index_2-projection)/max(reference_index_2-reference_index_1, 1)
        weight_2 = 1 - weight_1                      
        
        # Calculate flat field
        reference_image = weight_1*flatfield_image[0] + weight_2*flatfield_image[1]

        # Calculate flat and dark field corrected projecion                                  
        projection_image = (projection_image-dark_image) / (reference_image-dark_image)

        # Zoom images to same scale, corresponding to highest magnification
        if position != self.reference_position and magnification: #TODO: should implement distance_number, but YAGN?
            projection_image = ndimage.zoom(projection_image, [1/self.magnification_y[position], 1/self.magnification_x[position]], mode='nearest') # Zoom to highest magnification TODO: combine w/ alignment?     
        
        # Correct alignment
        if aligned and self.aligned and position != self.reference_position:
            transform_parameters = self.get_alignment(projection=projection, position=position) # Get alignment parameters TODO: add zoom factor here?     
            projection_image = registrator.apply_transformation(projection_image, transform_parameters)
            
        if Fourier:
            projection_image = Utilities.resize(projection_image, (self.nfy, self.nfx))
            return np.fft.fft2(projection_image)
        else:
            if pad:
                return Utilities.resize(projection_image, (self.nfy, self.nfx))
            else:
                return Utilities.resize(projection_image, (self.ny, self.nx))
        
class NanomaxPreprocessed(Dataset):
    """
    Class for preprocessed Nanomax datasets.
    
    Attributes
    ----------
    ztot : float
        Focus to detector distance in m.
    path : str
        Path to dataset.
    version : str
        Version of Dataset backend structure, e.g. 'TIEHOM'
    name : str
        Name of dataset
    correct_alignent : int
        Flag to correct alignment.
    aligned: int
        Flag for dataset aligned.
    phase_prefix : str
        HDF5 path prefix to phase maps
    attenuation_prefix : str
        HDF5 path prefix to attenuation images
    projection_prefix : str
        HDF5 path to corrected projections
    reference_position : int
        Position number as reference for alignment (usually highest resolution)
    energy : float
        Energy in keV
    nx, ny : int
        Number of pixels, horizontal, vertical
    padding : int
        Padding factor
    detector_pixel_size : float
        The detector pixel size in µm
    alpha : tuple of floats
        Regularization parameter in the form [LF, HF]   
    """
 
    def __init__(self, name: str, path: str='.', version: str='master'):
        self.ztot = 1.123 # Focus to detector distance in m
        self.path = path
        self.version = version
        self.name = name
        
        self.correct_alignment = 1 #TODO: should be bool
        self.aligned = 0 # TODO: inelegant solution, and should be bool
        # TODO: Pre-defined parameters for the moment. Should be stored in h5 files
        self.projection_prefix = '/data/cpr/'
        self.phase_prefix = '/phase/'
        self.attenuation_prefix = '/attenuation/'
        self.cpr_prefix = '/data/cpr/'
        self.reference_position = 0
        self.energy = 13
        self.nx = 2048 
        self.ny = 2048
        self.padding=2
        self.detector_pixel_size = 0
        self.alpha = [-8, -10] #TODO: HOw should this be handled really? The reconstruction parameters (alg, params, ... should definitely be in dataset, associated to a version even!)

    @property
    def nfx(self):
        """Number of pixels in Fourier domain, horizontal, calculated from 
        nx and padding"""

        return self.padding * self.nx
    
    @property
    def nfy(self):
        """Number of pixels in Fourier domain, vertical, calculated from 
        ny and padding"""

        return self.padding * self.ny

    @property
    def Lambda(self):
        """Wavelength in m, calculated from energy."""
        return 12.4e-10 / self.energy

class NanomaxPreprocessedTomo(NanomaxPreprocessed):
    """
    Class for preprocessed Nanomax tomographic datasets.


    Attributes
    ----------
    data_basename : str
        File prefix to projection data files.
    magnification_x,y : numpy array
        Magnification factors for each position, horizontal, vertical
    pixel_size_x,y : numpy array
        Effective pixel size for each position
    pixel_size : numpy array
         Effective pixel size at reconstruction position  
    
    """

    def __init__(self, name: str, path: str='.', version: str='master'):
        #TODO: Pixel sizes must be calculated on the fly, but not possible at the moment
        
        self.path = path
        self.version = version
        self.name = name
                         
        self.data_basename = self.path + '/' + self.name + '_dist'
        self.dataset_filename = pyphase_path / (self.name + '_' + self.version + '.h5')

        super().__init__(name, path, version)

        with h5py.File(self.dataset_filename, 'a') as f:    
            if 'scan_parameters' in f.keys():
                #TODO: initialise must override reading of parameter file
                #self.ReadParameters()
                data = f['scan_parameters']
                self.energy = data['energy'][()]
                self.effective_distance = data['effective_distance'][()]
                self.pixel_size_x = data['pixel_size_x'][()]
                self.pixel_size_y = data['pixel_size_y'][()]
                self.number_of_projections = data['number_of_projections'][()]
            else:
                #self.Initialise()
                #TODO: Stub. Should take values from the data h5 file
                self.pixel_size_filename = path + '/' 'pixelsizes.txt'
                with open(self.pixel_size_filename) as pixel_size_file:
                    lines = pixel_size_file.readlines()
        
                del lines[0] # remove header line
                self.pixel_size_x = np.zeros(len(lines))
                self.pixel_size_y = np.zeros_like(self.pixel_size_x)
                z1 = np.zeros_like(self.pixel_size_x) # Focus to sample distance (given in um in text files)
        
                for n, line in enumerate(lines):
                    line = line.split()
                    z1[n] = float(line[0])/1e6 # Focus to sample distance (given in um in text files)
                    self.pixel_size_x[n] = float(line[1])
                    self.pixel_size_y[n] = float(line[2])
                    pass

                z2 = self.ztot - z1 # Sample to detector distance
                self.effective_distance =  z1*z2 / (z1+z2)

                g = f.create_group('scan_parameters')
                g.create_dataset('energy', data=self.energy)
                g.create_dataset('effective_distance', data=self.effective_distance)
                g.create_dataset('pixel_size_x', data=self.pixel_size_x)
                g.create_dataset('pixel_size_y', data=self.pixel_size_y)
                with h5py.File(self.data_basename+str(0)+'.h5', 'r') as data_file:
                    self.number_of_projections = len(data_file[self.projection_prefix])
                g.create_dataset('number_of_projections', data=self.number_of_projections)
                
            if 'alignment' in f.keys():
                self.aligned = 1

        self.position=np.arange(len(self.effective_distance))

        self.magnification_x = self.pixel_size_x[0] / self.pixel_size_x
        self.magnification_y = self.pixel_size_x[0] / self.pixel_size_x
        self.pixel_size = np.array([self.pixel_size_x[0], self.pixel_size_y[0]])

    def get_projection(self, *, projection, position, pad=True, magnification=True, aligned=True, Fourier=False):
        """
        Read one recorded image.
        
        Parameters
        ----------
        projection : int
            Number of projection to read.
        position : int
            Number of positiion ("distance") to read.
        pad : bool, optional
            Pads the image.
        magnification : bool, optional
            Brings the image to the magnification of reference_position.
        aligned : bool, optional
            Corrects alignment (requires alignment of projections).
        Fourier : bool
            Returns the Fourier transform of the image.
            
        """

        #TODO: consider renaming just projection() to get data.projection(1,23)
        
        data_filename = self.data_basename + str(position) + '.h5'
        projection_name = str(projection).zfill(6)
        
        with h5py.File(data_filename, 'r') as f:    
            projection_image = f[self.projection_prefix+projection_name][:] # Read projection
        
        if position != self.reference_position and magnification: #TODO: should implement distance_number, but YAGN?
            projection_image = ndimage.zoom(projection_image, [1/self.magnification_y[position], 1/self.magnification_x[position]], mode='nearest') # Zoom to highest magnification TODO: combine w/ alignment?     
        
        if aligned and self.aligned and position != self.reference_position:
          
            # TODO: Currently using similarity transform. It would be nice to be able to change transforms:
            # TODO: Translation, Euler, Similarity, Affine... and handle the different length of parameter vector
            transform_parameters = self.get_alignment(projection=projection, position=position) # Get alignment parameters TODO: add zoom factor here?     
            #projection_image = ndimage.zoom(projection_image, magnification_factor/transform_parameters[0])
            #projection_image = ndimage.rotate(projection_image, np.degrees(transform_parameters[1]))
            #projection_image = ndimage.shift(projection_image, [-transform_parameters[3], -transform_parameters[2]])
            projection_image = registrator.apply_transformation(projection_image, transform_parameters)
            
        if Fourier:
            projection_image = Utilities.resize(projection_image, (self.nfy, self.nfx))
            return np.fft.fft2(projection_image)
        else:
            if pad:
                return Utilities.resize(projection_image, (self.nfy, self.nfx))
            else:
                return Utilities.resize(projection_image, (self.ny, self.nx))
             
class NanomaxPreprocessed2D(NanomaxPreprocessed):
    """
    Class for single projection Nanomax data.
    """

    def __init__(self, name: str, path: str='.', version: str='master'):        
        self.number_of_projections = 1 # Single projection data set

        super().__init__(name, path, version)

        # Define file access parameters
        self.data_filename = self.path + '/' + self.name + '.h5'
        self.dataset_filename = pyphase_path / (self.name + '_' + self.version + '.h5')             
        self.projection_prefix = '/data/cpr/'
        self.phase_prefix = '/phase/'
        self.attenuation_prefix = '/attenuation/'

        # Read scan parameters
        with h5py.File(self.dataset_filename, 'a') as f:    
            # If the dataset is already created, read from the file
            if 'scan_parameters' in f.keys():
                #TODO: initialise must override reading of parameter file
                #self.ReadParameters()
                data = f['scan_parameters']
                self.energy = data['energy'][()]
                self.effective_distance = data['effective_distance'][()]
                self.pixel_size_x = data['pixel_size_x'][()]
                self.pixel_size_y = data['pixel_size_y'][()]
                self.number_of_projections = data['number_of_projections'][()]
            else:
                #Else, read from the pixelsizes.txt (temporary solution? To be defined with Nanomax)
                #TODO: Stub. Should take values from the data h5 file
                self.pixel_size_filename = path + '/' 'pixelsizes.txt'
                with open(self.pixel_size_filename) as pixel_size_file:
                    lines = pixel_size_file.readlines()
        
                del lines[0] # remove header line
                self.pixel_size_x = np.zeros(len(lines))
                self.pixel_size_y = np.zeros_like(self.pixel_size_x)
                z1 = np.zeros_like(self.pixel_size_x) # Focus to sample distance (given in um in text files)
        
                for n, line in enumerate(lines):
                    line = line.split()
                    z1[n] = float(line[0])/1e6 # Focus to sample distance (given in um in text files)
                    self.pixel_size_x[n] = float(line[1])
                    self.pixel_size_y[n] = float(line[2])
                    pass
                
                # Calculate effective distances
                z2 = self.ztot - z1 # Sample to detector distance
                self.effective_distance =  z1*z2 / (z1+z2)
                
                # Create the hdf5 data
                g = f.create_group('scan_parameters')
                g.create_dataset('energy', data=self.energy)
                g.create_dataset('effective_distance', data=self.effective_distance)
                g.create_dataset('pixel_size_x', data=self.pixel_size_x)
                g.create_dataset('pixel_size_y', data=self.pixel_size_y)
                g.create_dataset('number_of_projections', data=self.number_of_projections)
        
            if 'alignment' in f.keys():
                self.aligned = 1

        # Calculate magnifications and effective pixel size
        self.magnification_x = self.pixel_size_x[0] / self.pixel_size_x
        self.magnification_y = self.pixel_size_x[0] / self.pixel_size_x
        self.pixel_size = np.array([self.pixel_size_x[0], self.pixel_size_y[0]])
        self.position=np.arange(len(self.effective_distance))
        
        pass # __init__

    def get_projection(self, *, projection, position, pad=True, magnification=True, aligned=True, Fourier=False):
        """
        Read one recorded image.
        
        Parameters
        ----------
        projection : int
            Number of projection to read.
        position : int
            Number of positiion ("distance") to read.
        pad : bool, optional
            Pads the image.
        magnification : bool, optional
            Brings the image to the magnification of reference_position.
        aligned : bool, optional
            Corrects alignment (requires alignment of projections).
        Fourier : bool
            Returns the Fourier transform of the image.
            
        """
        
        #TODO: consider renaming just projection() to get data.projection(1,23)
        #TODO: Possibility non-corrected?
        
        projection_name = str(position).zfill(6) # Constitute projection name in Nanomax dataset
        
        with h5py.File(self.data_filename, 'r') as f:    
            projection_image = f[self.projection_prefix+projection_name][:] # Read projection
        
        if position != self.reference_position and magnification: #TODO: should implement distance_number, but YAGN?
            projection_image = ndimage.zoom(projection_image, [1/self.magnification_y[position], 1/self.magnification_x[position]], mode='nearest') # Zoom to reference_plane magnification TODO: combine w/ alignment?     
        
        if aligned and position != self.reference_position: # COrrect alignment
            # TODO: Currently using similarity transform. It would be nice to be able to change transforms
            # TODO: Translation, Euler, Similarity, Affine... and handle the different length of parameter vector
            transform_parameters = self.get_alignment(projection=projection, position=position) # Get alignment parameters TODO: add zoom factor here?     
            #projection_image = ndimage.zoom(projection_image, magnification_factor/transform_parameters[0])
            #projection_image = ndimage.rotate(projection_image, np.degrees(transform_parameters[1]))
            #projection_image = ndimage.shift(projection_image, [-transform_parameters[3], -transform_parameters[2]])
            projection_image = registrator.apply_transformation(projection_image, transform_parameters)
            
        if Fourier: # Pad and Fourier transform
            projection_image = Utilities.resize(projection_image, (self.nfy, self.nfx))
            return np.fft.fft2(projection_image)
        else:
            if pad: # Pad
                return Utilities.resize(projection_image, (self.nfy, self.nfx))
            else:
                return Utilities.resize(projection_image, (self.ny, self.nx))

class ESRF():
    """Legacy code for ESRF datasets.
    
    Note
    ----
    Will be aligned with Nanomax functionality.
    
    """
    def __init__(self, name: str, path: str='.', version: str='master'):
        self.registrator = Utilities.ElastixRigid()
        self.path=Path(path) # TODO: OR shoud it take a path object already?
        self.name=name  # name of sample
        self.version=version
        self.axis_position = np.array([])
        self.position = np.array([]) # in m
        self.distance_source_sample= 0
        self.energy = 0 # in keV
        self.nD = 0
        self.nbprojections = 0
        self.position_number = []
        self.nx = 0
        self.ny = 0
        self.detector_pixel_size = 0
        self.reference_interval = 0
        self.reference_position = 0
        self.scan_range = 0
        self.delta_beta = 0       
        #TODO: Some of these things should surely go in some kind of configuration file
        self.phase_tag='phase'
        self.delta_tag='delta'
        self.beta_tag='beta'
        self.measured_beta_tag='measured_beta'
        self.forward_tag='forward'
        self.attenuation_tag='attenuation'
        self.measured_attenuation_tag='measured_attenuation'
        self.propagated_tag='propagated'
        self.update_tag='update'
        self.prior_tag='prior'
        self.attenuation_update_tag='attenuation_update'
        self.difference_tag='difference'    
        self.parameter_filename='parameters.yml'
        self.dataset_parameter_filename = self.path / pyphase_path / name / self.version / self.parameter_filename
        if (self.dataset_parameter_filename.exists()):
            #TODO: initialise must override reading of parameter file
            self._read_parameter_file()
        else:
            self.dataset_parameter_filename.parents[0].mkdir(parents=True)
            self.dataset_parameter_filename.touch()
            self._initialise()
        self.preprocessing=1        
        self.correct_alignment=1
        self.correct_motion=1
        self.padding=2 #TODO: should be a parameter/in parameter file (in each version?)
        self.magnification = 1 #TODO: should be a parameter/in parameter file (in each version?)        
        #  Pixel size & magnification (from parameter file)
        self.magnification_factor = (self.distance_source_sample+self.position*1e-3) / self.distance_source_sample
        self.max_magnification_factor = max(self.magnification_factor)
        self.pixel_size = self.detector_pixel_size / self.max_magnification_factor               
        # File access
        # Basix paths
        self.data_prefix = self.path / self.name
        self.result_prefix = self.path / pyphase_path / self.name / self.version
#        self.vol_prefix = self.path+'/volfloat/' #TODO: Should change, volumes should be stored either directly under version path or in their own directories

        # Phase
        self.phase_name = self.name+'_'+self.version # TODO necessary? Good to mark data set & version in file names?
        self.phase_prefix = self.result_prefix / self.phase_tag / self.phase_name

        # Delta
        volume_suffix = '.vol' #TODO: This will be dependent on the internal format
        self.delta_filename = (self.result_prefix / self.delta_tag / self.phase_name).with_suffix(volume_suffix)

        # Projected Delta
        self.phase_forward_name = self.phase_name+'_'+self.delta_tag+'_'+self.forward_tag # Projection of delta
        self.phase_forward_prefix = self.result_prefix / (self.delta_tag+'_'+self.forward_tag) / self.phase_forward_name

        # Measured attenuation / shortest distance
        self.measured_attenuation_name = self.name+'_'+str(self.position_number[0])
        self.measured_attenuation_prefix = self.path / self.measured_attenuation_name / self.measured_attenuation_name
        
        # measured beta
        self.measured_beta_filename = (self.result_prefix / self.measured_attenuation_name).with_suffix(volume_suffix)
        
        # Projected measured beta
        self.measured_attenuation_forward_name = self.measured_attenuation_name+'_'+self.forward_tag
        self.measured_attenuation_forward_prefix = self.result_prefix / (self.measured_beta_tag+'_'+self.forward_tag) / self.measured_attenuation_forward_name

        # retrieved attenuation 
        self.attenuation_name = self.phase_name+'_'+self.attenuation_tag
        self.attenuation_prefix = self.result_prefix / self.attenuation_tag / self.attenuation_name
        
        # reconstructed beta
        self.beta_filename = (self.result_prefix / self.attenuation_name).with_suffix(volume_suffix)

        # projected reconstructed beta
        self.attenuation_forward_name = self.attenuation_name+'_'+self.forward_tag
        self.attenuation_forward_prefix = self.result_prefix / (self.beta_tag+'_'+self.forward_tag) / self.attenuation_forward_name

        # prior
        self.prior_name=self.phase_name+'_'+self.prior_tag
        self.prior_prefix= self.result_prefix / self.prior_name

        self.prior_vol_filename = (self.result_prefix / self.prior_tag / self.prior_name).with_suffix(volume_suffix)

        self.prior_forward_name = self.prior_name+'_'+self.forward_tag
        self.prior_forward_prefix = self.result_prefix / (self.prior_tag+'_'+self.forward_tag) / self.prior_forward_name
        
        # Propagated, updates
        self.propagated_prefix=self.result_prefix / self.propagated_tag / (self.phase_name+'_'+self.propagated_tag)
        
        self.update_name = self.phase_name+'_'+self.update_tag
        self.update_prefix = self.result_prefix / self.update_tag /  self.update_name
        self.update_vol_filename = (self.result_prefix / self.update_tag / self.update_name).with_suffix(volume_suffix)
        
        self.attenuation_update_name=self.phase_name+'_'+self.attenuation_update_tag
        self.attenuation_update_prefix= self.result_prefix / self.attenuation_update_tag / self.attenuation_update_name
        self.attenuation_update_vol_filename = (self.result_prefix / self.attenuation_update_tag / self.attenuation_update_name).with_suffix(volume_suffix)

        self.difference_name = self.phase_name + self.difference_tag
        self.difference_prefix = self.result_prefix / self.difference_tag / self.difference_name 

        self.lengthscale=10e-6 #TODO: Dataset should surely not have lengthscale?
        # TODO: there should be the possibility to use different parameter files ofc

        # TODO: Add argument for parameter file\
        # TODO: If no parameter file is given, check if the standard one exists, otherwise populate from info files
 
        #
        self.angle_increment=self.scan_range/self.nbprojections
        self.Lambda = 12.4e-10 / self.energy # wavelength (in m?)
        self.effective_distance = (self.distance_source_sample * self.position) / (self.distance_source_sample + self.position)
        self.Fresnel_number = self.Lambda * self.effective_distance / (self.lengthscale ** 2)
        self.position_number_no_reference_position = self.position_number.copy()
        self.position_number_no_reference_position.remove(self.position_number[self.reference_position])
        
        # Check if shifts are done (i.e. file shifts.pickle exists)
        self.shift_filename = self.result_prefix / 'shifts.pickle'
        self.shift_polynomial_x=[0 for x in range (self.nD)]
        self.shift_polynomial_y=[0 for x in range (self.nD)]
        self.shift_polynomial_y=[0 for x in range (self.nD)]
        self.shift_polynomial_order=[4, 2, 2]
        
        if os.path.isfile(self.shift_filename):
            self.aligned=1
            with open(self.shift_filename, 'rb') as f:
                 self.shifts=pickle.load(f)
            self.fit_alignment()
        else:
            self.shifts=SortedList() # entry: []
            self.aligned=0
        
        self.motion_filename=path+'/'+name+'_/motion.pickle'
        
        if os.path.isfile(self.motion_filename):
            self.motion_corrected=1
            with open(self.motion_filename, 'rb') as f:
                self.motion_correlations = pickle.load(f)
        else:
            self.motion_correlations = [ [] for i in range(self.nD) ]
            self.motion_corrected=0

    def _initialise(self):
        self.motion_corrected = 0
        self.aligned = 0
        # TODO: Should erase files also
        self.populate()
        self._write_parameter_file()

    @property
    def nfx(self):
        return self.padding * self.nx
    
    @property
    def nfy(self):
        return self.padding * self.ny

    def get_projection(self, *, projection, position, pad=True, Fourier=False, aligned=True, magnification=True):
        #TODO: consider renaming just projection() to get data.projection(1,23)
        # For integration with rest, necessary(?) conversion from position 0 to 1
        
        projection_prefix = self.path / (self.name + '_' + str(self.position_number[position]) + '_/')
        projection_filename = projection_prefix / (self.name + '_' + str(self.position_number[position]) + '_' + str(projection).zfill(4) + '.edf')
        projection_image = EdfFile(projection_filename).GetData(0)

        xpad = (self.nfx-self.nx)//2
        ypad = (self.nfy-self.ny)//2
        
        if self.preprocessing:
            dark_filename = projection_prefix / 'dark.edf'
            #print(dark_filename)
            dark_image = EdfFile(dark_filename).GetData(0)
            # Get previous ff
            reference_index_1 = (self.reference_interval * (projection // self.reference_interval))
            reference_index_2 = min(reference_index_1+self.reference_interval, self.nbprojections)
            # Weights flat fields
            weight_1 = (reference_index_2-projection)/max(reference_index_2-reference_index_1, 1)
            weight_2 = 1 - weight_1                      
            # Open 'previous' flat
            reference_filename_1 = projection_prefix / ('refHST' + str(reference_index_1).zfill(4) + '.edf')
            #print(reference_filename_1)
            reference_image_1 = EdfFile(reference_filename_1).GetData(0) - dark_image
            # open 'next' flat
            reference_filename_2 = projection_prefix / ('refHST' + str(reference_index_2).zfill(4) + '.edf')
            #print(reference_filename_2)
            reference_image_2 = EdfFile(reference_filename_2).GetData(0) - dark_image
            # reference for current projection
            reference_image = weight_1*reference_image_1 + weight_2*reference_image_2
            # Open projection image                                  
            projection_image = (projection_image-dark_image) / reference_image

        if self.magnification and magnification:
            if self.magnification_factor[position] != self.max_magnification_factor: # Unless its the highest magnification
                projection_image = ndimage.zoom(projection_image, self.max_magnification_factor/self.magnification_factor[position])
            #TODO: Simple approach first. Magnification and center should perhaps be calculated in shifts?
        
        correction = np.zeros(3)

        projection_image = np.pad(projection_image, ((ypad, ypad), (xpad, xpad)), 'edge') #TODO: watch out with magnification

        if self.correct_motion and self.motion_corrected and aligned:
            #Apply motion
            #TODO: Is there any reason to not apply motion correction to reconstructed image like this?
            #TODO: Should be better to just motion correct the reference plane no?
            motion = self.get_motion(projection, self.reference_position)
            correction[0] = motion[0] * self.max_magnification_factor/self.magnification_factor[self.reference_position]
            correction[1] = motion[1] * self.max_magnification_factor/self.magnification_factor[self.reference_position]

        if self.correct_alignment and self.aligned and aligned:
            correction = correction + self.get_alignment(projection, position)
            
            # TODO: Why am I not using freq. variables correctly? investigate
            # TODO: Handle padding somehow/correcly
        if sum(correction*correction)**0.5 > 0.001:
            f = np.arange(0, self.nfx//2, 1)
            f = f * ((correction[0])*2*pi / self.nfx)
            f = np.exp(1j*f)
            #f = np.cos(f)-1j*np.sin(f)
            f = np.hstack((f, np.conj(f[::-1])))
            f = np.tile(f, (self.nfy, 1))
            
            g = np.arange(0, self.nfy//2, 1)
            g = g * ((correction[1]) * 2 * pi / self.nfy)
            g = np.exp(1j*g)
            #g = np.cos(g)-1j*np.sin(g)
            g = np.hstack((g, np.conj(g[::-1])))
            g = np.tile(g, (self.nfx, 1))
            g = g.transpose()
            t = f * g

            projection_image_FT = np.fft.fft2(projection_image) * t
        else:
            projection_image_FT = np.fft.fft2(projection_image)
        
        projection_image = np.real(np.fft.ifft2(projection_image_FT))
        projection_image = projection_image[self.ny//2:-self.ny//2, self.nx//2:-self.nx//2]
            # watch out how to keep track... corrected image is already in FD
            #fig = pyplot.imshow(np.real(np.fft.ifft2(projection_image)))
            #pyplot.colorbar()

        if Fourier:
            return projection_image_FT
        else:
            if pad:
                projection_image = np.pad(projection_image, ((ypad, ypad), (xpad, xpad)), 'edge')
            return projection_image
        
    def preprocess(self):
        """Runs all preprocessing on the dataset: axis positions, motion, 
        reference position, alignment"""
        self.calculate_axis_positions() # TODO: Verify order of shifts, motion, rotaxs
        self.calculate_motion()
        self.calculate_reference_position()
        self.align()
        self.fit_alignment()
        #self.display_alignment() #TODO: does not work for the moment
        self._write_parameter_file()
            
    def _write_parameter_file(self):
        """Writes current state to the associated parameter file"""
        parameters = dict(reference_position=self.reference_position,
                          nD=self.nD,
                          nbprojections=self.nbprojections,
                          position_number=self.position_number,
                          position=self.position.tolist(),
                          dss=self.distance_source_sample,
                          energy=self.energy,
                          detector_pixel_size=self.detector_pixel_size,
                          nx=self.nx,
                          ny=self.ny,
                          axis_position=self.axis_position.tolist(),
                          reference_interval=self.reference_interval,
                          scan_range=self.scan_range,
                          version=self.version,
                          alpha=self.alpha.tolist(),
                          delta_beta=self.delta_beta)
        
        with open(self.dataset_parameter_filename, 'w') as f:
            yaml.dump(parameters, f)
        pass
    
    def _read_parameter_file(self):
        """Reads the associated parameter file."""
        with open(self.dataset_parameter_filename, 'r') as f:
            parameters=yaml.safe_load(f)
        
        self.reference_position=parameters['reference_position']
        self.nD=parameters['nD']
        self.nbprojections=parameters['nbprojections']
        self.position_number=parameters['position_number']
        self.position=np.array(parameters['position'])
        self.distance_source_sample=parameters['dss']
        self.energy=parameters['energy']
        self.detector_pixel_size=parameters['detector_pixel_size']
        self.nx=parameters['nx']
        self.ny=parameters['ny']
        self.axis_position=np.array(parameters['axis_position'])
        self.reference_interval=parameters['reference_interval']
        self.scan_range=parameters['scan_range']
        self.version=parameters['version']
        self.alpha=np.array(parameters['alpha'])
        self.delta_beta=parameters['delta_beta']
        pass
    
    def populate(self):
        """Tries to find the dataset parameters in the accompanying *.info* and *xml* files.
        It is called when a parameter file ``pyphase_parameter.yml`` is not found.

        Returns
        -------
        self
           Returns *self* with the updated attributes.
        """

        # self.N  # projection angles

        # Checks number of directories named 'name_%d_'
        print(self.path)
        type(self.path)
        directories = os.listdir(self.path)
        name_regex = r"%s(\_)(\d)(\_)" % self.name #TODO: Too unreadable. Avoid cryptic code
        data_dirs = [dir for dir in directories if re.match(name_regex, dir)]

        self.alpha = np.array([-8, -10]) #TODO: default value... should this really go here?

        self.nD = len(data_dirs)  # number of distances
        self.position_number = np.arange(1, self.nD + 1).tolist()

        #TODO: Not nice, should be based on the generic paths
        fname_info = self.path / (self.name + '_1_') / (self.name + '_1_' + '.info')
        fname_xml = fname_info.with_suffix('.xml')

        print(fname_info)
        print(fname_xml)
        
        if fname_xml.exists():
            xml_file = minidom.parse(str(fname_xml))
            dists = []
            for dist in range(1, self.nD + 1):
                fname = self.path / (self.name + '_' + str(dist) + '_') / (self.name + '_' + str(dist) + '_' + '.xml')
                print(fname)
                xml_f = minidom.parse(str(fname))
                distance = float(xml_f.getElementsByTagName('distance')[0].firstChild.data)
                dists.append(distance)
            self.reference_position = 0
            self.nbprojections = int(xml_file.getElementsByTagName('tomo_N')[0].firstChild.data)
            self.position = np.array(dists) / 1000
            self.distance_source_sample = 144
            self.energy = float( xml_file.getElementsByTagName('energy')[0].firstChild.data)
            self.detector_pixel_size = float(xml_file.getElementsByTagName('pixelSize')[0].firstChild.data)
            self.nx = int(xml_file.getElementsByTagName('DIM_1')[0].firstChild.data)
            self.ny = int(xml_file.getElementsByTagName('DIM_2')[0].firstChild.data)
            self.reference_interval = int(xml_file.getElementsByTagName('ref_On')[0].firstChild.data)
            self.scan_range = int(xml_file.getElementsByTagName('scanRange')[0].firstChild.data)
        elif fname_info.exists():
            # Adding file information into a dictionary
            parameter_dict = {}
            with open(fname_info) as file:
                for line in file:
                    listedline = line.strip().split('=')
                    parameter_dict[listedline[0]] = listedline[1].strip()
            #self.reference_distance = 1  # stationary plane in registration #D# TODO: Where to find this?
            self.nbprojections = int(parameter_dict['TOMO_N'])
            # Reading distances from the .info files
            dists = []
            for dist in range(1, self.nD + 1):
                fname = elf.path / (self.name + '_' + str(dist) + '_') / (self.name + '_' + str(dist) + '_' + '.info')
                with open(fname) as file:
                    for i, line in enumerate(file):
                        if i == 1:  # 'distance' is in the second line
                            listedline = line.strip().split('=')
                            distance = float(listedline[1].strip()) * 0.01
                            dists.append(distance)
            self.position = np.array(dists) / 1000
            self.distance_source_sample = 145  # Distance Source - Sample TODO: Where to find this?
            self.energy = float(parameter_dict['Energy'])  # keV
            self.detector_pixel_size = float(parameter_dict['Optic used'])
            self.nx = int(parameter_dict['Dim_1'])
            self.ny = int(parameter_dict['Dim_2'])
            self.reference_interval = int(parameter_dict['REF_ON'])
            self.scan_range = int(parameter_dict['ScanRange'])

        else: raise Exception('Supporting files '+ self.name +'_1_.info and ' + self.name +'_1_.xml not found.')

    def GetSinogram(self, distance=1):
        """Return projections as sinogram

        Parameters
        ----------
        distance : int
           Number of the distance of the projection.
           
        Returns
        -------
        sinogram
           Returns projections stacked as sinogram
        """
        #TODO: RIT/YAGNI. Add functionality as needed
        #TODO: Probably has to be able to split over lines for serialisation/parallellisation
        
        #Allocate a 3D numpy array of correct dimensions
        sinogram = np.empty([self.nx, self.ny, self.nbprojections])
        #Load projections one by one
        for projection in range(self.nbprojections):
            sinogram[:,:,projection]=self.get_projection(projection, self.position_number[distance])
        #return
        return sinogram

    def get_image(self, *, projection, image_type='phase', Fourier=False, pad=False):
        """Get reconstructed images (phase by default)

        Parameters
        ----------
        projection : int
           Number of the projection.
        image_type : str, optional
            Can be phase or attenuation.
        Fourier : bool
            Return Fourier transform of image.
        pad : bool
            Return padded image.

        Returns
        -------
        image : ndarray
        """
        if image_type.lower() == 'attenuation':
            fname = self.attenuation_prefix.as_posix()+'_'+str(projection).zfill(4)+'.edf'
#        elif 'propagated' in argv:        
#            fname = self.propagated_prefix.as_posix()+'_'+str(distance)+'_'+str(projection).zfill(4)+'.edf'
#        elif 'difference' in argv:
#            fname = self.difference_prefix.as_posix()+'_'+str(distance)+'_'+str(projection).zfill(4)+'.edf'
#        elif 'prior' in argv:
#            fname = self.prior_forward_prefix.as_posix()+('_'+str(projection).zfill(4)+'.edf')
        elif image_type.lower() == 'phase':
            fname = self.phase_prefix.as_posix()+'_'+str(projection).zfill(4)+'.edf' #TODO: separate file type...
        else:
            raise TypeError("proj_type must be one of the following: 'phase', 'attenuation'" )

        image = EdfFile(fname).GetData(0)
        xpad = (self.nfx-self.nx)//2
        ypad = (self.nfy-self.ny)//2
        
        if Fourier:
            #TODO: Refactor into an FT function in utilities (f ex)
            return np.fft.fft2(np.pad(image, ((ypad, ypad), (xpad, xpad)), 'edge'))
        else:    
            if pad:
                image = np.pad(image, ((ypad, ypad), (xpad, xpad)), 'edge')
            return image
    
    def calculate_axis_position(self, distance):
        """Calculates the axis position at one position"""
        # During/after scan
        # 180/360

        correct_alignment = self.correct_alignment 
        self.correct_alignment = 0
        
        correct_motion = self.correct_motion
        self.correct_motion = 0
        
        if self.scan_range == 360:
            projection_0 = self.nbprojections + 4 
            projection_180 = self.nbprojections + 2
        else:
            projection_0 = self.nbprojections + 2
            projection_180 = self.nbprojections
        
        im0 = self.get_projection(projection=projection_0, position=distance)  
        im180 = np.fliplr(self.get_projection(projection=projection_180, position=distance))
        
        field, image_transormed, transform_parameters = self.registrator.register(im0, im180)
        offset = field[0][0][0]/2
        print("Rotation axis offset: {}".format(offset))
        self.axis_position[distance] = self.nx/2 + offset
        print("Axis position: {}".format(self.axis_position[distance]))   
        
        self.correct_alignment = correct_alignment
        self.correct_motion = correct_motion
        
    def calculate_axis_positions(self):
        """Calculates all axis positions."""
        self.axis_position = np.zeros(self.nD)
        for distance in range(self.nD):
            self.calculate_axis_position(distance)
        self._write_parameter_file()

    def calculate_motion(self):
        """Estimates motion in the scans based on supplementary images."""
        # define reference positions
        self.motion_correlations = [ [] for i in range(self.nD) ]
        self.motion_corrected = 0

        correct_alignment = self.correct_alignment 
        self.correct_alignment = 0
        
        reference_positions = range(0,self.scan_range, 90)
        reference_position_indices = range(self.nbprojections+int(self.scan_range/90), self.nbprojections, -1) 
               
        image = []        
        reference = []
        
        for Distance in range(self.nD):
            print("Calculating motion position {}".format(str(Distance)))
            for index, position in enumerate(reference_positions):
                # Get reference projecions and scan projections at reference positions
                image_index = self.nbprojections * (position/self.scan_range)
                reference = self.get_projection(projection=reference_position_indices[index], position=Distance) 
                image = self.get_projection(projection=floor(image_index), position=Distance)           
                if image_index != floor(image_index):
                    print("Interpolating projection for motion estimation")
                    image2 = self.get_projection(floor(image_index)+1, Distance)
                    image = (floor(image_index)+1-image_index)*image + (image_index-floor(image_index))*image2           
                
                # calculate correlations
                field, aligned_image, transform_parameters = self.registrator.register(image, reference)
                # Will create a list of lists with [distance][angle (corresponding to ref positions, 0=0, 1=90, ...)][dim (0=x,1=y,...)]
                self.motion_correlations[Distance].append([field[0][0][0], field[1][0][0]])

        # save motions
        with open(self.motion_filename, 'wb') as f:
            pickle.dump(self.motion_correlations, f, pickle.HIGHEST_PROTOCOL)
        
        self.motion_corrected=1        
        self.correct_alignment = correct_alignment
        pass
            
    def get_motion(self, projection, distance):
        """Returns estimation motion at one projection and position."""
        #TODO: disgusting, needs cleaning
        # fit motion with sinusoids (?)
        #TODO: Should be called from get_projection no???
        if self.correct_motion and self.motion_corrected:
            reference_positions = np.arange(0, self.scan_range, 90)
            projections_all = np.arange(self.nbprojections+len(reference_positions)+1)
            polynomial_order_h = min(floor(len(reference_positions)/2), 3) 
            polynomial_order_v = min(len(reference_positions), 3)
            arg = np.array(reference_positions)/self.scan_range - 1 # TODO: to be extended in case of more images a la id16a?
            arg_all = projections_all/self.nbprojections - 1
            corr=np.array(self.motion_correlations[distance-1]) 
            # Vertical fit
            M = np.tile(arg, [polynomial_order_v, 1])
            for k in range(1, polynomial_order_v):
                M[k,:] = M[k,:]*M[k-1,:]
            M = -np.flipud(M) #TODO: could certainly be more elegant...
                
            coefficients = np.dot(np.linalg.pinv(M).T, corr[:,1])
            # subtract mean from vertical motion
            ref_v = np.polyval(np.append(coefficients,0), arg_all)
            ref_v_mean = np.mean(ref_v);
            ref_v -= ref_v_mean;
            ref_v_meas = -corr[:,1]-ref_v_mean
            # Horizontal fit
            
            M = np.tile(arg, [polynomial_order_h, 1])
            for k in range(1, polynomial_order_h):
                M[k,:] = M[k,:]*M[k-1,:]
            M = -np.flipud(M) #TODO: could certainly be more elegant...
            
            M = np.append(M*np.tile(np.cos((np.pi/180)*reference_positions),[polynomial_order_h,1]), 
                          M*np.tile(np.sin((np.pi/180)*reference_positions),[polynomial_order_h,1]),axis=0)

            coefficients = np.dot(np.linalg.pinv(M).T, corr[:,0])
            sx_cal = np.polyval(np.append(coefficients[0:polynomial_order_h], 0), arg_all);
            sy_cal = np.polyval(np.append(coefficients[polynomial_order_h:], 0), arg_all);
            arg = projections_all*(self.scan_range/180*pi/self.nbprojections);
            reference_position_indices = self.nbprojections + np.arange(int(self.scan_range/90), 0, -1)        
            arg[reference_position_indices] = np.array(reference_positions)*(pi/180)
            M = np.vstack((np.cos(arg), np.sin(arg)));
            ref_h = sx_cal*M[0,:] + sy_cal*M[1,:]
            #Fit with sinusoid
            coefficients = np.dot(np.linalg.pinv(M).T,ref_h)
            sx_cal -= coefficients[0]; 
            sy_cal -= coefficients[1];
            ref_h -= np.dot(coefficients,M);
            index_rots = np.round(np.array(reference_positions)*(self.nbprojections/self.scan_range))
            ref_h_meas=-corr[:,0]
            ref_h_meas -= np.dot(coefficients, M[:,index_rots.astype(int)])
            #pyplot.figure()
            #pyplot.plot(projections_all, ref_v, 'b', index_rots, ref_v_meas, 'bx', projections_all, ref_h, 'r', index_rots, ref_h_meas, 'rx')
            return ref_h[projection], ref_v[projection]
        else:
            return 0, 0
        
    def calculate_reference_position(self):
        """Calculates and sets reference position based on motion""" 
        #TODO: add check if motion is calculated (self.motion_correlations)
        #TODO: manipulation of the motion list is quite roundabout
        correlation_amplitudes = []
        for index in range(self.nD):
            tmp = [item for sublist in self.motion_correlations[index] for item in sublist] #flatten list
            tmp = [x**2 for x in tmp] #square elements
            correlation_amplitudes.append(sqrt(sum(tmp)))
        self.reference_position = correlation_amplitudes.index(min(correlation_amplitudes))
        print('Reference position: {}'.format(self.reference_position))
        self._write_parameter_file()
#    def IsAligned(self):
#        return self.aligned

    # def ClearShifts(self, D):
    #     """clear shifts. Either 1 distance or D='all'"""
    #     pass

    def fit_alignment(self):
        """Fits measured alignment parameters with polynomials"""
        nb_shifts=len(self.shifts)
        
        projection=[0 for x in range(nb_shifts)]
        x=[[0 for x in range(nb_shifts)] for y in range(self.nD)] 
        y=[[0 for x in range(nb_shifts)] for y in range(self.nD)]
        theta=[[0 for x in range(nb_shifts)] for y in range(self.nD)]
        
        for p_i in range(nb_shifts):
            projection[p_i]=self.shifts[p_i][0]
        
        for d_i in range(self.nD):
            for p_i in range(nb_shifts):
                x[d_i][p_i] = self.shifts[p_i][1][d_i][0]
                y[d_i][p_i] = self.shifts[p_i][1][d_i][1]
                #theta[dndx][projndx]=self.shifts[projndx][1][dndx][2]
        
        for d_i in range(self.nD):
            self.shift_polynomial_x[d_i]=np.polyfit(projection, x[d_i], self.shift_polynomial_order[0])
            self.shift_polynomial_y[d_i]=np.polyfit(projection, y[d_i], self.shift_polynomial_order[1])
        pass

    def display_alignment(self, projection=0):
        """Displays alignment and fit."""
        # TODO: quite a lot of duplicate code with FitShifts
        nb_shifts=len(self.shifts)
        projection_range=[0 for x in range(nb_shifts)]
        x=[[0 for x in range(nb_shifts)] for y in range(self.nD)] 
        y=[[0 for x in range(nb_shifts)] for y in range(self.nD)] 
        theta=[[0 for x in range(nb_shifts)] for y in range(self.nD)]
        X=np.linspace(0, self.nbprojections-1, self.nbprojections)

        for d_i in range(self.nD):
            for p_i in range(nb_shifts):
                x[d_i][p_i] = self.shifts[p_i][1][d_i][0]
                y[d_i][p_i] = self.shifts[p_i][1][d_i][1]
                #theta[dndx][projndx]=self.shifts[projndx][1][dndx][2]

        for p_i in range(nb_shifts):
            projection_range[p_i]=self.shifts[p_i][0]

        pyplot.figure()
        for index, distance in enumerate(self.position_number_no_reference_position):
            pyplot.subplot(len(self.position_number_no_reference_position), 2, 2*(index)+1)
            pyplot.plot(projection_range, x[distance-1], 'rx', X, np.polyval(self.shift_polynomial_x[distance-1], X), 'b')
            pyplot.title('Distance '+str(distance))
        
        for index, distance in enumerate(self.position_number_no_reference_position):
            #yi=np.polyfit(projection_range, y[distance], self.shift_polynomial_order[1])
            pyplot.subplot(len(self.position_number_no_reference_position), 2, 2*(index)+2)
            pyplot.plot(projection_range, y[distance-1], 'rx', X, np.polyval(self.shift_polynomial_y[distance-1], X), 'b')
            pyplot.title('Distance '+str(distance))
                
        pyplot.subplots_adjust(hspace=.6, left=.07, bottom=0.05, right=.98, top=.95)
        
        projection_images = np.zeros([self.ny, self.nx, self.nD])
        for distance in range(1, self.nD+1):
            projection_images[:,:,distance-1] = self.get_projection(projection=projection, position=distance)
        viewer.display_stack(projection_images)

    # def SetShift(self, projection, distance, x, y, theta=0):
        
    #     pass
    
    def get_alignment(self, projection, distance):
        """Returns the transform parameters for alignment at one projection
        and position"""
        #TODO: stub
#        if distance == self.reference_distance or self.correct_alignment == 0:
#            shift_x = 0
#            shift_y = 0
#            shift_theta = 0
#        else:
        shift_x = np.polyval(self.shift_polynomial_x[distance-1], projection)
        shift_y = np.polyval(self.shift_polynomial_y[distance-1], projection)
        shift_theta = 0
        
        return [shift_x, shift_y, shift_theta] 

    def align(self, interval=100):
        """
        Align(self, RA, interval=100)
        Align complete data set.

        TODO: Should probably check for the last projection also so that it is
        always included
        """ 
        self.position_number_no_reference_position = self.position_number.copy()
        self.position_number_no_reference_position.remove(self.reference_position)

        correct_alignment = self.correct_alignment
        self.correct_alignment = 0
        correct_motion = self.correct_motion
        self.correct_motion=0
        
        self.shifts = SortedList()

        for projection in range(0, self.nbprojections, interval):
            text = "\r" + "Aligning projection {}".format(projection)
            sys.stdout.write(text)
            sys.stdout.flush()
            self.align_projection(projection)
                    
        # Pickle the 'data' dictionary using the highest protocol available.
        with open(self.shift_filename, 'wb') as f:
            pickle.dump(self.shifts, f, pickle.HIGHEST_PROTOCOL)
                
        self.aligned=1
        self.correct_alignment = correct_alignment
        self.correct_motion = correct_motion
                
    def align_projection(self, projection):
        """
        projection = number of projection to register
        """
        correct_alignment = self.correct_alignment 
        self.correct_alignment = 0
                 
        im=[]

        # read images
        for position in range(self.nD):
            im.append(self.get_projection(projection=projection, position=position))

        # register (should be done with RA)
        tmp=[]
        
        for ndx in range(self.nD):
            if ndx != self.reference_position:
                field, aligned_image, transform_parameters = registrator.register(im[ndx], im[self.reference_position])
                tmp.append([field[0][0][0], field[1][0][0]])
            else: 
                tmp.append([0, 0])            
                   
        # store shifts
        self.shifts.add([projection, tmp])
        
        with open(self.shift_filename, 'wb') as f:
            pickle.dump(self.shifts, f, pickle.HIGHEST_PROTOCOL)      
        
        self.correct_alignment = correct_alignment

#TODO: Difference just shouldn't be here..
    # def DifferenceProjection(self, projection, distance):
    #     # Calculate difference measured/calculated intensity
    #     # TODO: Should eventually be an option if it's tomo or not?
    #     # TODO: Should this really go here?

    #     # Open measured
    #     measured_projection = self.get_projection(projection, distance)
        
    #     # Get calculated
    #     calculated_projection = self.get_image(projection, 'propagated', distance)        
        
    #     # Calculate difference (correct order?)
    #     difference_image = measured_projection - calculated_projection
        
    #     # Save difference
    #     #TODO: IMPORTANT: USE DATASET TO READ/WRITE IMAGES!!!
    #     fname = self.path + '/' + self.name + '_1_/' + self.name + '_1_0000.edf' # TODO:Fulhack deluxe to get edfheader. Whatever. Stub/spike. Fix later
    #     imEDF = EdfFile(fname)
    #     EDFHeader = imEDF.GetHeader(0)
    #     time.sleep(1) #TODO: This is ridiculous, how can this be necessary? What is causing this writing problem?

    #     # TODO: the difference files should be handled in constructor as well I suppose
    #     fname = self.difference_prefix + '_' + str(distance) + '_' + str(projection).zfill(4) + '.edf'
    #     EDF = EdfFile(fname)
    #     # TODO: verify what should go into the header...
    #     EDF.WriteImage(EDFHeader, difference_image, 0, "Float")
        
    # def DifferenceProjections(self, projections, distance):
    #     if projections[0] > projections[1]:
    #         projections[1]=projections[0]
        
    #     for projection in range(projections[0], projections[1]+1):
    #         self.DifferenceProjection(projection, distance)
        
    # def Difference(self):
    #     #TODO: Needs to be parallelised... how? Come up with something...
    #     parallelizer = Parallelizer.OAR()
    #     for distance in range(1, self.nD+1):
    #         parallelizer.Launch(self, 'difference', distance=distance)

    # def UpdatePhase(self):
    #     # add update to current solution
    #     current_solution = np.memmap(self.phase_vol_filename, dtype=np.float32)
    #     update = np.memmap(self.update_vol_filename, dtype=np.float32)
    #     current_solution += update
    #     del current_solution, update

    # def UpdateAttenuation(self):
    #     current_solution = np.memmap(self.retrieved_attenuation_vol_filename, dtype=np.float32)
    #     update = np.memmap(self.attenuation_update_vol_filename, dtype=np.float32)
    #     current_solution -= update
    #     del current_solution, update

    # def write_image(self, data, proj_type, projection_number, *args):
    def write_image(self, *, image, projection, projection_type='phase'):
        # TODO: raise error where needed! number of args, args type...
        ''' Saves images into file.

        Parameters
        ----------
        data : ndarray
           An ndarray of the image to save.
        projection_number : int
           The number of the projection to be saved.
        proj_type : str
           A string containing the prefix for the name of the file to be saved.
        args: int
           Number of the distance of the projection.
        Returns
        -------
        None
           Saves the images into the file ``[prefix]_[projection_number].edf`` or
           ``[prefix]_[distance]_[projection_number].edf`` into the ``name_`` directory.
        '''
        if projection_type.lower() == 'phase':
            pref = self.phase_prefix          
        elif projection_type.lower() == 'attenuation':
            pref =  self.attenuation_prefix
        # elif projection_type.lower() == 'phase update':
        #     pref = self.update_prefix
        # elif projection_type.lower() == 'attenuation update':
        #     pref = self.attenuation_update_prefix
        # elif projection_type.lower() == 'intensity':
        #     pref = self.propagated_prefix
        else:
            raise TypeError("proj_type must be one of the following: 'phase', 'attenuation'" )
        
        if not pref.parent.exists():
            pref.parent.mkdir(parents=True)
        
#        if len(args) != 0:
#            distance = str(args[0]) + '_'
#        else:
#            distance = ''
#            
        fname = self.path / (self.name + '_1_/' + self.name + '_1_0000.edf')  # TODO:Fulhack deluxe
        
        imEDF = EdfFile(fname)
        EDFHeader = imEDF.GetHeader(0)
        #fname = pref.as_posix() + '_' + position + str(projection).zfill(4) + '.edf'
        fname = pref.as_posix() + '_' + str(projection).zfill(4) + '.edf'
        EDF = EdfFile(fname)
        EDF.WriteImage(EDFHeader, image, 0, "Float")

    #D# TODO: These methods are to be used in propagator.py (for later)
    #TODO: Seems like a mistake, since get_image already exists... I guess both could exist? Ev. Refactor
    # def get_attenuation(self, projection_number):
    #     ''' Reads attenuation information from file and retruns it as an ndarray.

    #     Parameters
    #     ----------
    #     projection_number : int
    #        The number of the projection to get.

    #     Returns
    #     -------
    #     ndarray
    #        Array containing attenuation information.
    #     '''

    #     projection = projection_number
    #     fname = self.retrieved_attenuation_forward_prefix + '_' + str(projection).zfill(4) + '.edf'
    #     imEDF = EdfFile(fname)
    #     attenuation = imEDF.GetData(0)
    #     return attenuation

    # def get_phase(self, projection_number):
    #     ''' Reads phase information from file and retruns it as an ndarray.

    #     Parameters
    #     ----------
    #     projection_number : int
    #        The number of the projection to get.

    #     Returns
    #     -------
    #     ndarray
    #         Array containing phase information.
    #     '''
    #     projection = projection_number
    #     fname = self.phase_forward_prefix + '_' + str(projection).zfill(4) + '.edf'
    #     imEDF = EdfFile(fname)
    #     phase = imEDF.GetData(0)
    #     return phase
#class ESRF(Dataset):
    #maybe necessary at some point
#    pass
        
#class ESRFID19(Dataset):
    
#    pass

#class ESRFID16A(Dataset):
#    def __init__(self):
#        super().__init__()
        
        
#    def populate(self):
        
#    pass
