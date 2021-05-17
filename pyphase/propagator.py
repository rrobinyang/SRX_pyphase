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

import numpy as np
from vendor.EdfFile import EdfFile #TODO: This should not be necessary here!
import scipy.ndimage as ndimage
from math import *
import pyphase.parallelizer as Parallelizer
import pyphase.utilities as Utilities

class Propagator:
    def __init__(self, *, dataset=None, shape=None, energy=None, pixel_size=None, pad=2, oversampling=4):
        """Takes either a dataset object or necessary parameters as keyword arguments
    
        Parameters
        ----------
        dataset : pyphase.Dataset, optional
            Dataset object for parameters
        shape : tuple of ints, optional
            Size of images (ny, nx) for creation of frequency variables etc
        pixel_size : float , optional
            In m
        distance : list of floats, optional
            Effective propagation distances in m
        energy : float, optional
            Effective energy in keV
        alpha : tuple of floats, optional
            Regularisation parameter
        pad : int, optional
            Padding factor for propagator
        oversampling : int, optional
            Oversampling factor of projections
        """
        self.padding=pad
        self.oversampling=oversampling
        self.lengthscale=10e-6 # where should this be...
        if dataset:
            self.pixel_size=dataset.pixel_size*1e-6
            self.nx=dataset.nx
            self.ny=dataset.ny
            self.Lambda=dataset.Lambda
        else:
            self.nx=shape[1]
            self.ny=shape[0]
            self.Lambda = 12.4e-10 / energy
            self.pixel_size=pixel_size
        
        self.nfx = self.nx*pad 
        self.nfy = self.ny*pad
        
        self.sample_frequency=self.lengthscale/self.pixel_size
        
        # creation of frequency variable
        # According to numpy fft convention; 
        # a[0] should contain the zero frequency term,
        # a[1:n//2] should contain the positive-frequency terms,
        # a[n//2 + 1:] should contain the negative-frequency terms, in increasing 
        # order starting from the most negative frequency.
    
        # self.x=0
        # self.x=np.append(self.x, np.linspace(self.sample_frequency/self.nfx, self.sample_frequency/2, self.nfx//2))
        # self.x=np.append(self.x, np.linspace(-self.sample_frequency/2 + self.sample_frequency/self.nfx, -self.sample_frequency/self.nfx, self.nfx//2-1 + (np.ceil(self.nfx/2)-self.nfx//2)))
        
        # self.y=0
        # self.y=np.append(self.y, np.linspace(self.sample_frequency/self.nfy, self.sample_frequency/2, self.nfy//2))
        # self.y=np.append(self.y, np.linspace(-self.sample_frequency/2 + self.sample_frequency/self.nfy, -self.sample_frequency/self.nfy, self.nfy//2-1 + (np.ceil(self.nfy/2)-self.nfy//2)))
        
        # self.fx, self.fy = np.meshgrid(self.x, self.y)
            
#    def Propagate(self, dataset, distance=''):
#        if distance:
#            parallelizer = Parallelizer.OAR()
#            parallelizer.Launch(dataset, 'propagate', algorithm=self.__class__.__name__, distance=distance)
#        else:            
#            for distance in range(1, dataset.nD+1):
#                parallelizer = Parallelizer.OAR()
#                parallelizer.Launch(dataset, 'propagate', algorithm=self.__class__.__name__, distance=distance)
                 
    # def PropagateProjections(self, dataset, projections, distance):
    #     if projections[0] > projections[1]:
    #         projections[1]=projections[0]
        
    #     for projection in range(projections[0], projections[1]):
    #         self.PropagateProjection(dataset, projection, distance)
        

class Fresnel(Propagator):
    """Propagator using Fresnel transform"""
    def PropagateProjection(self, *, dataset=None, position_number=None, projection=None, phase=None, attenuation=None, position=None, oversampled=False):
        """
        Propagate one projection.
        
        Arguments
        ---------
        dataset : pyphase.Dataset, optional
            Datset with projection data.
        position_number : int, optional
            Which position to propagate to
        projection : int, optional
            Which projection to propagate
        phase : ndarray, optional
            Phase of wave to propagate
        attenuation : ndarray, optional
            Amplitude of wave to propagate
        position : float
            Effective propagation distance
        oversampled : bool
            True if imput images are already oversampled
        
        """
        #TOOD: Should check if phase and attenuation are calculated
        #TODO: Allow propagation of a single projection and a range of projections
        #TODO: Split in image/projection like phaseretrieval

        length_scale=10e-6
                
        ps=self.pixel_size/self.oversampling
        fs=length_scale/ps
        nx = self.nx*self.oversampling
        ny = self.ny*self.oversampling
        x=np.linspace((-ps*nx/2), ps*(nx/2-1), nx)
        y=np.linspace((-ps*ny/2), ps*(ny/2-1), ny)        
        xx, yy = np.meshgrid(x, y)
        
        f=np.linspace(-fs/2, fs/2-fs/(nx*self.padding), nx*self.padding)
        g=np.linspace(-fs/2, fs/2-fs/(ny*self.padding), ny*self.padding)
        ff, gg = np.meshgrid(f, g)        
        
        if dataset:
            position = dataset.position[position_number]
            phase = dataset.GetImage(projection=projection, image_type='generated_phase') #TODO: generated or not should be choise. Not yet implemented
            attenuation = dataset.GetImage(projection=projection, image_type='generated_attenuation')
 
        if  self.oversampling and not oversampled:
            phase = ndimage.zoom(phase, self.oversampling)
            attenuation = ndimage.zoom(attenuation, self.oversampling)
        
        # TODO: Creation of propagator should probably be in constructor? So that one can get it out
       
        wave = np.exp(-attenuation+1j*phase) #TODO: Decide form of projection. If attenuation in mu, no square?
        # TODO: THIS PADDING IS DANGEROUS NO?! VERIFY!
        Utilities.resize(wave, [ny*self.padding, nx*self.padding])
        wave = np.fft.fft2(wave)
        P=np.fft.ifftshift(np.exp(-1j*pi*self.Lambda*position*(ff**2+gg**2)/(length_scale**2)))
        Id=np.fft.ifft2(wave*P)
        Id=np.abs(Id)**2
        Id = Utilities.resize(Id, [ny, nx])
        Id = ndimage.zoom(Id,1/self.oversampling)
        if dataset:
            dataset.WriteImage(Id, 'intensity', projection, position_number)
        
        return Id            
        
class CTF(Propagator):
    """Propagates using the CTF. Legacy code to be aligned with Fresnel.""" 
    def PropagateProjection(self, dataset, projection, distance):    
        length_scale=10e-6
        oversampling=1
        padding=2
        ps=(dataset.pixel_size/oversampling)*1e-6
        fs=1/ps
        nx = dataset.nx*oversampling
        ny = dataset.ny*oversampling
        x=np.linspace((-ps*nx/2), ps*(nx/2-1), nx)
        y=np.linspace((-ps*ny/2), ps*(ny/2-1), ny)        
        xx, yy = np.meshgrid(x, y)
        
        f=np.linspace(-fs/2, fs/2-fs/(nx*padding), nx*padding)
        g=np.linspace(-fs/2, fs/2-fs/(ny*padding), ny*padding)
        ff, gg = np.meshgrid(f, g)        

        #TODO: should have getters and setters for all the images. How? one per type? one for all w arg?
        fname = dataset.phase_forward_prefix+'_'+str(projection).zfill(4)+'.edf'
        imEDF = EdfFile(fname)
        phase = imEDF.GetData(0)
        EDFHeader = imEDF.GetHeader(0)
        phase = ndimage.zoom(phase, oversampling)
        phase = np.pad(phase, ((ny//2, ny//2), (nx//2, nx//2)), 'edge')

        fname = dataset.retrieved_attenuation_forward_prefix+'_'+str(projection).zfill(4)+'.edf'
        imEDF = EdfFile(fname)
        attenuation = imEDF.GetData(0)
        attenuation = ndimage.zoom(attenuation, oversampling)
        attenuation = np.pad(attenuation, ((ny//2, ny//2), (nx//2, nx//2)), 'edge')

        #P = [0 for q in range(DS.nD)]holosim_PP_prop_4_0000.edf
        #Id = [0 for q in range(DS.nD)]
        
        # TODO: Creation of propagator should probably be in constructor? So that one can get it out
       
        FN = dataset.Fresnel_number[distance-1]
        coschirp = np.cos((pi*FN) * (self.fx**2) + (pi*FN) * (self.fy**2))
        sinchirp = np.sin((pi*FN) * (self.fx**2) + (pi*FN) * (self.fy**2))

       
#       for n in DS.distance_number:
       
        FId = 2 * coschirp * np.fft.fft2(attenuation) - 2 * sinchirp * np.fft.fft2(phase)

        Id = 1 + np.real(np.fft.ifft2(FId))
        Id = Id[ny//2:-ny//2, nx//2:-nx//2]
        Id = ndimage.zoom(Id,1/oversampling)
        fname=dataset.propagated_prefix+'_'+str(distance)+'_'+str(projection).zfill(4)+'.edf'
        EDF = EdfFile(fname)
        # TODO: verify what should go into the header...
        EDF.WriteImage(EDFHeader, Id, 0, "Float")
    pass 
