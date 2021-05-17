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
from math import *
import pyphase.parallelizer as Parallelizer
import pyphase.propagator as Propagator
import pyphase.tomography as Tomography
#from vendor.EdfFile import EdfFile #TODO: This should not be necessary here!
import matplotlib.pyplot as pyplot
import scipy.ndimage
from pyphase.config import *
from matplotlib.pyplot import pause
from scipy import interpolate
import pickle #TODO: Is there a better way to handle imports? Centralised?
from scipy import ndimage
import pyphase.dataset as Dataset

class PhaseRetrievalAlgorithm2D:
    """Base class for 2D phase retrieval algorithms.

    Parameters
    ----------
    dataset : pyphase.Dataset, optional 
        A Dataset type object.
    shape : tuple of ints, optional
        Size of images (ny, nx) for creation of frequency variables etc.
    pixel_size : float, optional
        In m.
    distance : list of floats, optional
        Effective propagation distances in m.
    energy : float, optional
        Effective energy in keV.
    alpha : tuple of floats, optional
        Regularisation parameters. First entry for LF, second for HF. 
        Typically [1e-8, 1e-10].
    pad : int
        Padding factor (default 2).

    Attributes
    ----------
    nx : int
        Number of pixels in horizontal direction.
    ny : int
        Number of pixels in horizontal direction.
    pixel_size : tuple of floats
        Pixel size [x, y] in µm.
    ND : int
        Number of positions.
    energy : float
        Energy in keV.
    alpha : tuple of floats
        First entry for LF, second for HF. Tyically [1e-8, 1e-10].
    distance : numpy array
        Effective propagation distances in m.
    padding : int
        Padding factor.
    sample_frequency : float
        Reciprocal of pixel size in lengthscale.
    nfx : int
        Number of samples in Fourier domain (horizontal).
    nfy : int
        Number of samples in Fourier domain (vertical).
    fx : numpy array
        Frequency variable (horizontal), calculated by frequency_variable.
    fy : numpy array
        Freuency variable (vertical).
    alpha_cutoff : float
        Cutoff frequency for regularization parameter in normalized frequency.
    alpha_slope : float
        Slope in regularization parameter.
    
    Notes
    -----
    Takes either a Dataset type object with keyword dataset (which contains
    all necessary parameters), or parameters as above.
    
    """
    
    def __init__(self, dataset=None, **kwargs):
        self.lengthscale=10e-6 # TODO where should this be... Makes formulae dimensionless
        #print(kwargs)
        # Check whether usage is with dataset or external images 
        #TODO: Which parameters in constructor?
        if (isinstance(dataset, Dataset.Dataset) or isinstance(dataset, Dataset.ESRF)): #TODO: needs integration 
            #dataset = args[0]
            self.nx=dataset.nx # TODO: Probably make into an array [x, y]
            self.ny=dataset.ny
            self.pixel_size = dataset.pixel_size*1e-6
            self.ND = len(dataset.position)
            self.energy = dataset.energy
            self.distance = np.array(dataset.effective_distance) #TODO: wallah I've mixed up position, distance, effective distance. Sort out
        else:
            self.nx = kwargs["shape"][1]
            self.ny = kwargs["shape"][0]
            self.pixel_size = kwargs["pixel_size"]
            self.distance = np.array(kwargs["distance"])
            self.energy = kwargs["energy"]
            self.ND = len(self.distance)

        if 'pad' in kwargs:
            self.padding = kwargs['pad']

        self.alpha=[-8, -10]
 
        if (type(self.pixel_size) == float) or (type(self.pixel_size) == np.float64) : # If x,y pixelsizes are not given, assign them
            self.pixel_size = np.array([self.pixel_size, self.pixel_size])
        elif type(self.pixel_size) == list:
            self.pixel_size = np.array(self.pixel_size)
        self.sample_frequency = self.lengthscale/self.pixel_size #TODO: Should be attribute
        self.fx, self.fy = self.frequency_variable(self.nfx, self.nfy, self.sample_frequency)
        self._compute_factors()
        self.alpha_cutoff = .5
        self.alpha_cutoff_frequency=self.alpha_cutoff*self.sample_frequency # TODO: should be a property (dynamically calculated from alpha_cutoff)
        self.alpha_slope = .1e3

    def _algorithm(self, image, positions=None):
        """'Pure virtual' method containing purely the algorithm code part. Should be defined by each subclass. """

    def _compute_factors(self):
        """
        Computes factors used in phase retrieval. 
        
        Motivation is to save time and memory when reconstructing a series of 
        projections on a processor (usually in parallel with others). Default 
        are the sin and cos chirps in CTF. CTF and Mixed approaches extend 
        this method, while TIE methods override it.
        """
        
        self.coschirp = np.zeros((self.ND, self.nfy, self.nfx))
        self.sinchirp = np.zeros_like(self.coschirp)
        for distance in range(self.ND):
            self.coschirp[distance] = np.cos((pi*self.Fresnel_number[distance]) * (self.fx**2) + (pi*self.Fresnel_number[distance]) * (self.fy**2))
            self.sinchirp[distance] = np.sin((pi*self.Fresnel_number[distance]) * (self.fx**2) + (pi*self.Fresnel_number[distance]) * (self.fy**2))

    def __getstate__(self):
        """
        Used in parallel computing to override writing of voluminous variables
        when serializing using pickle by parallelizer. They are instead 
        re-calculated by each process (see __setstate__()).
        """
        state = self.__dict__.copy()
        del state['fx'], state['fy'], state['sinchirp'], state['coschirp']
        return state

    def __setstate__(self, state):
        """
        Used in parallel computing to re-calculate voluminous variables when
        serializing using pickle by parallelizer instead of saving them to disk
        (see __setstate__()).
        """
        self.__dict__.update(state)
        self._compute_factors()
             
    @property
    def nfx(self):
        return self.padding*self.nx
    
    @property
    def nfy(self):
            return self.padding*self.ny
                
    @property
    def Lambda(self):
        """Wavelength based on energy (float)"""
        return 12.4e-10 / self.energy

    @property
    def Fresnel_number(self):
        """Fresnel number at each position, calculated from energy and distance (float)"""
        return self.Lambda * self.distance / (self.lengthscale**2)

    @property
    def Alpha(self):
        """Image implementation of regularisation parameter (np.array)"""
        x=np.linspace(-1,1,self.nfx)
        y=np.linspace(-1,1,self.nfy)
        xv, yv = np.meshgrid(x,y)
        R=np.sqrt(np.square(xv) + np.square(yv))
        R=np.fft.fftshift(R)
        if self.alpha[0] > self.alpha[1]:
            # Logistic function instead of error function (to be seen)
            Alpha = self.alpha[0] - ((self.alpha[0] - self.alpha[1]) / (1 + np.exp(-self.alpha_slope * (R-self.alpha_cutoff))))
        elif self.alpha[0] < self.alpha[1]:
            Alpha = self.alpha[0] + ((self.alpha[1] - self.alpha[0]) / (1 + np.exp(-self.alpha_slope * (R-self.alpha_cutoff))))
        else:
            Alpha = self.alpha[0] * R**0 
        return 10**Alpha 

#    def Reconstruct(self, dataset):    
#        #@parallelize(cores) #ou meme pas, should be a configuration thing I suppose
#        parallelizer = Parallelizer.OAR()
#        parallelizer.Launch(dataset, 'retrieve', algorithm=self.__class__.__name__, parameter=self.alpha)
#        pass

#    def ReconstructDifference(self, dataset):    
#        #@parallelize(cores) #ou meme pas, should be a configuration thing I suppose
#        parallelizer = Parallelizer.OAR()
#        parallelizer.Launch(dataset, 'retrieve_difference', algorithm=self.__class__.__name__, parameter=self.alpha)
#        pass

    @Parallelize
    def reconstruct_projections(self, *, dataset, projections):
        """
        Reconstruct a range of projections (parallelized function).
        
        Parameters
        ----------
        dataset : pyphase.Dataset
            Dataset to reconstruct.
        projections : list of int
            In the form [start, end]
        
        """
        
        for projection in range(projections[0], projections[1]+1):
            #print("Projection: {}".format(projection))
            self.reconstruct_projection(dataset=dataset, projection=projection)

    def frequency_variable(self, nfx, nfy, sample_frequency):
        """
        Calculate frequency variables.
        
        Parameters
        ----------
        nfx : int
            Number of samples in x direction
        nfy : int
            Number of samples in y direction
        sample_frequency : float
            Reciprocal of pixel size in 1/m
            
        Returns
        -------
        nparray
            Frequency variables as an array of size  [nfy, nfx, 2]
            
        Notes
        -----
        Follows numpy FFT convention. Zero frequency at [0,0], [1:n//2] 
        contain the positive frequencies, [n//2 + 1:] n the negative 
        frequencies in increasing order starting from the most negative 
        frequency.
        """
        
        if type(sample_frequency) == int:
            sample_frequency = np.array([sample_frequency, sample_frequency]) #TODO: refactor
        
        x=0
        x=np.append(x, np.linspace(sample_frequency[0]/nfx, sample_frequency[0]/2, nfx//2))
        x=np.append(x, np.linspace(-sample_frequency[0]/2+sample_frequency[0]/nfx, -sample_frequency[0]/nfx, int(nfx//2-1+(np.ceil(nfx/2)-nfx//2))))
        
        y=0
        y=np.append(y, np.linspace(sample_frequency[1]/nfy, sample_frequency[1]/2, nfy//2))
        y=np.append(y, np.linspace(-sample_frequency[1]/2+sample_frequency[1]/nfy, -sample_frequency[1]/nfy, int(nfy//2-1+(np.ceil(nfy/2)-nfy//2))))
        return np.meshgrid(x, y)

    def simple_propagator(self, pxs, Lambda, z):
        """
        Creates a Fresnel propagator.
        
        Parameters
        ----------
        pxs : float
            Pixel size in µm.
        Lambda : float
            Wavelength in m.
        z : float
            Effective propagation distance in m.
            
        Returns
        -------
        H : nparray
            Fresnel propagator.
        
        Notes
        -----
        Temporary implementation by Y. Zhang. Will be integrated with the
        propagator module.
        
        """
        
        # TODO: need to be replaced by the propagator class
        # TODO: Needs refactoring to meet coding standards
        # Generates the ifftshift version of the propagators
        x = np.arange(-np.fix(self.nfx / 2), np.ceil(self.nfx / 2))
        y = np.arange(-np.fix(self.nfy / 2), np.ceil(self.nfy / 2))
        fx = np.fft.ifftshift(x / (self.nfx * pxs))
        fy = np.fft.ifftshift(y / (self.nfy * pxs))
        Fx, Fy = np.meshgrid(fx, fy)
        f2 = Fx ** 2 + Fy ** 2
        H = np.zeros([self.ND, self.nfy, self.nfx], dtype="complex_")
        for distance in range(self.ND):
            H[distance,:,:] = np.exp(-1j * np.pi * Lambda * z[distance] * f2) 
        return H

    def reconstruct_projection(self, dataset, projection=0, positions=None, pad=True):
        """
        Reconstruct one projection from a Dataset object and saves the result.
        
        Parameters
        ----------
        dataset : Dataset
            Dataset object to use (not necessarily the same as initialised).
        projection : int, optional
            Number of projection to reconstruct.
        positions : int or list of ints, optional
            Subset of positions to use for reconstruction
            
        Returns
        -------
        phase : np.array
            Reconstructed phase.
        attenuation : np.array
            Reconstructed attenuation.        
        """      
                    
        ID = np.zeros((self.ND, self.nfy, self.nfx))
        #dataset.padding = self.padding #TODO: needs proper sorting out...
        for position in range(self.ND):
#                if difference:
#                    FID[position] = dataset.get_image(projection=projection, difference=True, Fourier=True, position=position) #TODO: difference not really necessary, just give the difference image as input...
#                else:
                    ID[position] = dataset.get_projection(projection=projection, position=position, pad=pad)
#       elif image.any():
        
        # Call algorithm part
        phase, attenuation = self.reconstruct_image(ID, positions=positions, pad=pad)
        
        
#            if difference:
#                dataset.write_image(phase, 'phase update', projection)
#            else:
        dataset.write_image(image=Utilities.resize(phase, [self.ny, self.nx]), projection=projection)
        dataset.write_image(image=Utilities.resize(attenuation, [self.ny, self.nx]), projection=projection, projection_type='attenuation')

        return phase, attenuation

    def reconstruct_image(self, image, positions=None, pad=False):
        """
        Template for reconstructing an image given as argument. 
        
        Arguments
        ---------
        image : numpy.array
            A phase contrast image or an ndarray of images stacked along the
            first dimension.
        positions : int or list of ints, optional
            Subset of positions to use for reconstruction.
        
        Note
        ----
        Calls _algorithm (container purely for algorithm part).
        
        Returns
        -------
        phase : numpy.array
        attenuation : numpy.array
        """

        if len(image.shape) == 2: # If only one image is given, add 3rd dimention for compatibility with loops
            image = image[np.newaxis]

        if not positions:
            positions = list(range(self.ND))
        elif not isinstance(positions, list):
            positions = [ positions ]
            
        image = Utilities.resize(image, [self.nfy, self.nfx])
            
        phase, attenuation = self._algorithm(image, positions=positions)
        
        if not pad:
            phase = Utilities.resize(phase, [self.ny, self.nx])
            attenuation = Utilities.resize(attenuation, [self.ny, self.nx])
        
        return phase, attenuation
      
class TIEHOM(PhaseRetrievalAlgorithm2D):
    """
    Transport of Intensity Equation for homogeneous objects (or "Paganin's algorithm") [1]
    
    Parameters
    ----------
    delta_beta : float, optional
        Material dependent ratio delta over beta.
    
    References
    ----------
    [1] Paganin et al. J. Microsc. 206 (2002) 33
    """
        
    def __init__(self, dataset=None, delta_beta=500, **kwargs):        
        self._delta_beta=delta_beta
        self.padding = 2
        super().__init__(dataset, **kwargs)

    @property
    def delta_beta(self):
        """Material dependent ratio delta over beta (float)."""
        return self._delta_beta
    
    @delta_beta.setter
    def delta_beta(self, delta_beta):
        """Recalculates dependent factors on setting (float)."""
        if delta_beta != self._delta_beta:
            self._delta_beta = delta_beta
            self._compute_factors()

    def _compute_factors(self):
        """Calculate TIEHOM factors. Overrides PhaseRetrievalAlgorithm2D."""
        self.TIEHOM_factor = [0 for xxx in range(self.ND)] 
        for distance in range(self.ND):
            self.TIEHOM_factor[distance] = 1 + self.Fresnel_number[distance] * np.pi * self.delta_beta * ((self.fx ** 2) + (self.fy ** 2))

    def _algorithm(self, image, positions=None):
        """Reconstruct one image or a set of images using TIEHOM."""
        #TODO: Needs verification on simpler images

        FID=np.fft.fft2(image)

        numerator_TIEHOM = np.zeros((self.nfy, self.nfx))
        denominator_TIEHOM = numerator_TIEHOM.copy()
        
        if len(positions) == 1:
            phase = 1/2 * self.delta_beta * np.log(np.real(np.fft.ifft2(FID[positions[0]] / self.TIEHOM_factor[positions[0]])))
        else:  
            for position in positions:
                numerator_TIEHOM = numerator_TIEHOM + (self.TIEHOM_factor[position] * FID[position]) 
                denominator_TIEHOM = denominator_TIEHOM + (self.TIEHOM_factor[position])**2 
            phase = 1/2 * self.delta_beta * np.log(np.real(np.fft.ifft2(numerator_TIEHOM / denominator_TIEHOM)))
        
        attenuation = -1/(self.delta_beta) * phase

        return phase, attenuation

class CTFPurePhase(PhaseRetrievalAlgorithm2D):
    """
    Contrast Transfer Function for pure phase objects [1].
                
    References
    ----------
    [1] Cloetens et al. J. Phys. D: Appl. Phys 29 (1996) 133
    """
        
    def __init__(self, dataset=None, **kwargs):
        self.padding = 2
        super().__init__(dataset, **kwargs)
        
    def _algorithm(self, image, positions=None):

        image=np.fft.fft2(image) #Saves memory?
            
        # Generate CTF factors
        # TODO: should possibly be done in constructor
        sinCTFfactor = np.zeros((self.nfy, self.nfx))
        Delta = np.zeros_like(sinCTFfactor)
        
        for position in positions:
            sinCTFfactor = sinCTFfactor + self.sinchirp[position]*image[position]
            Delta = Delta + self.sinchirp[position]*self.sinchirp[position]   
            
        phase = sinCTFfactor / (2*Delta + self.Alpha)
        phase = np.real(np.fft.ifft2(phase))
        
        attenuation = np.zeros_like(phase)
        
        return phase, attenuation
    
class CTF(PhaseRetrievalAlgorithm2D):
    """
    Contrast Transfer Function [1].
    
    References
    ----------
    [1] Cloetens et al. Appl. Phys. Lett. 75 (1999) 2912
    
    """
    
    def __init__(self, dataset=None, **kwargs):
        self.padding = 2
        super().__init__(dataset, **kwargs)

    def _compute_factors(self):
        """Compute specific factors for CTF phase retrieval"""
        super()._compute_factors()
        self.A = np.zeros((self.nfy, self.nfx))
        self.B = self.A.copy()
        self.C = self.A.copy()

        for distance in range(self.ND):
            self.A += self.sinchirp[distance] * self.coschirp[distance]
            self.B += self.sinchirp[distance] * self.sinchirp[distance]
            self.C += self.coschirp[distance] * self.coschirp[distance]

        self.Delta = self.B*self.C - self.A**2

    def __getstate__(self):
        """Includes specific CTF factors"""
        state = super().__getstate__()
        del state['A'], state['B'], state['C'], state['Delta']
        return state

    def _algorithm(self, image, positions=None):
           
        FID=np.fft.fft2(image)
            
        # Generate CTF factors
        # TODO: should possibly be done in constructor
        sinCTFfactor = np.zeros((self.nfy, self.nfx))
        cosCTFfactor = np.zeros_like(sinCTFfactor)
        
        for distance in positions:
            sinCTFfactor = sinCTFfactor + self.sinchirp[distance]*FID[distance]
            cosCTFfactor = cosCTFfactor + self.coschirp[distance]*FID[distance]            
            # TODO: The removal of the delta is not explicit in the paper 
            # but should probably be done
            # s{k}(1,1) -= nf*mf; # remove 1 in real space
            # TODO: verify correct padding
                
        phase = (self.C*sinCTFfactor - self.A*cosCTFfactor) / (2*self.Delta + self.Alpha)
        attenuation = (self.A*sinCTFfactor - self.B*cosCTFfactor) / (2*self.Delta + self.Alpha)

        phase = np.real(np.fft.ifft2(phase))
        attenuation = np.real(np.fft.ifft2(attenuation))
    
        return phase, attenuation
    
class Mixed(PhaseRetrievalAlgorithm2D):
    """Mixed approach phase retrieval
    
    
    Note
    ----
    Legacy code to be aligned with current APIxs
    """
    # Guigay et al. Optics Letters 32, 1617, 2007, Langer et al. TIP 19, 2428, 2010
    # TODO: Regularisation should somehow be separated
    
    def __init__(self, dataset):
        self.padding = 2
        super().__init__(dataset)

        self.delta_beta = 500 #TODO: should there even be a default value? Where should delta_beta live?

        self.sumAD2 = np.zeros((self.nfy, self.nfx))
#            Cf = np.zeros((self.nfy, self.nfx))
#            Cg = np.zeros((self.nfy, self.nfx))
        self.coschirp_dfx = [0 for x in range(self.ND)]
        self.coschirp_dfy = [0 for x in range(self.ND)]
        
        for distance in range(1,self.ND):
            self.sumAD2 = self.sumAD2 + self.sinchirp[distance]*self.sinchirp[distance] # Denominator Eq. 17
            self.coschirp_dfx[distance] = self.coschirp[distance] * 1j*self.fx*self.Fresnel_number[distance] # First part of Delta
            self.coschirp_dfy[distance] = self.coschirp[distance] * 1j*self.fy*self.Fresnel_number[distance] # First part of Delta

        self.R = np.sqrt(np.square(self.fx) + np.square(self.fy))
        self.LP_cutoff = 0.5
        self.LP_slope = .5e3
        self.LPfilter = 1 - 1/(1 + np.exp(-self.LP_slope * (self.R-self.LP_cutoff))) #Logistic filter
        
        self.sigma_I0filter = 10

        self.iterations = 5

        self.prior = 'homogeneous' # 'homogeneous' 'functional' TODO: put in data file
        
    def get_prior(self, projection):
        """Generates a prior estimate of the phase"""
        if self.prior == 'forward':
            self.LPfilter * np.fft.fft2(np.log(I0)*I0/2)
        else:
            self.LPfilter * np.fft.fft2(np.log(I0)*I0/2)
        
    def Lcurve(self, dataset, projection):
        """Calculate the L-curve (for finding regularisation parameter)"""
        Lcurve_min = -9
        Lcurve_max = -3
        Lcurve_step = 1
        Lcurve_range = np.arange(Lcurve_min, Lcurve_max+1, Lcurve_step, dtype=float)
        
        alpha_HF = -10

        FID = [0 for x in range(self.ND)]
        for distance in range(self.ND):
            FID[distance] = dataset.get_projection(projection, distance+1, 'Fourier')              
            if not distance == 0:     
                FID[distance] = FID[distance] - FID[0] #Wouldn't it be better filtered as well?

        FI0_filtered = scipy.ndimage.fourier_gaussian(FID[0], self.sigma_I0filter)
        I0_filtered = np.real(np.fft.ifft2(FI0_filtered))
        dfxI0 = np.real(np.fft.ifft2(2j*np.pi*self.fx * FI0_filtered)) / I0_filtered
        dfyI0 = np.real(np.fft.ifft2(2j*np.pi*self.fy * FI0_filtered)) / I0_filtered

        I0 = np.real(np.fft.ifft2(FID[0]))
        #TODO: refactor doule code with reconstruct_projection

        if self.prior == 'forward':
            prior = dataset.get_image(projection, 'prior', 'Fourier')
            self.delta_beta = 1 #TODO: Necessary?
            print('forward')
        else:
            prior = np.fft.fft2(np.log(I0)*I0/2)

        prior = self.LPfilter * prior

        model_error = np.zeros(Lcurve_range.shape)
        regularisation_error = np.zeros(Lcurve_range.shape)        
        for index in range(Lcurve_range.shape[0]):
            print('Alpha: {}, Delta/Beta: {}'.format(Lcurve_range[index], self.delta_beta))
            self.alpha = np.array([10**Lcurve_range[index], 10**alpha_HF])

            phase = self.reconstruct_projection(dataset, projection) # TODO: reconstruct_projection could return image?
            
            # Propagate with mixed (move to propagator)
            
            phase_dfxI0 = np.fft.fft2(phase*dfxI0)
            phase_dfyI0 = np.fft.fft2(phase*dfyI0)
            
            phasef = np.fft.fft2(phase)
            
            for distance in range(1, self.ND):
                mixed_contrast = 2*self.sinchirp[distance]*phasef + self.coschirp_dfx[distance]*phase_dfxI0 + self.coschirp_dfy[distance]*phase_dfyI0
                model_difference = FID[distance] - mixed_contrast
                model_difference[0, 0] = 0 # Disregard offset (necessary?)
                model_difference_r = np.fft.ifft2(model_difference)
                model_difference_rc = np.real(model_difference_r[self.ny//2:-self.ny//2, self.nx//2:-self.nx//2])
                model_error[index] += np.sum(np.square(model_difference_rc)) / (self.nx*self.ny)
            
            model_error[index] = model_error[index] / (self.ND-1)
            regularisation_difference = phasef - self.delta_beta*prior
            regularisation_difference = np.fft.ifft2(regularisation_difference)
            regularisation_difference = np.real(regularisation_difference[self.ny//2:-self.ny//2, self.nx//2:-self.nx//2])
            regularisation_error[index] = np.sum(np.square(regularisation_difference)) / (self.nx*self.ny)
            print("ME: {} , RE: {}".format(model_error[index], regularisation_error[index]))
        
        model_error_log = np.log10(model_error)
        regularisation_error_log = np.log10(regularisation_error)

#        model_error_log = np.array([-4.2471, -4.2460, -4.2442, -4.2384, -4.2114, -4.1003, -3.8592])
#        regularisation_error_log = np.array([1.477485, 0.790983, 0.313785, -0.039234, -0.403655, -0.896721, -1.472450])
#        LR = np.array([-9, -8, -7, -6, -5, -4, -3])
        
        Loversamp=10
        #LR = np.log10(Lcurve_range)     
        t = np.linspace(0, 1, len(model_error_log))
        ts = np.linspace(0, 1, len(model_error_log)*Loversamp)
        
        M = interpolate.UnivariateSpline(t, model_error_log, s=0)    
        R = interpolate.UnivariateSpline(t, regularisation_error_log, s=0)    
        
        Mp = M.derivative()
        Rp = R.derivative()
        
        Mpp = Mp.derivative()
        Rpp = Rp.derivative()
        
        K = (Mp(ts)*Rpp(ts) - Rp(ts)*Mpp(ts)) / (Rp(ts)**2 + Mp(ts)**2)**1.5 # Langer 2010 eq. 25

        Mts = M(ts)
        Rts = R(ts)
        
#        Kmax = K[Mts.argmin():Rts.argmin()].argmax()
        Kmax = K.argmax()
        Mmin = Mts.argmin()
        
        Lts = np.linspace(Lcurve_range.min(), Lcurve_range.max(), len(Lcurve_range)*Loversamp)
        
        #TODO: Plot function in display

        self.alpha[0]=10**Lts[Kmax]
        dataset.alpha = self.alpha
        
        lcurve_filename=dataset.path+'/'+dataset.name+'_/lcurve.pickle' #TODO: refactor
        with open(lcurve_filename, 'wb') as f:
            pickle.dump([Lts, model_error_log, regularisation_error_log, Mts, Rts, K, Kmax, Mmin], f, pickle.HIGHEST_PROTOCOL)

        dataset.WriteParameterFile()
        self.display_Lcurve(dataset)

    def display_Lcurve(self, dataset):
        """Displays the L-curve"""
        lcurve_filename=dataset.path+'/'+dataset.name+'_/lcurve.pickle' #TODO: refactor
        with open(lcurve_filename, 'rb') as f:
            Lts, model_error_log, regularisation_error_log, Mts, Rts, K, Kmax, Mmin = pickle.load(f)
###
        Loversamp=10
        #LR = np.log10(Lcurve_range)     
        t = np.linspace(0, 1, len(model_error_log))
        ts = np.linspace(0, 1, len(model_error_log)*Loversamp)
        
        M = interpolate.InterpolatedUnivariateSpline(t, model_error_log, k=4)    
        R = interpolate.InterpolatedUnivariateSpline(t, regularisation_error_log, k=4)    
        
        
        Mp = M.derivative()
        Rp = R.derivative()
        
        Mpp = M.derivative(2)
        Rpp = R.derivative(2)
        
        K = (Mp(ts)*Rpp(ts) - Rp(ts)*Mpp(ts)) / (Rp(ts)**2 + Mp(ts)**2)**1.5 # Langer 2010 eq. 25

        Mts = M(ts)
        Rts = R(ts)
###
        
        pyplot.figure()
        pyplot.plot(Lts, K)
        pyplot.show()

        pyplot.figure()
        pyplot.plot(model_error_log, regularisation_error_log, 'rx', Mts, Rts, 'b-', Mts[Kmax], Rts[Kmax], 'go', Mts[Mmin], Rts[Mmin], 'co')
        pyplot.show()
       
#        pyplot.figure()
#        pyplot.title("Mts")
#        pyplot.plot(Lts,Mts)
#        pyplot.show()
#        
#        pyplot.figure()
#        pyplot.title("Rts")
#        pyplot.plot(Lts,Rts)
#        pyplot.show()
#        
#        pyplot.figure()
#        pyplot.title("Mp")
#        pyplot.plot(Lts,Mp(ts))
#        pyplot.show()
#
#        pyplot.figure()
#        pyplot.title("Mpp")
#        pyplot.plot(Lts,Mpp(ts))
#       
#       pyplot.show()
#
#        pyplot.figure()
#        pyplot.title("Rp")
#        pyplot.plot(Lts,Rp(ts))
#        pyplot.show()
#
#        pyplot.figure()
#        pyplot.title("Rpp")
#        pyplot.plot(Lts,Rpp(ts))
#        pyplot.show()

        print('Maximum curvature at: {}'.format(Lts[Kmax]))
        print('Minimum model error at: {}'.format(Lts[Mmin]))

        pass
                
    def _algorithm(self, image, positions=None):
        FID = np.fft.fft2(image)
        for distance in positions:
            if not distance == 0: #TODO: Handle case when position_number makes sense
                FID[distance] = FID[distance] - FID[0] #Wouldn't it be better filtered as well?
 
        # TODO: Need a utility function for filters        
        # TODO: Quick and dirty implementation to be refactored
        
        I0 = np.real(np.fft.ifft2(scipy.ndimage.fourier_gaussian(FID[0], self.sigma_I0filter)))
        
        # Gradient of attenuation image. 
        dfxI0 = np.real(np.fft.ifft2(2j*np.pi*self.fx * scipy.ndimage.fourier_gaussian(FID[0], self.sigma_I0filter))) / I0
        dfyI0 = np.real(np.fft.ifft2(2j*np.pi*self.fy * scipy.ndimage.fourier_gaussian(FID[0], self.sigma_I0filter))) / I0

        # TODO: I guess priors etc should be in ther own classes/functions? How to handle the Lcurve case...
        # Estimate of phase*absorption with delta/beta = 1
        if self.prior == 'forward':
            prior = dataset.get_image(projection, 'prior', 'Fourier')
            self.delta_beta = 1 #TODO: Necessary?
            print('forward')
        else:
            prior = np.fft.fft2(np.log(I0)*I0/2)

        prior = self.LPfilter * prior

        phase = np.zeros((self.nfy, self.nfx))
        
        for n in range(self.iterations):
            nominator_term = np.zeros((self.nfy, self.nfx))
            
            phase_dfxI0 = np.fft.fft2(phase*dfxI0)
            phase_dfyI0 = np.fft.fft2(phase*dfyI0)

            for distance in positions[1:]:
                nominator_term = nominator_term + self.sinchirp[distance] * (FID[distance] - self.coschirp_dfx[distance]*phase_dfxI0 - self.coschirp_dfy[distance]*phase_dfyI0)
                
            nominator_term = nominator_term / (self.ND-1)
          
            phase_n = (nominator_term + (self.Alpha*self.delta_beta*prior)) / (self.Alpha+self.sumAD2)
            phase_n = np.real(np.fft.ifft2(phase_n))
            
            print("Iteration: {} RMS: {}".format(n, np.sqrt(np.sum((phase_n-phase)*(phase_n-phase).conjugate())/(self.nfx*self.nfy))))
            phase = phase_n
            
        attenuation = np.zeros_like(phase)

        return phase, attenuation

    def create_multimaterial_prior(self, data):
        """Generate a multi-material prior from a tomographic reconstruction of a contact plane scan
        
        
        Note
        ----
        Legacy code to be refactored.
        """
        # Reconstruct attenuation (if not already reconstructed)
        
        #TODO: Put parameters in parameter file        
        threshold = 2.1
        delta_beta_soft = 1938
        delta_beta_hard = 310
        delta_beta_tmp = [2480, 711, 132] #TODO: refactor
        median_size=2 #TODO: How to 3d filter memmaps? I guess one could parallelize proper
        
        # Create/read prior volume
        print('Creating multi-material prior: db_soft: {}, db_hard: {}, threshold: {}'.format(delta_beta_soft, delta_beta_hard, threshold))
        segmentation_filename = '/mntdirect/_data_id19_bones01/bones3/max/holodata/knee/lbtoKneeWTOA8weekFeb15/18_OA/18_OA_8weeks_1_slice_pag_db0250_1400_Seg.raw' #TODO: Refactor
        #segmentation_filename = ''
        attenuation_volume = np.memmap(data.attenuation_vol_filename, dtype=np.float32, mode='r')
        prior = np.memmap(data.prior_vol_filename, dtype=np.float32, mode='w+', shape=attenuation_volume.shape)

        if segmentation_filename:
            print('Using pre-segmented volume')
            segmentation_volume = np.memmap(segmentation_filename, mode='r', dtype=np.uint8, shape=attenuation_volume.shape)
            #prior[:] = np.piecewise(segmentation_volume[:].astype(np.float32, copy=False), [segmentation_volume < 128, segmentation_volume == 127, segmentation_volume > 127], [delta_beta_tmp[0], delta_beta_tmp[1], delta_beta_tmp[2]])
            #del segmentation_volume, prior
            print('Median filtering attenuation')
            prior = np.memmap(data.prior_vol_filename, dtype=np.float32, mode='w+', shape=attenuation_volume.shape)
            ndimage.filters.median_filter(attenuation_volume, size=3, mode='nearest', output=prior)
            prior.flush()
            prior[:] = -prior[:] * np.piecewise(segmentation_volume[:].astype(np.float32, copy=False), [segmentation_volume < 128, segmentation_volume == 127, segmentation_volume > 127], [delta_beta_tmp[0], delta_beta_tmp[1], delta_beta_tmp[2]]) / 2
            prior.flush()
            
        else:            
            prior[:] = -attenuation_volume[:] * np.piecewise(attenuation_volume[:], [attenuation_volume < threshold, attenuation_volume >= threshold], [delta_beta_soft, delta_beta_hard]) / 2

        #TODO: Interface for selecting parameters?
        #TODO: Read segmented volume
        #TODO: functional

        del prior, attenuation_volume
        
class HIO_ER(PhaseRetrievalAlgorithm2D):
    """
    Sequence of Hybrid Input Output [1] and Error Reduction [2].
    
    Attributes
    ----------
    retriever : PhaseRetrievalAlgorithm2D
        Algorithm for initialisation
    iterations : int
        Number of global iterations
    iterations_hio : int
        Number of HIO iterations per iteration
    iterations_er : int
        Number of ER iterations per iteration
    step_size_phase : float
        Update step size for the phase
    step_size_attenuation : float
        Update step size for the attenuation
    
    
    References
    ----------
    [1] Fienup Appl. Opt. 21 (1982) 2758
    [2] Gerchberg & Saxton Optik 35 (1972) 237   
    """
    
    #TODO: Needs refactoring to conform to coding standards
    
    def __init__(self,dataset=None, **kwargs):
        self.padding=1
        self.retriever = CTFPurePhase(dataset, **kwargs)

        super().__init__(dataset, **kwargs)
            
        self.iterations = 4 # 4 in Yuhe code
        self.iterations_hio = 45 # 45 in Yuhe code
        self.iterations_er = 5 # 5 in Yuhe code
        self.step_size_phase = 0.2
        self.step_size_attenuation = 0.2
        #self.retriever=CTF(dataset=dataset)
        
        self.propagator = self.simple_propagator(self.pixel_size[0], self.Lambda, self.distance)  # TODO: But then it's passed around anyway   
        
    def reconstruct_projection(self, dataset, projection=0, positions=None, pad=False):
        super().reconstruct_projection(dataset=dataset, projection=projection, positions=positions, pad=pad)
        
    @property
    def retriever(self):
        return(self._retriever)

    @retriever.setter
    def retriever(self, retriever):
        self._retriever = retriever
        self._retriever_class = self._retriever.__class__
        
    @property
    def alpha(self):
        return self._alpha
    
    @alpha.setter    
    def alpha(self, value):
        self.retriever.alpha = value
        self._alpha = value
        
    def __getstate__(self):
         """Includes HIO specific variables """
         state = self.__dict__.copy()
         del state['fx'], state['fy'], state['sinchirp'], state['coschirp'], state['retriever'], state['propagator'] #TODO: Needs proper separation PRAlg2d
         return state

    def __setstate__(self, state):
        """Includes HIO specifics (retriever for initialisation notably)"""
        self.__dict__.update(state)
        self._compute_factors()
        self.retriever = self._retriever_class(dataset=self.dataset) # only dataset, reconstruct_image not parallellisable ?
        self.propagator = self.simple_propagator(self.pixel_size[0], self.Lambda, self.distance) #TODO: update to use propagator module

    def amplitude_constraint(self, wavefront, amplitude, propagator, mask=[]):
        """Apply amplitude constraint.
        
        Parameters
        ----------
        wavefront : complex np.array
            Wavefront to constrain.
        amplitude : np.array
            Amplitude to impose.
        propagator : complex np.array
            Propagator corresponding to effective distance of amplitude.
        mask : np.array, optional
            Zone to apply constraint.
        
        Returns
        -------
        wavefront_constrained : complex np.array
            Wavefront after applied constraint.
        
        """
        #TODO: Proper handling of padding
        if mask == []:
            mask = np.ones_like(amplitude)
        wavefront_aux = np.fft.ifft2(np.fft.fft2(wavefront) * propagator) # TODO: Should be done with propagator instead 
        wavefront_aux = np.where(mask != 0, amplitude * np.exp(1j * np.angle(wavefront_aux)), wavefront_aux) # Apply amplitude constraint
        wavefront_constrained = np.fft.ifft2(np.fft.fft2(wavefront_aux) * np.conj(propagator))
        return wavefront_constrained

    def hybrid_input_output(self, wavefront, initial_wavefront, support, step_size_attenuation, step_size_phase):
        """
        One iteration of the Hybrid Input Output algorithm.
        
        Parameters
        ----------
        wavefront : complex np.array
            Constrained wavefront.
        initial_wavefront : complex np.array
            Wavefront to update.
        support : np.array
            Support of object.
        step_size_attenuation : float
            Step size for attenuation update.
        step_size_phase : float
            Step size for phase update.
        
        Returns
        -------
        wavefront_updated : complex np.array
            Updated wavefront.
        
        """
        # TODO: Refactor spaghetti code, make step sizes members. Make methods reasonable to use from outside instead of arbitrary cut up blocks.
        # phase constraint
        # TODO: Should be split out into its own class
        phase = np.angle(wavefront)
        phase = np.where(support == 0, np.angle(initial_wavefront) - step_size_phase * np.angle(wavefront), phase)
        # abs constraint
        attenuation = np.where(np.abs(wavefront) < 1, np.abs(wavefront), 1)
        attenuation = np.where(np.abs(wavefront) > 0, attenuation, 0)
        attenuation = np.where(support == 0, np.abs(initial_wavefront) - step_size_attenuation * (np.abs(wavefront) - 1), attenuation)
        wavefront_updated = attenuation * np.exp(1j * phase)
        return wavefront_updated

    def error_reduction(self, wavefront, support):
        """
        One iteration of Error Reduction.
        
        Parameters
        ----------
        wavefront : complex np.array
            wavefront to update.
        support : np.array
            Support constraint.
            
        Returns
        -------
        wavefront_updated : complex np.array
            Updated wavefront.
        """
        # phase constraint
        phase = np.where(np.angle(wavefront) < 0, np.angle(wavefront), 0)
        phase = np.where(support == 0, 0, phase)
        # abs constraint
        attenuation = np.where(np.abs(wavefront) < 1, np.abs(wavefront), 1)
        attenuation = np.where(np.abs(wavefront) > 0, attenuation, 0)
        attenuation = np.where(support == 0, 1, attenuation)
        wavefront_updated = attenuation * np.exp(1j * phase)
        return wavefront_updated

    def error_estimate(self, wavefront, amplitude, propagator, mask=[]):
        """
        Estimate fit to data.
        
        Parameters
        ----------
        wavefront : complex np.array
            Wavefront for estimation.
        amplitude : np.array
            Amplitude from measured image.
        propagator : complex np.array
            Fresnel propagator corresponding to effective propagation distance
            in measured image.
        mask : np.array
            Restrict estimate to a region of interest.
            
        Returns
        -------
        error : float
            MSE calculated and measured amplitude.
        
        """
        
#        if mask == []:
#            mask = np.ones_like(amplitude)
#        mask = np.where(mask != 0, 1, 0)
#        N = np.sum(mask[:])
        wavefront_aux = np.fft.ifft2(np.fft.fft2(wavefront) * propagator)
        aux = (np.abs(wavefront_aux) - amplitude) ** 2
        error = np.sum(aux[:]) / (self.nfx*self.nfy)
        return error

    def _algorithm(self, image, positions=None, support=None):
        phase, attenuation = self.retriever.reconstruct_image(image) 
        amplitude = np.sqrt(image)

        #FID=np.fft.fft2(image)
      
#        phase = Utilities.resize(phase, (self.ny, self.nx))
#        attenuation = Utilities.resize(attenuation, (self.ny, self.nx))
        
        adjust_offset = True
        if adjust_offset:
            phase = phase - phase.max()

        if not support:
            support = np.ones_like(image[0])

        initial_guess = np.exp(-attenuation) * np.exp(1j * phase) 
        #mask = np.ones_like(initial_guess)
        amplitude = np.sqrt(image)
        wavefront = initial_guess  # Prepare wavefront (ML: What is this supposed to do?)
        #mask = np.where(mask != 0, 1, 0)  # Prepare mask TODO: This doesńt actually do anything?
        support = np.where(support != 0, 1, 0)  # Prepare support
        
        reconstruction = np.empty_like(image, dtype="complex_")
        for distance in positions:
            print(F'========== processing distance {distance+1} ==========')
            error_count = 0
            for ii in range(self.iterations):
                for jj in range(self.iterations_hio + self.iterations_er):
                    initial_wavefront = wavefront
                    # Amplitude constraint
                    wavefront = self.amplitude_constraint(wavefront, amplitude[distance], self.propagator[distance])
                    if jj < self.iterations_hio:
                        wavefront = self.hybrid_input_output(wavefront, initial_wavefront, support, self.step_size_attenuation, self.step_size_phase)
                    else:
                        wavefront = self.error_reduction(wavefront, support)
                    error_count += 1
                    error_it = self.error_estimate(wavefront, amplitude[distance], self.propagator[distance])
                    print('Iteration {:04d}, error: {:0.2g}'.format(error_count,np.real(error_it))) #TODO: print every interation?
#            object = np.fft.fftshift(self.phi) #ML: Why fftshift?!
#            object = wavefront #ML: Why fftshift?!
            reconstruction[distance] = wavefront
        
        phase = np.average(np.angle(reconstruction),axis=0)
        attenuation = np.average(np.abs(reconstruction),axis=0)
                 
        return phase, attenuation

class GradientDescent(PhaseRetrievalAlgorithm2D):
    """Gradient descent algorithm.

    Parameters
    ----------
    PSF: ndarray with the same shape as images
        Point spread function (optional)
    """
    
    def __init__(self, dataset=None, PSF=[], **kwargs):
        self.padding=2
        self.retriever = CTF(dataset, **kwargs)
        super().__init__(dataset, **kwargs)
        # self.dataset = dataset
        self.PSF = PSF
        self.step_size = 0.5
        self.iterations = 20
        self.propagator = self.simple_propagator(self.pixel_size[0], self.Lambda, self.distance)

    #TODO: iterative algorithms could maybe have a base class?        
    @property
    def retriever(self):
        return(self._retriever)

    @retriever.setter
    def retriever(self, retriever):
        self._retriever = retriever
        self._retriever_class = self._retriever.__class__

    @property
    def alpha(self):
        return self._alpha
    
    @alpha.setter    
    def alpha(self, value):
        self.retriever.alpha = value
        self._alpha = value     

    def __getstate__(self):
         """Includes HIO specific variables """
         state = self.__dict__.copy()
         del state['fx'], state['fy'], state['sinchirp'], state['coschirp'], state['retriever'], state['propagator'] #TODO: Needs proper separation PRAlg2d
         return state

    def __setstate__(self, state):
        """Includes HIO specifics (retriever for initialisation notably)"""
        self.__dict__.update(state)
        self._compute_factors()
        self.retriever = self._retriever_class(dataset=self.dataset) # only dataset, reconstruct_image not parallellisable ?
        self.propagator = self.simple_propagator(self.pixel_size[0], self.Lambda, self.distance) #TODO: update to use propagator module

    def _algorithm(self, image, positions=False):
        # TODO: Correct handling of inoutting images
        # TODO: Actually, this is steepest descent, not CG
        #FID=np.fft.fft2(image)
 
        print('Initialising')
        phase, attenuation = self.retriever.reconstruct_image(image, pad=True) 
        
        adjust_offset = True
        if adjust_offset:
            phase = phase - phase.max()

        initial_guess = np.exp(-attenuation) * np.exp(1j * phase) 

        object_FT = np.fft.fft2(initial_guess)

        current_object = object_FT

        for iteration in range(self.iterations):
            error = 0
            print(F'Iteration {iteration + 1} of {self.iterations}')
            for distance in range(self.ND):
                field = np.fft.ifft2(self.propagator[distance, :, :] * object_FT)
                intensity_calculated = np.real(field * np.conj(field))
                if self.PSF != []:  
                    intensity_calculated = np.real(np.fft.ifft2(np.fft.fft2(intensity_calculated) * self.PSF))
                intensity_calculated = intensity_calculated / np.sum(intensity_calculated[:]) * self.nfx * self.nfy
                intensity_difference = image[distance] - intensity_calculated
                error = error + np.std(intensity_difference * intensity_difference) / self.ND
                update = np.conj(self.propagator[distance]) * np.fft.fft2(intensity_difference * field) #TODO: This just can be right
                current_object = current_object + (self.step_size / self.ND) * update
            object_FT = current_object
            print(F'Error for iteration {iteration + 1} is {error:.4g}')
            field = np.fft.ifft2(object_FT)
            phase = np.angle(field)
            attenuation = np.abs(field)
        
        return phase, attenuation

# class Iterative():
#     def __init__(self):
#         pass
    
#     def Reconstruct(self, dataset, iterations, parameter, options=''):
#         #initialise
#         print('Starting iterative reconstruction, {} iterations'.format(iterations))
#         retriever = CTF(dataset)
#         retriever.alpha=parameter
#         propagator = Propagator.CTF(dataset)
#         #TODO: I guess these are the kind of things that should be configured rather
#         tomography = Tomography.PyHST()
        
#         if not 'no_reinitialisation' in options:
#             print('initialising')
#             #TODO: I guess tomo should be without filter?
#             retriever.Reconstruct(dataset)
#             tomography.Reconstruct(dataset, volume='phase')
#             tomography.Reconstruct(dataset, volume='retrieved_attenuation')
        
#         #TODO: should have a convergence criteron also
#         for iteration in range(iterations):
#             print('----- Iteration {} / {} -----'.format(iteration+1, iterations))
#             tomography.ForwardProject(dataset, 'phase')
#             tomography.ForwardProject(dataset, 'retrieved_attenuation')
#             propagator.Propagate(dataset)
#             dataset.Difference()
#             retriever.ReconstructDifference(dataset)
#             tomography.Reconstruct(dataset,volume='update')    
#             tomography.Reconstruct(dataset,volume='attenuation_update')    
#             dataset.UpdatePhase()
#             dataset.UpdateAttenuation()
