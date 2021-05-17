# -*- coding: utf-8 -*-

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

import sys

import vendor.pyelastix as PyElastix
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import scipy.ndimage as ndimage


# update() : Displays or updates a console progress bar
## Accepts a float between 0 and 1. Any int will be converted to a float.
## A value under 0 represents a 'halt'.
## A value at 1 or bigger represents 100%
def update(title, position, target):
    barLength = 20 # Modify this to change the length of the progress bar
    status = ""
    if position < 0:
        progress = 0
        status = "Halt...\r\n"
    if position >= target:
        progress = 1
        status = "Done...\r\n"
    else:
        progress = position / target
    block = int(round(barLength*progress))
    text = "\r" + title + ": [{}] {} / {} {:0.2%} {}".format( "#"*block + 
        "-"*(barLength-block), position, target, progress, status)
    sys.stdout.write(text)
    sys.stdout.flush()


def resize(image, shape):
    """
    Resizes an image by either cutting out the centre or padding.
    
    Assumes images are stored along fist dimension.
    """
    #TODO: special case of one dimension smaller and one larger not covered atm
    adjust_axis = False #TODO: To handle 2 or 3 dimensions. can probably be refactored
    if len(image.shape) == 2: 
        image=image[np.newaxis]
        adjust_axis = True
    
    if  image.shape[1] < shape[0]:
        pady=(shape[0]-image.shape[1])//2
        padx=(shape[1]-image.shape[2])//2
        image = np.pad(image, ((0, 0), (pady,shape[0]-image.shape[1]-pady), (padx, shape[1]-image.shape[2]-padx)), 'edge')
    elif image.shape[1] > shape[0]:
        cy = (image.shape[1]-shape[0])//2
        cx = (image.shape[2]-shape[1])//2
        image = image[:, cy:cy+shape[0], cx:cx+shape[1]]
    else:
        pass
    if adjust_axis:
        return image[0]
    else:
        return image

class RegistrationAlgorithm:
    """
    Abstract class for registration algorithms
    
    Properties
    ----------
    """
    def __init__(self):
        pass
    
    def register(self, moving_image, stationary_image):
        """
        Register moving_image to stationary_image.
        
        Parameters
        ----------
        moving_image : ndarray
            The image to register.
        stationary_image : ndarray
            The image to register to.
            
        Returns
        -------
        field : ndarray
            The calculated deformation field.
        transformed_moving_image : ndarray
            The deformed moving image.
        transform_parameters : array
            The calculated transform parameters. Lenght varies with the number
            of parameters in the chosen algorithm (number_of_parameters)
        
        """
        transformed_moving_image, field, transform_parameters = PyElastix.register(np.ascontiguousarray(moving_image), np.ascontiguousarray(stationary_image), self.parameters, verbose=0)
        return field, transformed_moving_image, transform_parameters
    pass

    def apply_transformation(self, image, transform_parameters, **kwargs):
        """ 
        Apply an image transform from image registration.
        
        Parameters
        ----------
        image : ndarray
            The image to transform.
        transform_parameters : array
            Parameters of the transform (length depends on registration
            algorithm used (number_of_parameters))
            
        Returns
        -------
        transformed_image
            The transformed image.
        
        """
        if len(transform_parameters) == 4: # Similarity transform
            transformed_image = ndimage.zoom(image, 1/transform_parameters[0], mode='nearest')
            transformed_image = ndimage.rotate(transformed_image, np.degrees(transform_parameters[1]), mode='nearest')
            transformed_image = ndimage.shift(transformed_image, [-transform_parameters[3], -transform_parameters[2]], mode='nearest')
        elif len(transform_parameters) == 3: # Euler/Rigid
            transformed_image = ndimage.rotate(image,np.degrees(transform_parameters[0]), mode='nearest')
            transformed_image = ndimage.shift(transformed_image, [-transform_parameters[2], -transform_parameters[1]], mode='nearest')

        return transformed_image

class ElastixRigid(RegistrationAlgorithm):
    """
    Rigid registration algorithm using Elastix.
    
    Attributes
    ----------
    parameters : Parameters
        Elastix standard parameters
    number_of_parameters : int, default=3
        Number of parameters in the transform
        
    """
    def __init__(self):
        self.parameters = PyElastix.get_default_params(type='RIGID')
        
        #self.params.Metric = 'NormalizedMutualInformation'
        #self.params.Metric = 'AdvancedMeanSquares'
        self.parameters.NumberOfResolutions = 6
        self.parameters.MaximumNumberOfIterations = 500
        self.number_of_parameters=3        

class ElastixAffine(RegistrationAlgorithm):
    """
    Affine registration algorithm using Elastix.
    
    Attributes
    ----------
    parameters : Parameters
        Elastix standard parameters
    number_of_parameters : int, default=6
        Number of parameters in the transform
    """
    def __init__(self):
        self.parameters = PyElastix.get_default_params(type='AFFINE')
        
        #self.params.Metric = 'NormalizedMutualInformation'
        #self.params.Metric = 'AdvancedMeanSquares'
        self.parameters.MaximumNumberOfIterations = 500
        self.number_of_parameters=6

class ElastixSimilar(RegistrationAlgorithm):
    """
    Similarity transform registration algorithm using Elastix
    
    Attributes
    ----------
    parameters : Parameters
        Elastix standard parameters
    number_of_parameters : int, default=4
        Number of parameters in the transform    
    """
    def __init__(self):
        self.parameters = PyElastix.get_default_params(type='SIMILAR')
#        self.params.NumberOfResolutions = 8
#        self.params.MaximumNumberOfIterations = 3000
        self.parameters.NumberOfResolutions = 6
        self.parameters.MaximumNumberOfIterations = 500
        #self.params.SP_A = 1
        #self.params.AutomaticParameterEstimation = True
        #self.params.AutomaticTransformInitialization = False
        #self.params.AutomaticScalesEstimation = False
        self.number_of_parameters=4
        #print(self.params.as_dict())
                    
class ImageDisplayer:
    """
    Wrapper class for image display.
    """
    # Should be an abstract class, right?!
    # Which should be inherited by the class for each displayer, right?!

    def __init__(self):
        
        pass

#    def displayImage(self, DS, n, d, correct):
#        pass
# TODO: one viewer class
class PyplotImageDisplayer(ImageDisplayer):
    """
    Interface to Pyplot for image display.
    
    Notes
    -----
    With the idea to make the choice of display package flexible.
    """
    # TODO: rename Pyplot? Develop interface ImageJ?

    def display(self, image, title='', vmin=None, vmax=None):
        '''
        Display an image
        
        Parameters
        ----------
        image : nparray
            The image to be displayed.
        title : str, optional
            Title of figure.
        vmin : optional
            Lower limit of contrast range.
        vmax : optional
            Upper limit of contrast range.
        '''
        # TODO: should include all the different possibilities I presume. Should be kwargs, not positional
        fig, ax = plt.subplots()
        fig.suptitle(title)
        im = ax.imshow(image, cmap='gray', vmin=vmin, vmax=vmax)
        ax.set_axis_off()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        
        fig.colorbar(im, cax=cax)
        
        im.axes.figure.canvas.draw()
        plt.show()
        
    def display_stack(self, stack):
        # TODO: merge into display
        fig, ax = plt.subplots(1, 1)
        viewer = StackViewer(stack, ax)
        fig.canvas.mpl_connect('key_press_event', viewer._on_key)
        plt.show()
    def close_all(self):
        plt.close('all')
        
        
# maybe should be displayShifts?
class StackViewer(object):
    """
    Functionality to browse stacks.
    """
    #TODO: Working?! Should be with an _ to not be imported (internal functionality?)
    def __init__(self, X, ax):

        self.ax = ax
        self.ax.set_title('use left/right arrow to navigate images')

        self.X = X
        self.slices, rows, cols  = X.shape
        self.ind = 0

        self.im = self.ax.imshow(self.X[self.ind, :, :], cmap='gray')
        self._update()

    def _on_key(self, event):
        if event.key == 'right':
            self.ind = (self.ind + 1) % self.slices
        elif event.key == 'left':
            self.ind = (self.ind - 1) % self.slices
        self.update()

    def _update(self):
        self.im.set_data(self.X[:, :, self.ind])
        self.ax.set_ylabel('image %s' % self.ind)
        self.im.axes.figure.canvas.draw()

# 