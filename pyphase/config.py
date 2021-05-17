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

from pathlib import Path
import sys, importlib
import pyphase.utilities as Utilities
import pyphase.tomography as Tomography
import pyphase.parallelizer as Parallelizer
importlib.reload(sys.modules['pyphase.parallelizer'])

pyphase_path = Path('./')
viewer = Utilities.PyplotImageDisplayer()
registrator = Utilities.ElastixSimilar()
tomography = Tomography.PyHST()
default_method = 'CTF'
Parallelize = Parallelizer.SLURM
number_of_cores = 24
number_of_nodes = 2 #TODO: Should be in config
memory_per_node = 102000
min_memory_per_core = 6000
memory_per_core = memory_per_node/number_of_cores
