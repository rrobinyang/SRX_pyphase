Multidistance CTF Phase Retreival 
---------------------------------

Here is an example on how to use the ``phaseretreival.CTF`` to retreive phase and attenuation from a multiple distance projections.



:ref:`install` pyPhase and then:

.. code:: python

   import pyphase

Define the projects name

.. code:: python

   name = 'holosim'

and set the path to the data

.. code:: python
   
   path = '/data/staff/tomograms/HoloSim'

Make sure that the data folders and files have the right [format and structure](structure.md).

Create a DataSet object

.. code:: python

   ds = dataset.DataSet(path,name)


Choose a phase retreival algorithm(link to list of algorithms) and create a retreiver for the dataset, with the appropriate arguments 

.. code:: python
   
   alpha =  1e-8
   retriever = phaseretrieval.CTF(ds,alpha)

Select to range to of projections for which to run the retrieval algorithm

.. code:: python

    start = 0   #default = 0
    end = 3     #default = 0

Run the phase retreival

.. code:: python

   retriever.ReconstructProjections(ds,start, end)


For each projection, a file named `holosim_[version_predix]_000[n].edf` with the retrieved phase and and file named `holosim_att_`, with the retrieved attenuation, is created in the `[path]/myProject_`

