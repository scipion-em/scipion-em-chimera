======================
Chimera scipion plugin
======================

This plugin allows to use chimeraX commands within the Scipion framework.

Chimera  is a program for interactive visualization and analysis of molecular structures and related data. It is developed by the Resource for Biocomputing, Visualization, and Informatics (see `ChimeraX home page <https://www.cgl.ucsf.edu/chimerax/>`_ for details).


===================
Install this plugin
===================

You will need to use `3.0.0 <https://scipion-em.github.io/docs/release-3.0.0/>`_ version of Scipion to run these protocols. To install the plugin, you have two options:

- **Stable version**  

.. code-block:: 

      scipion3 installp -p scipion-em-chimera
      
OR

  - through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**
      
- **Developer's version**

1. Download repository:

.. code-block::

            git clone https://github.com/scipion-em/scipion-em-chimera.git

2. Install:

.. code-block::

            scipion3 installp -p path_to_scipion-em-chimera --devel

- **Binary files**

Chimera binaries could be installed automatically with the plugin after accepting ChimeraX licence terms,
but you can also link an existing installation. Default installation path assumed is *software/em/chimerax-1.0,
if you want to change it, set *CHIMERA_HOME* in *scipion.conf* file to the folder where ChimeraX is installed
or link your chimerax folder to *software/em/chimerax-1.0*.

- **Tests**

Tested with ChimeraX version: 1.X. where X is in the range 0-6

To check the installation, simply run the following Scipion tests: 

* scipion test chimera.tests.test_protocol_chimera_operate
* scipion test chimera.tests.test_protocol_chimera_fit
* scipion test chimera.tests.test_protocol_modeller_search
* scipion test chimera.tests.test_protocol_contact
* scipion test chimera.tests.test_protocol_chimera_map_subtraction

NOTE: the fact that chimerax-plugins supports the above mentioned chimerax versions does not
mean that chimerax supports your operating system. For example, if you want to install 
chimerax 1.4 or higher you need at least ubuntu 22 version.


=========
Protocols
=========

* rigit fit: Fits an atomic structure (PDB/CIF) to a 3D map using a rigid fit algorithm.
* model from template: Models three-dimensional structures of proteins using `Modeller <https://salilab.org/modeller/manual/node7.html>`_.
* contacts: Identifies interatomic `clashes and contacts <https://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/findclash/findclash.html>`_ based on van der Waals radii. 
* operate: Provides access to Chimera and allows to save the result in Scipion framework.
* restore: This protocol opens Chimera and restores the session previously saved with commands *scipionwrite* or *scipionss*. 
* map subtraction: Protocol to identify remnant densities in a density map by subtracting two maps or masking one of them.
* alphafold prediction: finds and retrieves existing models from the AlphaFold Database, runs new AlphaFold predictions using Google Colab, executes a local implementation of alphafold. 


========
Examples
========

See `Model Building Tutorial <https://scipion-em.github.io/docs/release-3.0.0/docs/user/user-documentation.html#model-building>`_

===
FAQ
===

 check the FAQ <https://github.com/scipion-em/scipion-em-chimera/blob/devel/FAQ.rst> for known problems


===============
Buildbot status
===============

Status devel version:

.. image:: http://scipion-test.cnb.csic.es:9980/badges/chimera_devel.svg

..
    Status production version: 

.. 
    image:: http://scipion-test.cnb.csic.es:9980/badges/chimera_prod.svg


======================
Chimera Extra commands
======================
A set of comamnd to allow interaction between scipion and chimerax has been implemented.
They may be executed from ChimeraX command line:
  
* scipionwrite: saves model file
* scipionss: saves chimera session
* scipionrs: restores chimera sesion
* scipioncombine: combines two models [No longer needed since chimerax now implements combine]
* scipionshell: shows the density of a 3D Map on a spherical shell of the map at a given radius
* scipionshellcrown: shows a shell of a 3D map between two radii
* scipion: summary with all scipion related commands
