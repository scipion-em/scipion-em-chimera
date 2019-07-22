================================
Chimera scipion plugin
================================

This plugin allows to use chimera commands within the Scipion framework.

Chimera  is a program for interactive visualization and analysis of molecular structures and related data. It is developed by the Resource for Biocomputing, Visualization, and Informatics (see `Chimera home page <https://www.cgl.ucsf.edu/chimera/>`_ for details).


===================
Install this plugin
===================

You will need to use `2.0.0 <https://github.com/I2PC/scipion/releases/tag/v2.0>`_ version of Scipion to run these protocols. To install the plugin, you have two options:

- **Stable version**  

.. code-block:: 

      scipion installp -p scipion-em-chimera
      
OR

  - through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**
      
- **Developer's version** 

1. Download repository: 

.. code-block::

            git clone https://github.com/scipion-em/scipion-em-chimera.git

2. Install:

.. code-block::

            scipion installp -p path_to_scipion-em-chimera --devel

- **Binary files** 

Chimera binaries will be installed automatically with the plugin, but you can also link an existing installation. 
Default installation path assumed is *software/em/chimera-1.13.1*, if you want to change it, set *CHIMERA_HOME* in *scipion.conf* file to the folder where Chimera is installed.

- **Tests**

Tested with Chimera version: 1.13.1.

To check the installation, simply run the following Scipion tests: 

* scipion test chimera.tests.test_protocol_chimera_operate
* scipion test chimera.tests.test_protocol_chimera_fit
* scipion test chimera.tests.test_protocol_modeller_search
* scipion test chimera.tests.test_protocol_contact

- **Supported versions of Chimera**

1.13.1


=========
Protocols
=========

* rigit fit: Fits an atomic structure (PDB/CIF) to a 3D map using a rigid fit algorithm.
* model from template: Models three-dimensional structures of proteins using `Modeller <https://salilab.org/modeller/manual/node7.html>`_.
* contacts: Identifies interatomic `clashes and contacts <https://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/findclash/findclash.html>`_ based on van der Waals radii. 
* operate: Provides access to Chimera and allows to save the result in Scipion framework.
* restore: This protocol opens Chimera and restores the session previously saved with commands *scipionwrite* or *scipionss*. 


========
Examples
========

See `Model Building Tutorial <https://github.com/I2PC/scipion/wiki/tutorials/tutorial_model_building_basic.pdf>`_


===============
Buildbot status
===============

Status devel version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/chimera_devel.svg

Status production version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/chimera_prod.svg

