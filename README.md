# Chimera plugin

This plugin allows to use chimera commands within the Scipion framework.

Chimera  is a program for interactive visualization and analysis of molecular structures and related data. It is developed by the Resource for Biocomputing, Visualization, and Informatics [(see home page for details)](https://www.cgl.ucsf.edu/chimera/).

Tested with version: 1.13.1

## Installation

You will need to use [2.0](https://github.com/I2PC/scipion/releases/tag/v2.0) version of Scipion to be able to run these protocols. To install the plugin, you have two options:

   a) Stable version
   ```
   scipion installp -p scipion-em-chimera
   ```
   b) Developer's version
   * download repository 
   ```
    git clone https://github.com/scipion-em/scipion-em-chimera.git
   ```
   * install 
   ```
    scipion installp -p path_to_scipion-em-chimera --devel
   ```

Chimera binaries will be installed automatically with the plugin, but you can also link an existing installation. 
Default installation path assumed is `software/em/chimera-1.13.1`, if you want to change it, set *CHIMERA_HOME,* in `scipion.conf` file to the folder where Chimera is installed.

To check the installation, simply run the following Scipion test: 
```
    scipion test chimera.tests.test_protocol_chimera_operate
    scipion test chimera.tests.test_protocol_chimera_fit
    scipion test chimera.tests.test_protocol_modeller_search
    scipion test chimera.tests.test_protocol_contact
```

## Supported versions of Motioncor2

1.13.1

## Protocols

* rigit fit: fits an atomic structure (PDB/CIF) to a 3D map using a rigid fit algorithm.
* model from template: models three-dimensional structures of proteins using [Modeller](https://salilab.org/modeller/manual/node7.html)
* contacts: Identifies interatomic clashes and contacts based on van der Waals radii (https://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/findclash/findclash.html)
* operate: provides access to Chimera and allows to save the result in Scipion framework
* restore: each time a 3Dmap or an atomic structure is saved using `scipionwrite` or `scipionss` commad a chimera session is saved. This protocol opens Chimera and restores the session. 

## Examples
[See Model Building Tutorial](https://github.com/I2PC/scipion/wiki/tutorials/tutorial_model_building_basic.pdf)

## Status Buildbot
Status devel version: ![build status](http://arquimedes.cnb.csic.es:9980/badges/chimera_devel.svg)

Status production version: ![build status](http://arquimedes.cnb.csic.es:9980/badges/chimera_prod.svg)
