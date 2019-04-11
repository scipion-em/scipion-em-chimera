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
* operate: open chimera and 


## References
1.  Shawn Q Zheng, Eugene Palovcak, Jean-Paul Armache, Kliment A Verba, Yifan Cheng & David A Agard. MotionCor2: anisotropic correction of beam-induced motion for improved cryo-electron microscopy. Nature Methods volume 14, pages 331–332 (2017).
![build status](http://arquimedes.cnb.csic.es:9980/badges/motioncorr_devel.svg "Build status")
