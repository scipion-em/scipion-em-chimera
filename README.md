# Chimera plugin

This plugin allows to use chimera commands within the Scipion framework.

Chimera  is a program for interactive visualization and analysis of molecular structures and related data. It is developed by the Resource for Biocomputing, Visualization, and Informatics [(see home page for details)](https://www.cgl.ucsf.edu/chimera/).


Supported version  of Motioncor2: 1.0.2, 1.0.5, 1.1.0, 1.2.1

![build status](http://arquimedes.cnb.csic.es:9980/badges/motioncorr_devel.svg "Build status")

This plugin provide a wrapper around [motioncor2](http://msg.ucsf.edu/em/software/motioncor2.html).

## Installation

You will need to use [2.0](https://github.com/I2PC/scipion/releases/tag/v2.0) version of Scipion to be able to run these protocols. To install the plugin, you have two options:

   a) Stable version
   ```
   scipion installp -p scipion-em-motioncorr
   ```
   b) Developer's version
   * download repository 
   ```
    git clone https://github.com/scipion-em/scipion-em-motioncorr.git
   ```
   * install 
   ```
    scipion installp -p path_to_scipion-em-motioncorr --devel
   ```

Motioncor2 binaries will be installed automatically with the plugin, but you can also link an existing installation. 
Default installation path assumed is `software/em/motioncor2-1.2.1`, if you want to change it, set *MOTIONCOR2_HOME* in `scipion.conf` file to the folder where the Motioncor2 is installed. Depending on your CUDA version you might want to change the default binary from `MotionCor2_1.2.1-Cuda80` to a different one by explicitly setting *MOTIONCOR2* variable. If you need to use CUDA different from the one used during Scipion installation (defined by CUDA_LIB), you can add *MOTIONCOR2_CUDA_LIB* variable to the config file. Various binaries can be downloaded from the official UCSF website.

To check the installation, simply run the following Scipion test: `scipion test motioncorr.tests.test_protocols_motioncor2.TestMotioncor2AlignMovies`

## Supported versions of Motioncor2

1.0.2, 1.0.5, 1.1.0, 1.2.1

## Protocols

* [movie alignment](ProtMotionCorr)

## References
1.  Shawn Q Zheng, Eugene Palovcak, Jean-Paul Armache, Kliment A Verba, Yifan Cheng & David A Agard. MotionCor2: anisotropic correction of beam-induced motion for improved cryo-electron microscopy. Nature Methods volume 14, pages 331–332 (2017).
