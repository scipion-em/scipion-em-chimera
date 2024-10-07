
# 1) scipionwrite command returns error: "Unkown command: scipionwrite"

    scipionwrite is a chimerax plugin. It is installed by the 
    scipion-em-chimera setup. If it does not work check that it 
    has been installed correctly, just open chimerax and execute in the command line

    `help scipionwrite`
    
    If installed, a help page will appear. If it is not installed
    type in the  command line:

    devel install /path_to_scipion3_plugins/scipion-em-chimera/chimera/Bundles/scipion
    
    where path_to_scipion3_plugins is the path to the directory with the Scipion3 plugins.
    
 The command

      `chimera --nogui --cmd "devel install /path_to_scipion3_plugins/scipion-em-chimera/chimera/Bundles/scipion; exit"`

may be useful the determine the cause of the error

# 2) Can not install scipion bundle because help documentation cannot be created

copy /path_to_scipion3_plugins/scipion-em-chimera/chimera/Bundles/scipion to a directory in which you have write/read permissions and install the Bundle from there

