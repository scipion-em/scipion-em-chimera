# When a test is executed the default working directory is tmp even if extra has been
# passed as ChimeraX argument. The following hack bypasses that problem.
from .constants import *

# devel install /home/roberto/Software/Plugins3/scipion-em-chimera/chimera/Bundles/scipion
from chimerax.core.commands import CmdDesc      # Command description
from chimerax.core.commands import TopModelsArg, ModelsArg
from chimerax.core.commands import StringArg, FloatArg, IntArg

from chimerax.core.commands import run # execute chimera cli command
from chimerax.map.volume import Volume # model type 3D map
from chimerax.atomic.structure import AtomicStructure # model type atomic structure

import os, time
# import numpy as np
import ntpath


def getConfig(session, key):
        return os.environ.get(key, False)


def scipionshellcrown(session,
                 model=None,
                 sphereRadius=None,
                 orientation='222r',
                 sphereFactor='1',
                 modelid=None,
                 crownwidth=10):

    """ Shows a shell of a 3D Map at radius sphereRadius
    and with width=width

    :param session:  chimera session
    :param model: 3D map to be analyzed
    :param sphereRadius: 3D density obtained at this radius
    :param orientation: virus symmetry, i.e: 222, 222r, etc.
    :param sphereFactor: allows generating a shape that is an
    interpolation between an icosahedron and a sphere of equal radius.
    :param modelid: output model_id. keep always the same value to
     overwrite the previous surface.
    :param crownwidth: shell width
    :return: none
    example: scipionshellcrown #1 220  10  orientation 222r   sphereFactor 1 crownwidth 60

    """
    # set model style as surface
    mapModelId = model[0].id_string # model is a tuple, id = 1
    command = "volume #%s style surface " % (mapModelId)
    run(session, command)

    # very likely divisions is irrelevant here
    command = "shape icosahedron " \
              "radius %s " \
              "divisions 2000 " \
              "orientation %s " \
              "sphereFactor %s" % (sphereRadius,
                                   orientation,
                                   sphereFactor)
    outIcosahedronId = run(session, command).id_string # id = 2

    # delete output model
    command = "close #%s;" % modelid
    run(session, command)

    # mask  in icosahedron
    command = "volume mask #%s surfaces #%s " \
              "fullMap true modelId %s slab %s" % \
              (mapModelId, outIcosahedronId, modelid, crownwidth) # model id
    run(session, command)

    command = "color #%s gray all; " \
              "lighting soft; " \
              "cofr coordinateSystem #%s; " \
              "lighting shadows true; " \
              'ui tool show "Side View"' \
              % (modelid, modelid)
    run(session, command)

    command = 'hide #%s;' \
              'close #%s' \
              %(mapModelId,
                outIcosahedronId)
    run(session, command)

scipionshellcrown_desc = CmdDesc(
    required=[('model', TopModelsArg),
              ('sphereRadius', StringArg),
              ('modelid', StringArg)
              ],
    optional= [
               ('orientation', StringArg),
               ('sphereFactor', StringArg),
               ('crownwidth', StringArg),
               ]
)


def scipionshell(session,
                 model=None,
                 sphereRadius=None,
                 orientation='222r',
                 sphereFactor='1',
                 modelid=None):

    """ Shows the density of a 3D Map on a spherical shell
    of the map at radius sphereRadius

    :param session:  chimera session
    :param model: 3D map to be analized
    :param sphereRadius: 3D density optained at this radius
    :param orientation: virus symmetry, i.e: 222, 222r, etc.
    :param sphereFactor: allows generating a shape that is an
    interpolation between an icosahedron and a sphere of equal radius.
    :param modelid: output model_id. keep always the same value to
     overwrite the previous surface.
    :return: none
    """

    mapModelId = model[0].id_string # model is a tuple
    [(x0, y0, z0), (x1, y1, z1)] = model[0].ijk_bounds()
    command = "volume #%s style image " \
              "positionPlanes %d,%d,%d " \
              "orthoplanes xyz;" \
              "ui mousemode right 'move planes' " % (mapModelId,
                                   x1//2,  y1//2, z1//2
                                   )
    run(session, command)
    command = "shape icosahedron " \
              "mesh false " \
              "divisions 2000 " \
              "radius %s " \
              "orientation %s " \
              "modelid %s " \
              "sphereFactor %s" % (sphereRadius,
                                   orientation,
                                   modelid,
                                   sphereFactor)
    icosahedronId = run(session, command).id_string
    session.logger.info("icosahedronId: " + str(icosahedronId))
    command = "color sample #%s map #%s palette gray; " \
              "lighting soft; cofr coordinateSystem #%s" \
              %(icosahedronId, mapModelId, mapModelId)
    run(session, command)

scipionshell_desc = CmdDesc(
    required=[('model', TopModelsArg),
              ('sphereRadius', StringArg),
              ('modelid', StringArg)
              ],
    optional= [
               ('orientation', StringArg),
               ('sphereFactor', StringArg),
               ]
)


def scipioncombine(session, models=None, modelid=None):
    """Emulate copy/combine chimera command. It will not handle
    properly models with submodels but it is better than nothing.
    scipioncombine #1,2 modelid 55
    will combine model 1 and 2 and produce an output model with
    id = 55
    arg modelid is optional
    """

    if not checkBundleEnabled(session, cmd="scipioncombine"):
        return

    # list with all models
    modelFileName = []

    for model in models:
        modelName = model._get_name()
        if isinstance(model, AtomicStructure) or\
            modelName[-4:] == '.cif' or\
            modelName[-4:] == '.pdb':
            # files starting with "tmp_" will not be converted in scipion objects
            modelFileName.append(getConfig(session, CHIMERA_PDB_TEMPLATE_FILE_NAME).
                                 replace("__","_in_%s_" % (model.id)[0]).
                                 replace("Atom_struct_", "tmp_Atom_struct_"))
        else:
            session.logger.error("I do not know how to combine model %s\n" % modelName)
            continue
        # save each model to be combined
        command = 'save %s #%s'%(modelFileName[-1],
                                 str((model.id)[0]))
        run(session, command)

    # create add script for scipion atomatructurils
    outFileName = getConfig(session, CHIMERA_PDB_TEMPLATE_FILE_NAME).\
        replace("__", "_out_%s_" % (model.id)[0])
    scriptFileName = outFileName[:-4] + ".py"

    f = open(scriptFileName, "w")
    f.write("""from pwem.convert.atom_struct import AtomicStructHandler
# recover list with all models to be combined    
modelFileName = eval("%s")
# read first model
aStruct1 = AtomicStructHandler(modelFileName[0])
# read rest of models but last one
for fileName in modelFileName[1:-1]:
    aStruct1.addStruct(fileName)
# read last model and save the sum of all
aStruct1.addStruct(modelFileName[-1], '%s')
"""% (str(modelFileName), outFileName))
    f.close()
    cmd = getConfig(session, SCIPIONPYTHON)
    cmd += " %s" % scriptFileName
    os.system(cmd)
    newModel = run(session, "open %s" % outFileName )
    # TODO: we do not fully understand how
    # models work in chimerax. The follwong two
    # way of getting modelID should be equivalent
    # but sometimes one of them fail.
    try:
        modelID = str(newModel[0][0].id[0])
    except:
        modelID = str(session.models[-1].id[0])
    session.logger.info("model --> " + modelID)
    run(session, "style #%s stick" % modelID )
    if modelid is not None:
        run(session, "rename  #%s id #%s" % (modelID, modelid))

scipioncombine_desc = CmdDesc(
    required = [('models', TopModelsArg)],
    optional= [('modelid', StringArg)],
)


scipion_desc = CmdDesc()

def checkBundleEnabled(session, cmd="scipion"):
    """ Checks if bundle is enabled, this should happens under protocol execution process
    and not during visualization"""
    if not getConfig(session, PROTID):
        session.logger.error("%s command can not be executed in this context" % cmd)
        return False
    else:
        return True


def scipionwrite(session, model, prefix=None):
    # models is a tuple with all selected models but we are only
    # process the first one

    if not checkBundleEnabled(session, cmd="scipionwrite"):
        return

    model = model[0] # model is a tuple, let us get the major model id

    # remove "." from prefix string
    if prefix is not None:
        prefix = prefix.replace(".", "_dot_")

    # session.logger.info("modelID %s" % str(dir(model)))
    modelName = model._get_name()
    if isinstance(model, AtomicStructure) or \
        (not isinstance(model, Volume) and modelName.find('cif')!= -1):
        modelFileName = getConfig(session, CHIMERA_PDB_TEMPLATE_FILE_NAME).replace("__", "__%s_" % (model.id)[0])
    elif (isinstance(model, Volume) or modelName.find('mrc')!= -1):
        modelFileName = getConfig(session, CHIMERA_MAP_TEMPLATE_FILE_NAME).replace("__", "__%s_" % (model.id)[0])
    else:
        session.logger.error("I do not know how to save model %s\n" % modelName)
        return

    if prefix is not None:
        modelFileName = os.path.join(ntpath.dirname(modelFileName),
                                   prefix + ntpath.basename(modelFileName))
    command = 'save %s #%s'%(modelFileName, str((model.id)[0]))
    session.logger.info(command)
    run(session, command)

    if not (prefix == "DONOTSAVESESSION_"):
        session.logger.info("Saving session")
        command = 'save %s' % getConfig(session, SESSIONFILE)
        run(session, command)

scipionwrite_desc = CmdDesc(
    required = [('model', TopModelsArg)],
    optional= [('prefix', StringArg)],
)

def scipionss(session):
    # All command functions are invoked with ``session`` as its
    # first argument.  Useful session attributes include:
    #   logger: chimerax.core.logger.Logger instance
    #   models: chimerax.core.models.Models instance
    session.logger.info("Saving session")

    if not checkBundleEnabled(session, cmd="scipionss"):
        return

    command = 'save %s' % getConfig(session, SESSIONFILE)
    run(session, command)

scipionss_desc = CmdDesc()

def scipionrs(session):
    # All command functions are invoked with ``session`` as its
    # first argument.  Useful session attributes include:
    #   logger: chimerax.core.logger.Logger instance
    #   models: chimerax.core.models.Models instance
    session.logger.info("Restoring session")
    if not checkBundleEnabled(session, cmd="scipionrs"):
        return
    command = 'open %s' % getConfig(SESSIONFILE)
    run(session, command)

scipionrs_desc = CmdDesc()

def scipion(session):
    print("""List of Scipion commands:
    scipionwrite: saves model file
    scipionss: saves chimera session
    scipionrs: restores chimera sesion
    scipioncombine: combines two models
    scipionshell: shows the density of a 3D Map on a spherical shell of the map at a given radius
    scipionshellcrown: cut a shell from a 3D Map at a given radius
    
    type "help command_name" for more information
""")