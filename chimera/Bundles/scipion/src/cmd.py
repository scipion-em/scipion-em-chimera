# When a test is executed the default working directory is tmp even if extra has been
# passed as ChimeraX argument. The following hack bypasses that problem.
from .constants import *

# devel install /home/roberto/Software/Plugins3/scipion-em-chimera/chimera/Bundles/scipion
from chimerax.core.commands import CmdDesc      # Command description
from chimerax.core.commands import TopModelsArg, ModelsArg
from chimerax.core.commands import StringArg
from chimerax.core.commands import run # execute chimera cli command
from chimerax.map.volume import Volume # model type 3D map
from chimerax.atomic.structure import AtomicStructure # model type atomic structure

import os
import ntpath


def getConfig(session, key):
        return os.environ.get(key, False)

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
