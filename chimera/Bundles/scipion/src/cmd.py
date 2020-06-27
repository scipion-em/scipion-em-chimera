# devel install /home/roberto/Software/Plugins3/scipion-em-chimera/chimera/Bundles/scipion
import configparser
from chimerax.core.commands import CmdDesc      # Command description
from chimerax.core.commands import TopModelsArg
from chimerax.core.commands import StringArg
from chimerax.core.commands import run # execute chimera cli command
from chimerax.map.volume import Volume # model type 3D map
from chimerax.atomic.structure import AtomicStructure # model type atomic structure

import os
import ntpath
chimeraConfigFileName = "chimera.ini"

def readConfigFile(session, configFileName):
    config = configparser.ConfigParser()
    config.read(configFileName)
    d = {}
    if not config.has_section('chimerax'):
        session.logger.info("No sectipon chimerax")
        d['enablebundle'] = False
    else:
        session.logger.info("Yes sectipon chimera")
        d['chimerapdbtemplatefilename'] = config.get("chimerax", "chimeraPdbTemplateFileName")
        d['chimeramaptemplatefilename'] = config.get("chimerax", "chimeraMapTemplateFileName")
        d['sessionfile'] = config.get("chimerax", "sessionFile")
        d['protid'] = config.getint("chimerax", "protId")
        d['enablebundle'] = config.getboolean("chimerax", "enablebundle")
    session.logger.info("d=%s" % str(d))
    return d

def scipionwrite(session, model, prefix=None):
    # models is a tuple with all selected models but we are only
    # process the first one

    d = readConfigFile(session, chimeraConfigFileName)
    if not d['enablebundle']:
        session.logger.error("scipionwriter cannot be called from Analyze or Viewers")
        return

    model = model[0] # model is a tuple, let us get the major model id

    # remoce "." from prefix string
    if prefix is not None:
        prefix = prefix.replace(".", "_dot_")

    # session.logger.info("modelID %s" % str(dir(model)))
    if isinstance(model, AtomicStructure):
        modelFileName = d['chimerapdbtemplatefilename'].replace("__","__%s_" % (model.id)[0] )
    elif isinstance(model, Volume):
        modelFileName = d['chimeramaptemplatefilename'].replace("__","__%s_" % (model.id)[0] )
    else:
        session.logger.error("I do not know how to save model %s\n" %
                             model._get_name())
        return

    if prefix is not None:
        modelFileName = os.path.join(ntpath.dirname(modelFileName),
                                   prefix + ntpath.basename(modelFileName))
    command = 'save %s #%s'%(modelFileName, str((model.id)[0]))
    session.logger.info(command)
    run(session, command)

    session.logger.info("Saving session")
    command = 'save %s' % d['sessionfile']
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
    d = readConfigFile(session, chimeraConfigFileName)
    if not d['enablebundle']:
        session.logger.error("scipionwriter cannot be called from Analyze or Viewers")
        return

    command = 'save %s' % d['sessionfile']
    run(session, command)

scipionss_desc = CmdDesc()

def scipionrs(session):
    # All command functions are invoked with ``session`` as its
    # first argument.  Useful session attributes include:
    #   logger: chimerax.core.logger.Logger instance
    #   models: chimerax.core.models.Models instance
    session.logger.info("Restoring session")
    d = readConfigFile(session, chimeraConfigFileName)
    if not d['enablebundle']:
        session.logger.error("scipionwriter cannot be called from Analyze or Viewers")
        return

    command = 'open %s' % d['sessionfile']
    run(session, command)

scipionrs_desc = CmdDesc()
