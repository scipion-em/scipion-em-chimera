from pyworkflow.em.protocol import EMProtocol
from pyworkflow.em.constants import (SYM_I222, SYM_I222r, SYM_In25, SYM_In25r,
                                     SYM_CYCLIC, SYM_DIHEDRAL, SYM_TETRAHEDRAL,
                                     SYM_OCTAHEDRAL, SCIPION_SYM_NAME)
from pyworkflow.protocol.params import (EnumParam,
                                        IntParam,
                                        PointerParam,
                                        StringParam,
                                        LEVEL_ADVANCED)
import sqlite3
import json
import collections
import os
from pyworkflow.em.viewers.viewer_chimera import Chimera
from chimera import Plugin
from pyworkflow.em.viewers.viewer_chimera import (sessionFile)

class ChimeraProtContacts(EMProtocol):
    _label = 'contacts'
    _program = ""
    commandDropView = """DROP view IF EXISTS {viewName}"""

    def _defineParams(self, form):
        form.addSection(label='Input')
        #pdbFileToBeRefined name is needed by the wizard. Do not change it
        form.addParam('pdbFileToBeRefined', PointerParam, pointerClass="AtomStruct",
                       label='Atomic Structure:', allowsNull=True,
                       important=True,
                       help="Input atomic structure.")
        form.addParam('chainStructure', StringParam, default="",
                      label='Chain Labeling',
                      help="Dictionary that maps chains to labels.\n"
                           "Example: {'A':'h1', 'B':'h1', 'E':'h2'}\n"
                           "Contacts are calculated between two chains with distinct "
                           "labels. Two chains with the same label are considered as "
                           "a group. Contacts will be computed between any chain included "
                           "in this group and any other group/chain. However, no contacts "
                           "among members of the group will be calculated.")
        form.addParam('symmetryGroup', EnumParam,
                      choices=[SCIPION_SYM_NAME[SYM_CYCLIC],
                               SCIPION_SYM_NAME[SYM_DIHEDRAL],
                               SCIPION_SYM_NAME[SYM_TETRAHEDRAL],
                               SCIPION_SYM_NAME[SYM_OCTAHEDRAL],
                               SCIPION_SYM_NAME[SYM_I222],
                               SCIPION_SYM_NAME[SYM_I222r],
                               SCIPION_SYM_NAME[SYM_In25],
                               SCIPION_SYM_NAME[SYM_In25r] ],
                      default=SYM_I222r,
                      label="Symmetry",
                      help="See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/"
                           "Symmetry for a description of the symmetry groups "
                           "format in Xmipp.\n"
                           "If no symmetry is present, use _c1_."
                      )
        form.addParam('symmetryOrder', IntParam, default=1,
                    condition='symmetryGroup<=%d' % SYM_DIHEDRAL,
        label='Symmetry Order',
        help='Select the order of cyclic or dihedral symmetry.')
        # some empty space so if symmetryOrder can be seem
        # without resizing the window
        form.addLine('')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        # connect to database, delete table and recreate it
        #execute chimera findclash
        self._insertFunctionStep('chimeraClashesStep')
        self._insertFunctionStep('postProcessStep')
        #self._insertFunctionStep('',
        #                         self.pdbFileToBeRefined.get().getFileName(),
        #                         self.chainStructure
        #                         )

    def postProcessStep(self):
        c, conn = connectDB(self.getDataBaseName(), None)
        self.removeDuplicates(c)

    def chimeraClashesStep(self):
        labelDict = json.loads(self.chainStructure.get(),
                               object_pairs_hook=collections.OrderedDict)
        pdbFileName = self.pdbFileToBeRefined.get().getFileName()
        sym = Chimera._symmetryMap[self.symmetryGroup.get()]
        # first elemet of dictionary
        firstValue = labelDict[list(labelDict)[0]]
        outFiles = []
        f = open(self.getChimeraScriptFileName(), "w")
        f.write("from chimera import runCommand\n")
        f.write("runCommand('open {}')\n".format(pdbFileName))
        f.write("runCommand('sym #0 group i,222r contact 3')\n")  # apply symmetry

        protId = firstValue
        chains = ""
        comma = ''
        for k, v in labelDict.iteritems():
            if protId == v:
                chains += "{}.{}".format(comma, k)
                comma = ','
                outFileBase = v
            else:
                outFile = os.path.abspath(self._getExtraPath("{}.over".format(outFileBase)))
                outFiles.append(outFile)
                f.write(
                    """runCommand('echo {}')\nrunCommand('findclash  #0:{} test other savefile {} overlap -0.4 hbond 0.0 namingStyle simple')\n""".format(
                        chains, chains, outFile))
                protId = v
                chains = ".{}".format(k)
                outFileBase = v
        outFile = os.path.abspath(self._getExtraPath("{}.over".format(outFileBase)))
        outFiles.append(outFile)
        f.write(
            """runCommand('echo {}')\nrunCommand('findclash  #0:{} test other savefile {} overlap -0.4 hbond 0.0 namingStyle simple')\n""".format(
                chains, chains, outFile))
        #f.write("runCommand('save %s')\n" % os.path.abspath(self._getExtraPath(sessionFile)))
        f.close()
        args = " --nogui --script " + self.getChimeraScriptFileName()
        self._log.info('Launching: ' + Plugin.getProgram() + ' ' + args)
        Chimera.runProgram(Plugin.getProgram(), args)

        # parse all files created by chimera
        c, conn = self.prepareDataBase()
        self.parseFiles(outFiles, c)
        conn.commit()
        conn.close()
        #return outFiles

    def prepareDataBase(self, drop=True):
        if drop:
            return  connectDB(self.getDataBaseName(), self.getTableName())
        else:
            return  connectDB(self.getDataBaseName())


    def parseFiles(self, outFiles, c):
        labelDict = json.loads(self.chainStructure.get(),
                               object_pairs_hook=collections.OrderedDict)
        d = {}
        d1 = {}
        d2 = {}
        for inFile in outFiles:
            print "processing file", inFile
            counter = 0
            for line in open(inFile):
                if counter < 8:
                    # print "skip line", line
                    counter += 1
                else:
                    info = line.split()  # ['#2', 'TYR', '851.B', 'CE1', '#0.7', 'PRO', '78.N', 'CA', '3.059', '0.581']
                    d1['modelId'] = "'" + info[0] + "'"  # '#2'
                    d1['aaName'] = "'" + info[1][0] + info[1][1:].lower() + "'"  # "'TYR'"
                    info2 = info[2].split(".")  # ['851', 'B']
                    d1['aaNumber'] = info2[0]  # 851
                    d1['chainId'] = "'" + info2[1] + "'"  # B
                    d1['atomId'] = "'" + info[3] + "'"  # CE1
                    d1['protId'] = "'" + labelDict[info2[1]] + "'"

                    d2['modelId'] = "'" + info[4] + "'"  # '#0.7'
                    d2['aaName'] = "'" + info[5][0] + info[5][1:].lower() + "'"  # PRO
                    info2 = info[6].split(".")
                    d2['aaNumber'] = info2[0]  # 78
                    d2['chainId'] = "'" + info2[1] + "'"  # N
                    d2['protId'] = "'" + labelDict[info2[1]] + "'"
                    d2['atomId'] = "'" + info[7] + "'"  # CA

                    d['overlap'] = info[8]  # 3.059
                    d['distance'] = info[9]  # 0.5

                    if d1['modelId'] == d2['modelId']:
                        if d1['protId'] <= d2['protId']:
                            for k in d1.keys():
                                d[k + '_1'] = d1[k]
                                d[k + '_2'] = d2[k]
                        else:
                            for k in d1.keys():
                                d[k + '_1'] = d2[k]
                                d[k + '_2'] = d1[k]
                    else:
                        if d1['modelId'] <= d2['modelId']:
                            for k in d1.keys():
                                d[k + '_1'] = d1[k]
                                d[k + '_2'] = d2[k]
                        else:
                            for k in d1.keys():
                                d[k + '_1'] = d2[k]
                                d[k + '_2'] = d1[k]

                    command = "INSERT INTO contacts "
                    keys = "("
                    values = " ("
                    for key, value in d.iteritems():
                        keys += key + ", "
                        values += str(value) + ", "
                    keys = keys[:-2] + ")"
                    values = values[:-2] + ")"

                    command += keys + " VALUES " + values
                    ##print command
                    c.execute(command)

    #    --------- util functions -----

    def getDataBaseName(self):
        return self._getExtraPath("overlaps.sqlite")

    def getTableName(self):
        return "contacts"

    def getChimeraScriptFileName(self):
        return self._getTmpPath("chimera.cmd")

    def removeDuplicates(self, c):
        # Remove duplicate contacts
        # that is, given chains A,B
        # we have contact A-B and B-A
        commandEliminateDuplicates = """CREATE VIEW {} AS
        SELECT DISTINCT modelId_1,
             protId_1,
             chainId_1,
             aaName_1,
             aaNumber_1,
             atomId_1,
             modelId_2,
             protId_2,
             chainId_2,
             aaName_2,
             aaNumber_2,
             atomId_2,
             overlap,
             distance
        FROM {}

        """
        commandEliminateDuplicates2 = """
        CREATE VIEW {} AS
        SELECT *
        FROM {}

        EXCEPT

        SELECT ca.*
        FROM {} ca, {} cb
        WHERE
              ca.protId_1    = cb.protId_2
          AND cb.protId_1    = ca.protId_2
          AND ca.chainId_1   = cb.chainId_2
          AND cb.chainId_1   = ca.chainId_2
          AND ca.aaNumber_1  = cb.aaNumber_2
          AND cb.aaNumber_1  = ca.aaNumber_2
          AND ca.atomId_1  = cb.atomId_2
          AND cb.atomId_1  = ca.atomId_2
          AND ca.modelId_2   > cb.modelId_2

        """
        # # Remove duplicate contacts
        # that is, given chains A,B
        # we have contact A.a-B.b and B.b-A.a
        c.execute(self.commandDropView.format(viewName="view_ND_1"))
        c.execute(commandEliminateDuplicates.format("view_ND_1",
                                                        "contacts",
                                                        "contacts"))

        # remove duplicate contacts due to symmetry
        # h1-h1p, h1-h2p
        c.execute(self.commandDropView.format(viewName="view_ND_2"))
        c.execute(commandEliminateDuplicates2.format("view_ND_2", "view_ND_1", "view_ND_1", "view_ND_1"))


def connectDB(sqliteFN, tableName=None):
    conn = sqlite3.connect(sqliteFN)
    c = conn.cursor()
    if tableName is not None:
        commandDropTable = """DROP TABLE IF EXISTS {}"""
        commandCreateTable = """
        CREATE TABLE {}(
             id integer primary key autoincrement,
             modelId_1  char(8),
             protId_1   char(8),
             chainId_1  char(8),
             aaName_1   char(3),
             aaNumber_1 int,
             atomId_1   char(8),
             modelId_2  char(8),
             protId_2   char(8),
             chainId_2  char(8),
             aaName_2   char(3),
             aaNumber_2 int,
             atomId_2   char(8),
             overlap float,
             distance float
             );"""

        c.execute(commandDropTable.format(tableName))
        c.execute(commandCreateTable.format(tableName))
    return c, conn
