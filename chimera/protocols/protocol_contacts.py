from pyworkflow.em.protocol import EMProtocol, Boolean
from pyworkflow.em.constants import (SYM_I222, SYM_I222r, SYM_In25, SYM_In25r,
                                     SYM_CYCLIC, SYM_DIHEDRAL, SYM_TETRAHEDRAL,
                                     SYM_OCTAHEDRAL, SCIPION_SYM_NAME, SYM_I2n3,
                                     SYM_I2n3r, SYM_I2n5, SYM_I2n5r)
from pyworkflow.protocol.params import (EnumParam,
                                        IntParam,
                                        PointerParam,
                                        StringParam,
                                        FloatParam,
                                        LEVEL_ADVANCED, BooleanParam)
import sqlite3
import json
import collections
import os
from pyworkflow.em.viewers.viewer_chimera import Chimera
from chimera import Plugin
from operator import itemgetter

from pyworkflow.utils import red


class ChimeraProtContacts(EMProtocol):
    """Identifies interatomic clashes and contacts based on van der Waals radii
    """
    _label = 'contacts'
    _program = ""
    commandDropView = """DROP view IF EXISTS {viewName}"""
    TetrahedralOrientation = ['222', 'z3']

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.SYMMETRY = Boolean(True)

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
        form.addParam('applySymmetry', BooleanParam,
                         label="Apply symmetry:", default=True,
                         help="'Symmetry = Yes' indicates that symmetry will be applied, and then"
                              " contacts will be computed between any two chains of the "
                              "atomic structure (the unit cell) "
                              "and between a chain of the unit cell and another chain of a "
                              "neigbour unit cell. Output results will show only non "
                              "redundant contatcs, i.e., contacts than you can infer by"
                              " symmetry will not be shown.\n'Symmetry = No' indicates that "
                              "symmetry will not be applied, and then  "
                              "contacts will only be calculated between chains within the "
                              "atomic structure. Output results will show all contacts between"
                              " any couple of interacting chains.\n")
        form.addParam('symmetryGroup', EnumParam,
                      choices=[SCIPION_SYM_NAME[SYM_CYCLIC],
                               SCIPION_SYM_NAME[SYM_DIHEDRAL],
                               SCIPION_SYM_NAME[SYM_TETRAHEDRAL],
                               SCIPION_SYM_NAME[SYM_OCTAHEDRAL],
                               SCIPION_SYM_NAME[SYM_I222],
                               SCIPION_SYM_NAME[SYM_I222r],
                               SCIPION_SYM_NAME[SYM_In25],
                               SCIPION_SYM_NAME[SYM_In25r],
                               SCIPION_SYM_NAME[SYM_I2n3],
                               SCIPION_SYM_NAME[SYM_I2n3r],
                               SCIPION_SYM_NAME[SYM_I2n5],
                               SCIPION_SYM_NAME[SYM_I2n5r]],
                      default=SYM_CYCLIC,
                      label="Symmetry",
                      condition='applySymmetry',
                      help="See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/"
                           "Symmetry for a description of the symmetry groups "
                           "format in Xmipp.\n"
                           "If no symmetry is present, use _c1_."
                      )
        form.addParam('symmetryOrder', IntParam, default=1,
                    condition='applySymmetry and symmetryGroup<=%d' %
                              (SYM_DIHEDRAL),
                    label='Symmetry Order',
                    help='Select the order of cyclic or dihedral symmetry.')
        form.addParam('tetrahedralOrientation', EnumParam,
                      choices=self.TetrahedralOrientation, default=0,
                      condition='applySymmetry and symmetryGroup==%d' %
                                (SYM_TETRAHEDRAL),
                      label='Tetrahedral Orientation',
                      help='Select the orientation of the tetrahedron:\n'
                           '222, by default: Two-fold symmetry axes along '
                           'the X, Y, and Z axes, '
                           'a three-fold along axis (1,1,1).\n'
                           'z3: A three-fold symmetry axis along Z and another three-fold '
                           'axis in the YZ plane.\n'
                           'More information: \n'
                           'https://www.cgl.ucsf.edu/chimera/current/docs/UsersGuide/midas/sym.html')
        # some empty space so if symmetryOrder can be seem
        # without resizing the window
        group = form.addGroup('Fit params for clashes and contacts')
        group.addParam('cutoff', FloatParam,
                       label="cutoff (Angstroms): ", default=-0.4,
                       expertLevel=LEVEL_ADVANCED,
                       help="Large positive cutoff identifies the more severe clashes, "
                            "whereas negative cutoff indicates favorable contacts:\n"
                            "default contact rule: -0.4 (from 0.0 to -1.0)\n"
                            "default clash rule: 0.6 (from 0.4 to 1.0)\n"
                            'More information: \n'
                            'https://www.cgl.ucsf.edu/chimera/current/docs/UsersGuide/midas/sym.html'
                       )
        group.addParam('allowance', FloatParam,
                       label="allowance (Angstroms): ", default=0.0,
                       expertLevel=LEVEL_ADVANCED,
                       help="default contact rule: 0.0\n"
                            "default clash rule: 0.4\n"
                            'More information: \n'
                            'https://www.cgl.ucsf.edu/chimera/current/docs/UsersGuide/midas/sym.html'
                       )
        form.addLine('')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self.sym = Chimera._symmetryMap[self.symmetryGroup.get()]
        self.symOrder = self.symmetryOrder.get()
        if not self.applySymmetry:
            self.sym = "Cn"
            self.symOrder = 1
            self.SYMMETRY = Boolean(False)
        elif ((self.sym == "Cn" or self.sym == "Dn") and self.symOrder == 1):
            self.SYMMETRY = Boolean(False)
        # connect to database, delete table and recreate it
        #execute chimera findclash
        self._insertFunctionStep('chimeraClashesStep')
        self._insertFunctionStep('postProcessStep')

        self._store()

    def postProcessStep(self):
        c, conn = connectDB(self.getDataBaseName(), None)
        self.removeDuplicates(c)

    def chimeraClashesStep(self):
        labelDictAux = json.loads(self.chainStructure.get(),
                               object_pairs_hook=collections.OrderedDict)
        labelDict = collections.OrderedDict(sorted(labelDictAux.items(), key=itemgetter(1)))

        pdbFileName = self.pdbFileToBeRefined.get().getFileName()
        # first element of dictionary
        firstValue = labelDict[list(labelDict)[0]]
        outFiles = []
        f = open(self.getChimeraScriptFileName1(), "w")
        f.write("from chimera import runCommand\n")
        f.write("runCommand('open {}')\n".format(pdbFileName))

        if self.sym == "Cn" and self.symOrder != 1:
            f.write("runCommand('sym #0 group C%d contact 3')\n" % self.symOrder)
        elif self.sym == "Dn" and self.symOrder != 1:
            f.write("runCommand('sym #0 group d%d contact 3')\n" % self.symOrder)
        elif self.sym == "T":
            f.write("runCommand('sym #0 group t,%s contact 3')\n" %
                    self.TetrahedralOrientation[self.tetrahedralOrientation.get()])
            # Look at: https://www.cgl.ucsf.edu/chimera/current/docs/UsersGuide/midas/sym.html
        elif self.sym == "O":
            f.write("runCommand('sym #0 group O contact 3')\n")
        elif self.sym == "222" or self.sym =="222r" or self.sym == "n25" or \
             self.sym =="n25r" or self.sym=="2n3" or self.sym=="2n3r" or \
             self.sym =="2n5" or self.sym=="2n5r":
            f.write("runCommand('sym #0 group i,%s contact 3')\n" % self.sym)
        self.SYMMETRY = self.SYMMETRY.get()
        if self.SYMMETRY:
            f.write("runCommand('write #1 {symmetrizedModelName}')\n".format(
                symmetrizedModelName=self.getSymmetrizedModelName()))

        self.endChimeraScript(firstValue, labelDict, outFiles, f)
        f.close()
        args = " --nogui --script " + self.getChimeraScriptFileName1()
        self._log.info('Launching: ' + Plugin.getProgram() + ' ' + args)
        Chimera.runProgram(Plugin.getProgram(), args)

        if self.SYMMETRY and not os.path.exists(self.getSymmetrizedModelName()):
            # When self.SYMMETRY = TRUE and no one neighbor unit cell has not been
            # generated at less than 3 Angstroms, probably because the symmetry
            # center is not equal to the origin of coordinates, at least we have the
            # contacts that are within the unit cell.
            print (red("Error: No neighbor unit cells are available. "
                       "Is the symmetry center equal to the origin of "
                       "coordinates?"))
            self.SYMMETRY = False
            f = open(self.getChimeraScriptFileName2(), "w")
            f.write("from chimera import runCommand\n")
            f.write("runCommand('open {}')\n".format(pdbFileName))
            self.endChimeraScript(firstValue, labelDict, outFiles, f)
            f.close()
            args = " --nogui --script " + self.getChimeraScriptFileName2()
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
        labelDictAux = json.loads(self.chainStructure.get(),
                               object_pairs_hook=collections.OrderedDict)
        labelDict = collections.OrderedDict(sorted(labelDictAux.items(), key=itemgetter(1)))
        d = {}
        d1 = {}
        d2 = {}
        anyResult = False
        for inFile in outFiles:
            print "processing file", inFile
            if not os.path.exists(inFile):
                continue
            else:
                anyResult = True
            counter = 0
            # parse contact files. Note that C1 symmetry file is different from the rest
            for line in open(inFile):
                if counter < 8:
                    # print "skip line", line
                    counter += 1
                else:
                    if not self.SYMMETRY:
                        info = line.split() # ['HIS', '87.A', 'NE2', 'HEM', '1.A002', 'ND', '0.620', '2.660']
                        d1['modelId'] = "'" + "#0" + "'"
                        d1['aaName'] = "'" + info[0][0] + info[0][1:].lower() + "'"  # "'HIS'"
                        info2 = info[1].split(".")  # ['87', 'A']
                        d1['aaNumber'] = info2[0]  # 87
                        d1['chainId'] = "'" + info2[1] + "'"  # A
                        d1['atomId'] = "'" + info[2] + "'"  # NE2
                        d1['protId'] = "'" + labelDict[info2[1]] + "'"

                        d2['modelId'] = "'" + "#0" + "'"
                        d2['aaName'] = "'" + info[3][0] + info[3][1:].lower() + "'"  # HEM
                        info2 = info[4].split(".")
                        d2['aaNumber'] = info2[0]  # 1
                        d2['chainId'] = "'" + info2[1] + "'"  # A002
                        d2['protId'] = "'" + labelDict[info2[1]] + "'"
                        d2['atomId'] = "'" + info[5] + "'"  # ND
                        d['overlap'] = info[6]  # 0.620
                        d['distance'] = info[7]  # 2.660

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

        return anyResult

    #    --------- util functions -----

    def getDataBaseName(self):
        return self._getExtraPath("overlaps.sqlite")

    def getSymmetrizedModelName(self):
        return os.path.abspath(self._getExtraPath("symModel.pdb"))

    def getTableName(self):
        return "contacts"

    def getChimeraScriptFileName1(self):
        return self._getTmpPath("chimera1.cmd")

    def getChimeraScriptFileName2(self):
        return self._getTmpPath("chimera2.cmd")

    def endChimeraScript(self, firstValue, labelDict, outFiles, f):
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
                    """runCommand('echo {}')\nrunCommand('findclash  #0:{} test other savefile {} overlap {} hbond {} namingStyle simple')\n""".format(
                        chains, chains, outFile, self.cutoff, self.allowance))
                protId = v
                chains = ".{}".format(k)
                outFileBase = v
        outFile = os.path.abspath(self._getExtraPath("{}.over".format(outFileBase)))
        outFiles.append(outFile)

        f.write(
            """runCommand('echo {}')\nrunCommand('findclash  #0:{} test other savefile {} overlap {} hbond {} namingStyle simple')\n""".format(
                chains, chains, outFile, self.cutoff, self.allowance))
        # f.write("runCommand('save %s')\n" % os.path.abspath(self._getExtraPath(sessionFile)))


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

    def _validate(self):
        errors = []
        if (self.symmetryOrder.get() <= 0):
            errors.append("Error: Symmetry Order should be a positive integer" )

        return errors

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
  

