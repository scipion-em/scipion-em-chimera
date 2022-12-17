# **************************************************************************
# *
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *
# * L'Institut de genetique et de biologie moleculaire et cellulaire (IGBMC)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from pwem.protocols import EMProtocol
from pyworkflow.object import Boolean
from pwem.constants import (SYM_DIHEDRAL_X)
from ..convert import CHIMERA_LIST
from ..constants import (CHIMERA_SYM_NAME, CHIMERA_I222)

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
from pwem.viewers.viewer_chimera import Chimera
from .. import Plugin
from operator import itemgetter

from pyworkflow.utils import red


class ChimeraProtContacts(EMProtocol):
    """Identifies interatomic clashes and contacts based on van der Waals radii
    """
    _label = 'contacts'
    _program = ""
    commandDropView = """DROP view IF EXISTS {viewName}"""
    TetrahedralOrientation = ['222', 'z3']

    @classmethod
    def getClassPackageName(cls):
        return "chimerax"

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.SYMMETRY = Boolean(True)

    def _defineParams(self, form):
        form.addSection(label='Input')
        # pdbFileToBeRefined name is needed by the wizard. Do not change it
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
                      choices=CHIMERA_LIST,
                      default=CHIMERA_I222,
                      label="Symmetry",
                      condition='applySymmetry',
                      help="https://scipion-em.github.io/docs/release-2.0.0/docs/developer/symmetries.html?highlight=symmetry"
                           "Symmetry for a description of the symmetry groups "
                           "format in CHIMERA.\n"
                           "If no symmetry is present, use _c1_."
                           'More information: \n'
                           'https://www.cgl.ucsf.edu/chimera/current/docs/UsersGuide/midas/sym.html'
                      )
        form.addParam('symmetryOrder', IntParam, default=1,
                      condition='applySymmetry and symmetryGroup<=%d' % SYM_DIHEDRAL_X,
                      label='Symmetry Order',
                      help='Select the order of cyclic or dihedral symmetry.')

        group = form.addGroup('Fit params for clashes and contacts')
        group.addParam('cutoff', FloatParam,
                       label="cutoff (Angstroms): ", default=-0.4,
                       expertLevel=LEVEL_ADVANCED,
                       help="Large positive cutoff identifies the more severe clashes, "
                            "whereas negative cutoff indicates favorable contacts:\n"
                            "default contact rule: -0.4 (from 0.0 to -1.0)\n"
                            "default clash rule: 0.6 (from 0.4 to 1.0)\n"
                            'More information: \n'
                            'https://www.cgl.ucsf.edu/chimerax/docs/user/commands/clashes.html#top'
                       )
        group.addParam('allowance', FloatParam,
                       label="allowance (Angstroms): ", default=0.0,
                       expertLevel=LEVEL_ADVANCED,
                       help="default contact rule: 0.0\n"
                            "default clash rule: 0.4\n"
                            'More information: \n'
                            'https://www.cgl.ucsf.edu/chimerax/docs/user/commands/clashes.html#top'
                       )
        form.addLine('')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self.sym = CHIMERA_SYM_NAME[self.symmetryGroup.get()]
        self.symOrder = self.symmetryOrder.get()
        if not self.applySymmetry:
            self.sym = "Cn"
            self.symOrder = 1
            self.SYMMETRY = Boolean(False)
        elif (self.sym == "Cn" or self.sym == "Dn") and self.symOrder == 1:
            self.SYMMETRY = Boolean(False)
        # connect to database, delete table and recreate it
        # execute chimera findclash
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
        # labelDict = collections.OrderedDict(sorted(list(labelDictAux.items()), key=itemgetter(1)))
        pdbFileName = os.path.abspath(self.pdbFileToBeRefined.get().getFileName())
        # first element of dictionary
        firstValue = labelDict[list(labelDict)[0]]
        outFiles = []
        f = open(self.getChimeraScriptFileName1(), "w")
        f.write("from chimerax.core.commands import run\n")
        f.write("run(session, 'open {}')\n".format(pdbFileName))
        if self.sym == "Cn" and self.symOrder != 1:
            f.write("run(session,'sym #1 C%d copies t')\n" % self.symOrder)
        elif self.sym == "Dn" and self.symOrder != 1:
            f.write("run(session,'sym #1 d%d copies t')\n" % self.symOrder)
        elif self.sym == "T222" or self.sym == "TZ3":
            f.write("run(session,'sym #1 t,%s copies t')\n" % self.sym[1:])
        elif self.sym == "O":
            f.write("run(session,'sym #1 O copies t')\n")
        elif self.sym == "I222" or self.sym == "I222r" or self.sym == "In25" or \
                self.sym == "In25r" or self.sym == "I2n3" or self.sym == "I2n3r" or \
                self.sym == "I2n5" or self.sym == "I2n5r":
            f.write("run(session,'sym #1 i,%s copies t')\n" % self.sym[1:])
        self.SYMMETRY = self.SYMMETRY.get()
        if self.SYMMETRY:
            f.write("run(session,'delete #2 & #1 #>3')\n")
            f.write("run(session,'save {symmetrizedModelName} #2')\n".format(
                symmetrizedModelName=self.getSymmetrizedModelName()))
            f.write("run(session, 'close #1')\n")
            f.write("run(session, 'rename #2 id #1')\n")
        self.endChimeraScript(firstValue, labelDict, outFiles, f)
        f.write("run(session, 'exit')\n")
        f.close()
        args = " --nogui --script " + self.getChimeraScriptFileName1()
        self._log.info('Launching: ' + Plugin.getProgram() + ' ' + args)
        Chimera.runProgram(Plugin.getProgram(), args)

        if self.SYMMETRY and not os.path.exists(self.getSymmetrizedModelName()):
            # When self.SYMMETRY = TRUE and no one neighbor unit cell has not been
            # generated at less than 3 Angstroms, probably because the symmetry
            # center is not equal to the origin of coordinates, at least we have the
            # contacts that are within the unit cell.
            print(red("Error: No neighbor unit cells are available. "
                      "Is the symmetry center equal to the origin of "
                      "coordinates?"))
            self.SYMMETRY = False
            f = open(self.getChimeraScriptFileName2(), "w")
            f.write("from chimerax.core.commands import run\n")
            f.write("session, run('open {}')\n".format(pdbFileName))
            self.endChimeraScript(firstValue, labelDict, outFiles, f)
            f.write("run(session, 'exit')\n")
            f.close()
            args = " --nogui --script " + self.getChimeraScriptFileName2()
            self._log.info('Launching: ' + Plugin.getProgram() + ' ' + args)
            Chimera.runProgram(Plugin.getProgram(), args)

        # parse all files created by chimera
        c, conn = self.prepareDataBase()
        self.parseFiles(outFiles, c)
        conn.commit()
        conn.close()
        # return outFiles

    def prepareDataBase(self, drop=True):
        if drop:
            return connectDB(self.getDataBaseName(), self.getTableName())
        else:
            return connectDB(self.getDataBaseName())

    def parseFiles(self, outFiles, c):
        labelDictAux = json.loads(self.chainStructure.get(),
                                  object_pairs_hook=collections.OrderedDict)
        labelDict = collections.OrderedDict(sorted(labelDictAux.items(), key=itemgetter(1)))
        # labelDict = collections.OrderedDict(sorted(list(labelDictAux.items()), key=itemgetter(1)))
        d = {}
        d1 = {}
        d2 = {}
        anyResult = False
        for inFile in outFiles:
            print("processing file", inFile)
            if not os.path.exists(inFile):
                continue
            else:
                anyResult = True
            counter = 0
            # parse contact files. Note that C1 symmetry file is different from the rest
            for line in open(inFile):
                if counter < 8:
                    # print ("skip line", line
                    counter += 1
                else:
                    # if not self.SYMMETRY:
                    if not self.SYMMETRY or line.split()[0].startswith("/"):
                        # Second option (line.split()[0].startswith("/") stands for
                        # cases in which the result of applying symmetry is identical
                        # to the starting structure (see test testContactsSymC2_b
                        # where after deleting the #2 submodel far more than 3 A from
                        # the input model, the resulting model is the same as the initial one.
                        info = line.split()  # ['/A002', 'HEM', '1', 'ND', '/A', 'HIS', '87', 'NE2', '0.620', '2.660']
                        d1['modelId'] = "'" + "#1" + "'"
                        d1['aaName'] = "'" + info[1][0] + info[1][1:].lower() + "'"  # 'Hem'
                        d1['aaNumber'] = info[2]  # '1'
                        d1['chainId'] = "'" + info[0].split("/")[1] + "'"  # 'A002'
                        d1['atomId'] = "'" + info[3] + "'"  # 'ND'
                        d1['protId'] = "'" + labelDict[info[0].split("/")[1]] + "'"

                        d2['modelId'] = "'" + "#1" + "'"
                        d2['aaName'] = "'" + info[5][0] + info[5][1:].lower() + "'"  # 'His'
                        d2['aaNumber'] = info[6]  # '87'
                        d2['chainId'] = "'" + info[4].split("/")[1] + "'"  # 'A'
                        d2['protId'] = "'" + labelDict[info[4].split("/")[1]] + "'"
                        d2['atomId'] = "'" + info[7] + "'"  # 'NE2'
                        d['overlap'] = info[8]  # '0.620'
                        d['distance'] = info[9]  # '2.660'
                    else:
                        info = line.split()
                        # 5ni1_unit_cell_HEM.cif #1.2/A002 HEM 1 ND   5ni1_unit_cell_HEM.cif #1.2/A HIS 87 NE2    0.620    2.660
                        d1['modelId'] = "'" + info[1].split("/")[0] + "'"  # '#1.2'
                        d1['aaName'] = "'" + info[2][0] + info[2][1:].lower() + "'"  # 'Hem'
                        d1['aaNumber'] = info[3]  # '1'
                        d1['chainId'] = "'" + info[1].split("/")[1] + "'"  # 'A002'
                        d1['atomId'] = "'" + info[4] + "'"  # 'ND'
                        d1['protId'] = "'" + labelDict[info[1].split("/")[1]] + "'" # 'HEM_A'

                        d2['modelId'] = "'" + info[6].split("/")[0] + "'"  # '#1.2'
                        d2['aaName'] = "'" + info[7][0] + info[7][1:].lower() + "'"  # 'His'
                        d2['aaNumber'] = info[8]  # '87'
                        d2['chainId'] = "'" + info[6].split("/")[1] + "'"  # N
                        d2['protId'] = "'" + labelDict[info[6].split("/")[1]] + "'" # 'chainA'
                        d2['atomId'] = "'" + info[9] + "'"  # 'NE2'

                        d['overlap'] = info[10]  # '0.620'
                        d['distance'] = info[11]  # '2.660'

                    AA_1 = d1['aaName']; AA_2 = d2['aaName']
                    if ( ( (AA_1 == "'Arg'") or (AA_1 == "'Lys'")) and ((AA_2 == "'Glu'") or (AA_2=="'Asp'"))) or \
                       ( ( (AA_2 == "'Arg'") or (AA_2 == "'Lys'")) and ((AA_1 == "'Glu'") or (AA_1=="'Asp'"))) :
                        d['salineBridge'] = 1
                    else:
                        d['salineBridge'] = 0

                    if d1['modelId'] == d2['modelId']:
                        if d1['protId'] <= d2['protId']:
                            for k in d1.keys():
                                # for k in list(d1.keys()):
                                d[k + '_1'] = d1[k]
                                d[k + '_2'] = d2[k]
                        else:
                            for k in d1.keys():
                                # for k in list(d1.keys()):
                                d[k + '_1'] = d2[k]
                                d[k + '_2'] = d1[k]
                    else:
                        if d1['modelId'] <= d2['modelId']:
                            for k in d1.keys():
                                # for k in list(d1.keys()):
                                d[k + '_1'] = d1[k]
                                d[k + '_2'] = d2[k]
                        else:
                            for k in d1.keys():
                                # for k in list(d1.keys()):
                                d[k + '_1'] = d2[k]
                                d[k + '_2'] = d1[k]

                    command = "INSERT INTO contacts "
                    keys = "("
                    values = " ("
                    for key, value in d.items():
                        keys += key + ", "
                        values += str(value) + ", "
                    keys = keys[:-2] + ")"
                    values = values[:-2] + ")"

                    command += keys + " VALUES " + values
                    # print(command)
                    c.execute(command)

        return anyResult

    #    --------- util functions -----

    def getDataBaseName(self):
        return self._getExtraPath("overlaps.sqlite")

    def getSymmetrizedModelName(self):
        return os.path.abspath(self._getExtraPath("symModel.cif"))
        # return self._getExtraPath("symModel.pdb")

    def getTableName(self):
        return "contacts"

    def getView2Name(self):
        return "view_ND_2"

    def getView1Name(self):
        return "view_ND_1"

    def getChimeraScriptFileName1(self):
        return os.path.abspath(self._getTmpPath("chimera1.cxc"))

    def getChimeraScriptFileName2(self):
        return os.path.abspath(self._getTmpPath("chimera2.cxc"))

    def endChimeraScript(self, firstValue, labelDict, outFiles, f):
        protId = firstValue
        chains = ""
        comma = ''
        for k, v in labelDict.items():
            if protId == v:
                # chains += "{}/{}".format(comma, k)
                chains += "{}{}".format(comma, k)
                comma = ','
                outFileBase = v

            else:
                outFile = os.path.abspath(self._getExtraPath("{}.over".format(outFileBase)))
                # outFile = self._getExtraPath("{}.over".format(outFileBase))
                outFiles.append(outFile)
                f.write("run(session,'echo {}')\nrun(session, 'contacts  #1{} "
                         "intersubmodel true "
                         "intramol False "
                         "restrict any "
                         "saveFile {} overlapCutoff {} hbondAllowance {} namingStyle simple')\n".
                         format(chains, chains, outFile, self.cutoff, self.allowance))
                protId = v
                # chains = "/{}".format(k)
                chains = "{}".format(k)
                outFileBase = v

            chains = "/" + chains.split("/")[-1]
        outFile = os.path.abspath(self._getExtraPath("{}.over".format(outFileBase)))
        # outFile = self._getExtraPath("{}.over".format(outFileBase))
        outFiles.append(outFile)

        f.write(
            "run(session,'echo {}')\nrun(session, 'contacts  #1{} "
            "intersubmodel true "
            "intramol False "
            "restrict any "
            "savefile {} overlap {} hbond {} namingStyle simple')\n".format(
                chains, chains, outFile, self.cutoff, self.allowance))
        # f.write("run('save %s')\n" % os.path.abspath(self._getExtraPath(sessionFile)))

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
             distance,
             salineBridge
        FROM {}

        """
        commandEliminateDuplicates2 = """
        CREATE VIEW {} AS
        SELECT *
        FROM {}

        EXCEPT -- Each bound appears two times, delete one of them

        SELECT ca.*
        FROM {} ca, {} cb
        WHERE
                ca.protId_1    = cb.protId_2
            AND cb.protId_1    = ca.protId_2
            AND cb.chainId_1   = ca.chainId_2
            AND ca.aaNumber_1  = cb.aaNumber_2
            AND cb.aaNumber_1  = ca.aaNumber_2
            AND ca.atomId_1  = cb.atomId_2
            AND cb.atomId_1  = ca.atomId_2
            AND ca.modelId_2   > cb.modelId_2
        
        EXCEPT -- Interprotein bounds in the same model are not allowed

        SELECT ca.*
        FROM {} ca
        WHERE  ca.modelId_1 = ca.modelId_2 
           AND ca.protId_1 = ca.protId_2 
     
        """
        if self.SYMMETRY:
            sqlCommand = """
            SELECT count(*) FROM {} ca
            WHERE ca.modelId_1 = '#1.1'
            """.format(self.getTableName())
            c.execute(sqlCommand)
            row = c.fetchone()
            if int(row[0]) == 0:
                self.SYMMETRY = False
            else:
                commandEliminateDuplicates2 +="""
                EXCEPT -- One of the atoms must belong to the input unit cell
            
                SELECT ca.*
                FROM {} ca
                WHERE ca.modelId_1 != '#1.1'  AND 
                      ca.modelId_2 != '#1.1'
        """.format(self.getView1Name())
        # # Remove duplicate contacts
        # that is, given chains A,B
        # we have contact A.a-B.b and B.b-A.a
        c.execute(self.commandDropView.format(viewName="view_ND_1"))
        # TODO: remove second contacts
        c.execute(commandEliminateDuplicates.format("view_ND_1",
                                                    "contacts",
                                                    "contacts",
                                                    "contacts"))

        # remove duplicate contacts due to symmetry
        # h1-h1p, h1-h2p
        c.execute(self.commandDropView.format(viewName="view_ND_2"))
        c.execute(commandEliminateDuplicates2.format("view_ND_2", "view_ND_1",
                                                     "view_ND_1", "view_ND_1",
                                                     "view_ND_1"))

    def _validate(self):
        errors = []
        if self.symmetryOrder.get() <= 0:
            errors.append("Error: Symmetry Order should be a positive integer")

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
             distance float,
             salineBridge int default 0
             );"""
        c.execute(commandDropTable.format(tableName))
        c.execute(commandCreateTable.format(tableName))
    return c, conn