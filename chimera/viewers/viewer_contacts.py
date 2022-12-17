from pwem.viewers import Chimera
from pyworkflow.protocol.params import EnumParam, BooleanParam, \
    LabelParam, IntParam
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from pwem import Domain

from ..protocols.protocol_contacts import ChimeraProtContacts
from pyworkflow.gui.text import _open_cmd
import os


class ChimeraProtContactsViewer(ProtocolViewer):
    _label = 'Contacts Viewer'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ChimeraProtContacts]

    def __init__(self, **kwargs):

        ProtocolViewer.__init__(self, **kwargs)
        self.c, self.conn = self.protocol.prepareDataBase(drop=False)
        # compute all pairs of chains that interact
        # this information is needed for the menu
        self.pairChains = self._displayPairChains()

    def _defineParams(self, form):
        form.addSection(label="Display Results")
        group = form.addGroup('3D Visualization')
        group.addParam('displayModel', LabelParam,
                       label="View models in ChimeraX",
                       help="Display of input atomic structure and its respective"
                            " symmetrized models.")
        group = form.addGroup('Interacting chains')
        group.addParam('displayPairChains', LabelParam,
                       label="Summary list of all Interacting Chains",
                       help="Display the interacting chains, and "
                            "the number of atoms involved in these interactions.")
        group = form.addGroup('Contacts between interacting chains')
        group.addParam("doInvert", BooleanParam, label="Swap chain columns in the summary of contacts",
                       default=False,
                       help="Set to YES to swap the first chain by the second one.\n")
        group.addParam('aaDistance', IntParam,
                       label='Distance to group residues (Number of residues)',
                       default=4,
                       help='If two residues are closer than this distance (number of residues),'
                            ' then those two residues will be grouped.')
        group.addParam('chainPair', EnumParam,
                       choices=self.pairChains,
                       default=0,
                       label="Select two interacting chains and get the summary of contacts",
                       help="Format of the interacting chains in the display window:\n"
                            "# modelName1, chainLabelName1, chainName1 # modelName2, "
                            "chainLabelName2, chainName2\n\n"
                            "Output list format:\n"
                            "Title: 'RESULTS for: # modelName1, chainLabelName1, chainName1 "
                            "# modelName2, chainLabelName2, chainName2'\n"
                            "Below the title, several paragraphs grouping residues are"
                            " shown according to the distance to group residues selected by "
                            "the user. Columns in each paragraph:\nnumberOfatoms, "
                            "chainLabelName1, modelName1, chainName1, residueName1, "
                            "chainLabelName2, modelName2, chainName2, residueName2.\n"
                            "Meaning of question marks below each group:\n'?': First and last "
                            "residues are separated by more than 20 residues.\n'????':  First "
                            "and last residues are separated by less than 20 residues.")

    def _getVisualizeDict(self):
        return {
            'displayModel': self._displayModel,
            'chainPair': self._chainPair,
            'displayPairChains': self._visualizeChainPairFile
        }

    def _displayModel(self, e=None):
        # bildFileName = os.path.abspath(self.protocol._getTmpPath(
        #    "axis_output.bild"))
        bildFileName = self.protocol._getExtraPath("axis_output.bild")

        # Axis Dim
        dim = 150.
        sampling = 1.
        Chimera.createCoordinateAxisFile(dim,
                                         bildFileName=bildFileName,
                                         sampling=sampling)

        fnCmd = self.protocol._getExtraPath("chimera_output.cxc")
        f = open(fnCmd, 'w')
        # change to workingDir
        # If we do not use cd and the project name has an space
        # the protocol fails even if we pass absolute paths
        f.write('cd %s\n' % os.getcwd())
        # reference axis model = 0
        f.write("open %s\n" % bildFileName)
        f.write("cofr 0,0,0\n")  # set center of coordinates
        # f.write("open %s\n" % os.path.abspath(
        #     self.protocol.pdbFileToBeRefined.get().getFileName()))
        f.write("open %s\n" % self.protocol.pdbFileToBeRefined.get().getFileName())

        if self.protocol.SYMMETRY.get() and \
                os.path.exists(self.protocol.getSymmetrizedModelName()):
            f.write("open %s\n" % self.protocol.getSymmetrizedModelName())
        f.close()
        # run in the background
        chimeraPlugin = Domain.importFromPlugin('chimera', 'Plugin', doRaise=True)
        chimeraPlugin.runChimeraProgram(chimeraPlugin.getProgram(), fnCmd + "&",
                                        cwd=os.getcwd())
        return []

    def _visualizeChainPairFile(self, e=None):
        """Show file with the chains that interact."""
        _open_cmd(self.getPairChainsFileName(), self.getTkRoot())

    def _chainPair(self, e=None):
        if self.doInvert.get():
            commandDisplayInteractions = """
SELECT count(*), protId_2, modelId_2, chainId_2,
       aaName_2 || aaNumber_2, protId_1, modelId_1,
       chainId_1,  aaName_1 || aaNumber_1, salineBridge
FROM view_ND_2
WHERE modelId_1='{}' AND protId_1='{}' AND chainId_1='{}'
  AND modelId_2='{}' AND protId_2='{}' AND chainId_2='{}'
GROUP BY protId_1, modelId_1, chainId_1,
         aaNumber_1, aaName_1, protId_2,
         modelId_2, chainId_2, aaNumber_2,
         aaName_2, salineBridge
ORDER BY protId_2, modelId_2, chainId_2,
        protId_1, modelId_1, chainId_1,
        aaNumber_2, aaName_2,  aaNumber_1, aaName_1;
"""
        else:
            commandDisplayInteractions = """
SELECT count(*), protId_1, modelId_1, chainId_1,
       aaName_1 || aaNumber_1, protId_2, modelId_2,
       chainId_2,  aaName_2 || aaNumber_2, salineBridge
FROM   view_ND_2
WHERE modelId_1='{}' AND protId_1='{}' AND chainId_1='{}'
  AND modelId_2='{}' AND protId_2='{}' AND chainId_2='{}'
GROUP BY protId_1, modelId_1, chainId_1,
         aaNumber_1, aaName_1, protId_2,
         modelId_2, chainId_2, aaNumber_2,
         aaName_2, salineBridge
ORDER BY protId_1, modelId_1, chainId_1,
         protId_2, modelId_2, chainId_2,
         aaNumber_1, aaName_1,  aaNumber_2, aaName_2;
"""
        f = open(self.getInteractionFileName(), 'w')
        if len(self.all_pair_chains) == self.chainPair.get():
            f.write("No contacts found by applying symmetry: Is the symmetry "
                    "center equal to the origin of coordinates?")
        else:
            row = self.all_pair_chains[self.chainPair.get()]
            command = commandDisplayInteractions.format(row[1],
                                                        row[2],
                                                        row[3],
                                                        row[4],
                                                        row[5],
                                                        row[6]
                                                        )
            # print command
            rows_count = self.c.execute(command)
            all_rows = self.c.fetchall()

            f.write("RESULTS for: {}\n".format(', '.join(str(s) for s in row)))
            f.write("# atoms, prot_1, model_1, chain_1, AA_1, prot_2, model_2, chain2, AA_2 salineBridge\n")
            first = None
            last = None
            first2 = None
            last2 = None
            aaDistance = self.aaDistance.get()
            for row in all_rows:
                AA_1 = row[4]
                AA_1Int = int(AA_1[3:])
                AA_2 = row[8]
                AA_2Int = int(AA_2[3:])
                if first is None:
                    first = AA_1
                    last = first
                    lastInt = AA_1Int

                    first2 = AA_2
                    firstInt2 = AA_2Int
                    last2 = first2
                    lastInt2 = AA_2Int
                else:
                    if (AA_1Int - lastInt) > aaDistance:

                        f.write(">>>> {first}".format(first=first))
                        if last != first:
                            f.write("_{last}".format(last=last))
                        f.write(" ---- {first}".format(first=first2))
                        if last2 != first2:
                            f.write("_{last}".format(last=last2))
                        if (lastInt2 - firstInt2) > 20:
                            f.write("????\n\n")
                        else:
                            f.write("?\n\n")

                        first = AA_1
                        last = first
                        lastInt = AA_1Int
                        first2 = AA_2
                        firstInt2 = AA_2Int
                        last2 = AA_2
                        lastInt2 = AA_2Int
                    else:
                        last = AA_1
                        lastInt = int(AA_1[3:])
                        tmpLastInt2 = int(AA_2[3:])
                        if lastInt2 < tmpLastInt2:
                            last2 = AA_2
                            lastInt2 = tmpLastInt2
                        if firstInt2 > tmpLastInt2:
                            first2 = AA_2
                            firstInt2 = tmpLastInt2

                f.write(', '.join(str(s) for s in row) + "\n")
                # print first, last, first2, last2, firstInt2, lastInt2
            f.write(">>>> {first}".format(first=first))
            if last != first:
                f.write("_{last}".format(last=last))
            f.write(" ---- {first}".format(first=first2))
            if last2 != first2:
                f.write("_{last}".format(last=last2))
            if (lastInt2 - firstInt2) > 20:
                f.write("????\n\n")
            else:
                f.write("?\n\n")

        f.close()
        _open_cmd(self.getInteractionFileName(), self.getTkRoot())

    def _displayPairChains(self, ):
        # auxiliary view name
        viewPairChain = 'atoms_interacting_per_pair_of_chains'

        # drop auxiliary view if exists
        sqlCommand = self.protocol.commandDropView. \
            format(viewName=viewPairChain)
        self.c.execute(sqlCommand)

        # create view with pair of chains
        commandDisplayPairChains = """
CREATE VIEW {viewName} AS
SELECT count(*) as AAs,  modelId_1, protId_1, chainId_1, modelId_2,  protId_2,  chainId_2
FROM view_ND_2
GROUP BY modelId_1, protId_1, chainId_1, modelId_2, protId_2,  chainId_2
-- ORDER BY modelId_1, protId_1, chainId_1, modelId_2, protId_2,  chainId_2;
""".format(viewName=viewPairChain)
        self.c.execute(commandDisplayPairChains)

        # remove duplicates and execute command
        commandDisplayPairChainsNR = """
SELECT *
FROM {viewName}

EXCEPT

SELECT ca.*
FROM {viewName} ca, {viewName} cb
WHERE
      ca.protId_1    = cb.protId_2
  AND cb.protId_1    = ca.protId_2
  AND ca.chainId_1   = cb.chainId_2
  AND cb.chainId_1   = ca.chainId_2
  AND ca.AAs  = cb.AAs
  AND ca.protId_1 > cb.protId_1

ORDER BY modelId_1, protId_1, chainId_1, modelId_2, protId_2,  chainId_2;
""".format(viewName=viewPairChain)
        rows_count = self.c.execute(commandDisplayPairChainsNR)

        # create text file and list with pairs of chains
        f = open(self.getPairChainsFileName(), 'w')
        choices = []

        self.all_pair_chains = self.c.fetchall()

        formatted_row = '{:<4} {:>3} {:<11} {:<3} {:>4} {:<11} {:<3}\n'
        f.write("# atoms, model_1, prot_1, chain_1,  model_2, prot_2, chain_2\n")

        for row in self.all_pair_chains:
            f.write(formatted_row.format(*row))
            choices.append("{model_1},{prot_1},{chain_1}:"
                           "{model_2},{prot_2},{chain_2}".format(model_1=row[1],
                                                                 prot_1=row[2],
                                                                 chain_1=row[3],
                                                                 model_2=row[4],
                                                                 prot_2=row[5],
                                                                 chain_2=row[6]
                                                                 )
                           )
        if self.protocol.SYMMETRY.get() and not \
                os.path.exists(self.protocol.getSymmetrizedModelName()):
            f.write("No contacts found by applying symmetry: Is the symmetry "
                    "center equal to the origin of coordinates?\n")
            choices.append("No contacts found by applying symmetry: Is the symmetry "
                           "center equal to the origin of coordinates?")
        f.close()
        return choices  # list with pairs of chains

    def getPairChainsFileName(self):
        return self.protocol._getExtraPath('pairChainsFile.txt')

    def getInteractionFileName(self):
        return self.protocol._getExtraPath('interactions.txt')
