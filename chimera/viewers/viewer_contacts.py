from pyworkflow.protocol.params import EnumParam, BooleanParam, \
    LabelParam
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from chimera.protocols.protocol_contacts import ProtContacts
from pyworkflow.gui.text import _open_cmd

class ProtContactsViewer(ProtocolViewer):
    _label = 'Contacts Viewer'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ProtContacts]

    def __init__(self,  **kwargs):
        ProtocolViewer.__init__(self,  **kwargs)
        self.c, self.conn = self.protocol.prepareDataBase(drop=False)
        # compute all pairs of chains that interact
        # this information is needed for the menu
        self.pairChains = self._displayPairChains()

    def _defineParams(self, form):
        form.addSection(label="Chains")
        form.addParam("doInvert", BooleanParam, label="swap chain columns",
                      default=False,
                      help="Set to TRUE to swap the first by the second chain.\n")
        form.addParam('chainPair', EnumParam,
                      choices=self.pairChains,
                      default=1,
                      label="Chain Contacts",
                      help="syntax: modelName1, complexName1, chainName1"
                           "        modelName2, complexName2, chainName2"
                      )

        form.addParam('displayPairChains', LabelParam,
                      label="show file with Interacting Chains",
                      help="Display the chains that interact, and "
                           "the # of atoms that interact")

    def _getVisualizeDict(self):
        return{
            'chainPair': self._chainPair,
            'displayPairChains': self._visualizeChainPairFile
        }

    def _visualizeChainPairFile(self, e=None):
        """Show file with the chains that interact."""
        _open_cmd(self.getPairChainsFileNme(), self.getTkRoot())

    def _chainPair(self, e=None):
        if self.doInvert.get():
            commandDisplayInteractions = """
SELECT count(*), protId_2, modelId_2, chainId_2,
       aaName_2 || aaNumber_2, protId_1, modelId_1,
       chainId_1,  aaName_1 || aaNumber_1
FROM view_ND_2
WHERE modelId_1='{}' AND protId_1='{}' AND chainId_1='{}'
  AND modelId_2='{}' AND protId_2='{}' AND chainId_2='{}'
GROUP BY protId_1, modelId_1, chainId_1,
         aaNumber_1, aaName_1, protId_2,
         modelId_2, chainId_2, aaNumber_2,
         aaName_2
ORDER BY protId_2, modelId_2, chainId_2,
        protId_1, modelId_1, chainId_1,
        aaNumber_2, aaName_2,  aaNumber_1, aaName_1;
"""
        else:
            commandDisplayInteractions = """
SELECT count(*), protId_1, modelId_1, chainId_1,
       aaName_1 || aaNumber_1, protId_2, modelId_2,
       chainId_2,  aaName_2 || aaNumber_2
FROM   view_ND_2
WHERE modelId_1='{}' AND protId_1='{}' AND chainId_1='{}'
  AND modelId_2='{}' AND protId_2='{}' AND chainId_2='{}'
GROUP BY protId_1, modelId_1, chainId_1,
         aaNumber_1, aaName_1, protId_2,
         modelId_2, chainId_2, aaNumber_2,
         aaName_2
ORDER BY protId_1, modelId_1, chainId_1,
         protId_2, modelId_2, chainId_2,
         aaNumber_1, aaName_1,  aaNumber_2, aaName_2;
"""
        row = self.all_pair_chains[self.chainPair.get()]
        comamnd = commandDisplayInteractions.format(row[1],
                                                    row[2],
                                                    row[3],
                                                    row[4],
                                                    row[5],
                                                    row[6]
                                                    )
        #print comamnd

        self.c.execute(comamnd)
        all_rows = self.c.fetchall()
        f = open(self.getInteractionFileName(), 'w')
        f.write("RESULTS for: {}\n".format(', '.join(str(s) for s in row)))
        f.write("# atoms, prot_1, model_1, chain_1, AA_1, prot_2, model_2, chain2, AA_2\n")
        for row in all_rows:
            f.write(', '.join(str(s) for s in row) + "\n")
        f.close()
        _open_cmd(self.getInteractionFileName(), self.getTkRoot())


    def _displayPairChains(self,):
        # auxiliary view name
        viewPairChain = 'atoms_interacting_per_pair_of_chains'

        # drop auxiliary view if exists
        sqlCommand = self.protocol.commandDropView.\
            format(viewName=viewPairChain)
        self.c.execute(sqlCommand)

        # create view with pair of chains
        commandDisplayPairChains="""
CREATE VIEW {viewName} AS
SELECT count(*) as AAs,  modelId_1, protId_1, chainId_1, modelId_2,  protId_2,  chainId_2
FROM view_ND_2
GROUP BY modelId_1, protId_1, chainId_1, modelId_2, protId_2,  chainId_2
-- ORDER BY modelId_1, protId_1, chainId_1, modelId_2, protId_2,  chainId_2;
""".format(viewName=viewPairChain)
        self.c.execute(commandDisplayPairChains)

        # remove duplicates and execute command
        commandDisplayPairChainsNR="""
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
        self.c.execute(commandDisplayPairChainsNR)
        self.all_pair_chains = self.c.fetchall()

        # create text file and list with pairs of chains
        f = open(self.getPairChainsFileNme(), 'w')
        formatted_row = '{:<4} {:>3} {:<11} {:<3} {:>4} {:<11} {:<3}\n'
        f.write(        "# atoms, model_1, prot_1, chain_1,  model_2, prot_2, chain_2\n")

        choices = []
        for row in self.all_pair_chains:
            f.write(formatted_row.format(*row))
            choices.append("{model_1},{prot_1},{chain_1}:"
                           "{model_2},{prot_2},{chain_2}".format(model_1=row[1],
                                                       prot_1=row[2],
                                                       chain_1=row[3],
                                                       model_2=row[4],
                                                       prot_2=row[5],
                                                       chain_2=row[6])
                           )
        f.close()
        return choices # list with pairs of chains

    def getPairChainsFileNme(self):
        return self.protocol._getTmpPath('pairChainsFile.txt')

    def getInteractionFileName(self):
        return self.protocol._getTmpPath('interactions.txt')