# **************************************************************************
# *
# * Authors:     Pablo Conesa (pconesa@cnb.csic.es)
#                Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

from pwem.wizards import GetStructureChainsWizard, pwobj, emconv
from .protocols import ChimeraModelFromTemplate, ChimeraSubtractionMaps
from .editList import EntryGrid
from .protocols.protocol_contacts import ChimeraProtContacts
from pyworkflow.wizard import Wizard
from pyworkflow.gui.tree import ListTreeProviderString
from pyworkflow.gui import dialog


class GetStructureChainsWizardChimera(GetStructureChainsWizard):
    _targets = [(ChimeraModelFromTemplate, ['inputStructureChain']),
                (ChimeraSubtractionMaps, ['inputStructureChain'])]

class GetStructureChains2WizardChimera(GetStructureChainsWizard):
    _targets = [(ChimeraSubtractionMaps, ['selectStructureChain']),
                (ChimeraModelFromTemplate, ['selectStructureChain'])]

    def show(self, form, *params):
        protocol = form.protocol
        try:
            modelsLength, modelsSeq = self.getModelsChainsStep(protocol)
        except Exception as e:
            print("ERROR: ", e)
            return

        self.editionListOfChains(modelsLength)
        finalChainList = []
        for i in self.chainList:
            finalChainList.append(pwobj.String(i))
        provider = ListTreeProviderString(finalChainList)
        dlg = dialog.ListDialog(form.root, "Model chains", provider,
                                "Select one of the chains (model, chain, "
                                "number of chain residues)")

        form.setVar('selectStructureChain', dlg.values[0].get())


class ProtContactsWizardChimera(Wizard):
    """ Return a table with two columns. First one is the chain id second one may
    be used by the user to group merge chains into a single object. If two
    or more chains are merged, the contacts will NOT be computed between the chains
    belonging to the same group"""
    recibingAttribute = 'chainStructure'
    _targets = [(ChimeraProtContacts, [recibingAttribute])]

    def show(self, form, *args):
        cols = ['label']
        chainWizard = GetStructureChainsWizard()
        protocol = form.protocol
        models, modelsFirstResidue = chainWizard.getModelsChainsStep(protocol)
        rows = []
        for chainID, lenResidues in sorted(models[0].items()):
            rows.append(str(chainID))

        EntryGrid(cols, rows, form, self.recibingAttribute)

class GetChainResiduesWizardChimera(GetStructureChainsWizard):
    _targets = [(ChimeraSubtractionMaps, ['firstResidueToRemove'])]

    def editionListOfResidues(self, modelsFirstResidue, model, chain):
        self.residueList = []
        for modelID, chainDic in modelsFirstResidue.items():
            if int(model) == modelID:
                for chainID, seq_number in chainDic.items():
                    if chain == chainID:
                        for i in seq_number:
                            self.residueList.append(
                                '{"residue": %d, "%s"}' % (i[0], str(i[1])))

    def getResidues(self, form):
        protocol = form.protocol
        try:
            modelsLength, modelsFirstResidue = self.getModelsChainsStep(protocol)
        except Exception as e:
            print("ERROR: ", e)
            return
        if protocol.selectStructureChain.get() is not None:
            selection = protocol.selectStructureChain.get()
        else:
            selection = protocol.inputStructureChain.get()
        model = selection.split(',')[0].split(':')[1].strip()
        chain = selection.split(',')[1].split(':')[1].split('"')[1]
        self.editionListOfResidues(modelsFirstResidue, model, chain)
        finalResiduesList = []
        for i in self.residueList:
            finalResiduesList.append(pwobj.String(i))
        return finalResiduesList

    def show(self, form, *params):
        finalResiduesList = self.getResidues(form)
        provider = ListTreeProviderString(finalResiduesList)
        dlg = dialog.ListDialog(form.root, "Chain residues", provider,
                                "Select one residue (residue number, "
                                "residue name)")
        form.setVar('firstResidueToRemove', dlg.values[0].get())

class GetChainResidues2WizardChimera(GetChainResiduesWizardChimera):
    _targets = [(ChimeraSubtractionMaps, ['lastResidueToRemove'])]

    def show(self, form, *params):
        finalResiduesList = self.getResidues(form)
        provider = ListTreeProviderString(finalResiduesList)
        dlg = dialog.ListDialog(form.root, "Chain residues", provider,
                                "Select one residue (residue number, "
                                "residue name)")
        form.setVar('lastResidueToRemove', dlg.values[0].get())

