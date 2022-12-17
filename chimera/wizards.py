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

from pwem.wizards import SelectChainWizard, SelectResidueWizard
from .protocols import ChimeraModelFromTemplate, ChimeraSubtractionMaps
from .editList import EntryGrid
from .protocols.protocol_contacts import ChimeraProtContacts
from pyworkflow.wizard import Wizard


SelectChainWizard().addTarget(protocol=ChimeraModelFromTemplate,
                              targets=['inputStructureChain'],
                              inputs=['pdbFileToBeRefined'],
                              outputs=['inputStructureChain'])

SelectChainWizard().addTarget(protocol=ChimeraSubtractionMaps,
                              targets=['inputStructureChain'],
                              inputs=['pdbFileToBeRefined'],
                              outputs=['inputStructureChain'])

SelectChainWizard().addTarget(protocol=ChimeraModelFromTemplate,
                              targets=['selectStructureChain'],
                              inputs=['pdbFileToBeRefined'],
                              outputs=['selectStructureChain'])

SelectChainWizard().addTarget(protocol=ChimeraSubtractionMaps,
                              targets=['selectStructureChain'],
                              inputs=['pdbFileToBeRefined'],
                              outputs=['selectStructureChain'])


class ProtContactsWizardChimera(Wizard):
    """ Return a table with two columns. First one is the chain id second one may
    be used by the user to group merge chains into a single object. If two
    or more chains are merged, the contacts will NOT be computed between the chains
    belonging to the same group"""
    recibingAttribute = 'chainStructure'
    _targets = [(ChimeraProtContacts, [recibingAttribute])]

    def show(self, form, *args):
        cols = ['label']
        chainWizard = SelectChainWizard()
        protocol = form.protocol
        models, modelsFirstResidue = chainWizard.getModelsChainsStep(
            protocol, protocol.pdbFileToBeRefined.get())
        rows = []
        for chainID, lenResidues in sorted(models[0].items()):
            rows.append(str(chainID))

        EntryGrid(cols, rows, form, self.recibingAttribute)

SelectResidueWizard().addTarget(protocol=ChimeraSubtractionMaps,
                                 targets=['residuesToRemove'],
                                 inputs=['pdbFileToBeRefined', ['selectStructureChain', 'inputStructureChain']],
                                 outputs=['residuesToRemove'])
