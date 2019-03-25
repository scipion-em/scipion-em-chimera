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

from pyworkflow.em.wizards import GetStructureChainsWizard
from chimera.protocols import ChimeraModelFromTemplate
from editList import EntryGrid
from chimera.protocols.protocol_contacts import \
    ChimeraProtContacts
from pyworkflow.wizard import Wizard

class GetStructureChainsWizardChimera(GetStructureChainsWizard):
    _targets = [(ChimeraModelFromTemplate, ['inputStructureChain'])]

class ProtContactsWizardChimera(Wizard):
    recibingAttribute = 'chainStructure'
    _targets = [(ChimeraProtContacts, [recibingAttribute])
                ]

    def show(self, form):
        cols = ['label']
        chainWizard = GetStructureChainsWizard()
        protocol = form.protocol
        models = chainWizard.getModelsChainsStep(protocol)
        rows = []
        for chainID, lenResidues in sorted(models[0].iteritems()):
            rows.append(str(chainID))

        EntryGrid(cols, rows, form , self.recibingAttribute)
