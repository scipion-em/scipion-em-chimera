# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
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

from .protocol_base import ChimeraProtBase


class ChimeraProtRigidFit(ChimeraProtBase):
    # """Protocol to perform rigid fit using Chimera.
    #     Execute command *scipionwrite [model #n] [refmodel #p]
    #     [saverefmodel 0|1]* from command line in order to transferm fitted
    #     pdb to scipion. Default values are model=#0,
    #     refmodel =#1 and saverefmodel 0 (false).
    #     model refers to the pdb file. refmodel to a 3Dmap"""
    """Protocol to perform rigid fit using Chimera.
            Execute command *scipionwrite #n [prefix stringAddedToFilename]
            model refers to the pdb file"""
    _label = 'rigid fit'

    def _defineParams(self, form):
        super(ChimeraProtRigidFit, self)._defineParams(form)
        param = form.getParam('pdbFileToBeRefined')
        param.label.set('Atomic structure to be fitted')
        param.help.set('PDBx/mmCIF file to be fitted. ')
        param = form.getParam('inputPdbFiles')
        param.label.set('Other reference atomic structures')
        param.help.set('Other PDBx/mmCIF files used as reference.')
        param.allowsNull.set('True')

        # --------------------------- INSERT steps functions --------------------

    def prerequisitesStep(self):
        """
        """
        if self.inputVolume.get() is None:
            fnVol = self.pdbFileToBeRefined.get().getVolume()
            index, fn = fnVol.getLocation()
            print("Volume: Volume associated to atomic structure %s(%d)\n"
                  % (fn, index))
        else:
            fnVol = self.inputVolume.get()
            # print("Volume: Input volume %s\n" % fnVol)

    def _validate(self):
        errors = super(ChimeraProtRigidFit, self)._validate()

        # Check that the input volume exist
        if (not self.pdbFileToBeRefined.get().hasVolume()) and (
                    self.inputVolume.get() is None):
            errors.append("Error: You should provide a volume.\n")

        return errors
