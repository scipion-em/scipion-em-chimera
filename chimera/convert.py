# **************************************************************************
# *
# * Authors:    Marta Martinez (mmmtnez@cnb.csic.es)
# *             Roberto Marabini (roberto@cnb.csic.es)
# *
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
from pwem.constants import (SYM_DIHEDRAL_X,SCIPION_SYM_NAME)
from .constants import (CHIMERA_SYM_NAME, CHIMERA_TO_SCIPION, CHIMERA_CYCLIC,
                         CHIMERA_DIHEDRAL_X, CHIMERA_TETRAHEDRAL,
                         CHIMERA_TETRAHEDRALZ3, CHIMERA_OCTAHEDRAL,
                         CHIMERA_I222, CHIMERA_I222r, CHIMERA_In25, CHIMERA_In25r,
                         CHIMERA_I2n5, CHIMERA_I2n5r, CHIMERA_I2n3, CHIMERA_I2n3r,
                         )
CHIMERA_LIST = [CHIMERA_SYM_NAME[CHIMERA_CYCLIC] +
                               " (" + SCIPION_SYM_NAME[CHIMERA_TO_SCIPION[CHIMERA_CYCLIC]] + ")",
                               CHIMERA_SYM_NAME[CHIMERA_DIHEDRAL_X] +
                               " (" + SCIPION_SYM_NAME[CHIMERA_TO_SCIPION[CHIMERA_DIHEDRAL_X]] + ")",
                               CHIMERA_SYM_NAME[CHIMERA_TETRAHEDRAL] +
                               " (" + SCIPION_SYM_NAME[CHIMERA_TO_SCIPION[CHIMERA_TETRAHEDRAL]] + ")",
                               CHIMERA_SYM_NAME[CHIMERA_TETRAHEDRALZ3] +
                               " (" + SCIPION_SYM_NAME[CHIMERA_TO_SCIPION[CHIMERA_TETRAHEDRALZ3]] + ")",
                               CHIMERA_SYM_NAME[CHIMERA_OCTAHEDRAL] +
                               " (" + SCIPION_SYM_NAME[CHIMERA_TO_SCIPION[CHIMERA_OCTAHEDRAL]] + ")",
                               CHIMERA_SYM_NAME[CHIMERA_I222] +
                               " (" + SCIPION_SYM_NAME[CHIMERA_TO_SCIPION[CHIMERA_I222]] + ")",
                               CHIMERA_SYM_NAME[CHIMERA_I222r] +
                               " (" + SCIPION_SYM_NAME[CHIMERA_TO_SCIPION[CHIMERA_I222r]] + ")",
                               CHIMERA_SYM_NAME[CHIMERA_In25] +
                               " (" + SCIPION_SYM_NAME[CHIMERA_TO_SCIPION[CHIMERA_In25]] + ")",
                               CHIMERA_SYM_NAME[CHIMERA_In25r] +
                               " (" + SCIPION_SYM_NAME[CHIMERA_TO_SCIPION[CHIMERA_In25r]] + ")",
                               CHIMERA_SYM_NAME[CHIMERA_I2n3] +
                               " (" + SCIPION_SYM_NAME[CHIMERA_TO_SCIPION[CHIMERA_I2n3]] + ")",
                               CHIMERA_SYM_NAME[CHIMERA_I2n3r] +
                               " (" + SCIPION_SYM_NAME[CHIMERA_TO_SCIPION[CHIMERA_I2n3r]] + ")",
                               CHIMERA_SYM_NAME[CHIMERA_I2n5] +
                               " (" + SCIPION_SYM_NAME[CHIMERA_TO_SCIPION[CHIMERA_I2n5]] + ")",
                               CHIMERA_SYM_NAME[CHIMERA_I2n5r] +
                               " (" + SCIPION_SYM_NAME[CHIMERA_TO_SCIPION[CHIMERA_I2n5r]] + ")",
                               ]