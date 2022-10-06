# **************************************************************************
# *
# * Authors:    Yunior C. Fonseca Reyna (cfonseca@cnb.csic.es)
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

CHIMERA_HOME = 'CHIMERA_HOME'
ALPHAFOLD_HOME = 'ALPHAFOLD_HOME'
ALPHAFOLD_DATABASE_DIR = 'ALPHAFOLD_DATABASE_DIR'
CLUSTALO = 'clustalo'
MUSCLE = 'muscle'
CHIMERAX=True

# Supported Versions
V1_0 = '1.0'
V1_1 = '1.1'
V1_2_5 = '1.2.5'
V1_3 = '1.3'
V1_4 = '1.4'

chimeraTARs={V1_1: 'ChimeraX-1.1.tar.gz',
             V1_2_5: "chimerax-1.2.5-rc-2021.05.24",
             V1_3: "chimerax-1.3",
             V1_4: "chimerax-1.4",
}

CHIMERA_TO_SCIPION = {}
CHIMERA_CYCLIC = 0  # SYM_CYCLIC = 0
CHIMERA_DIHEDRAL_X = 1  # SYM_DIHEDRAL_X = SYM_DIHEDRAL = 1
CHIMERA_TETRAHEDRAL = 2  # SYM_TETRAHEDRAL = 3
CHIMERA_TETRAHEDRALZ3 = 3  # SYM_TETRAHEDRAL_Z3 = 4
CHIMERA_OCTAHEDRAL = 4  # SYM_OCTAHEDRAL = 5
CHIMERA_I222 = 5  # SYM_I222 = 6
CHIMERA_I222r = 6  # SYM_I222r = 7
CHIMERA_In25 = 7  # SYM_In25 = 8
CHIMERA_In25r = 8  # SYM_In25r = 9
CHIMERA_I2n3 = 9  # SYM_I2n3 = 10
CHIMERA_I2n3r = 10  # SYM_I2n3r = 11
CHIMERA_I2n5 = 11  # SYM_I2n5 = 12
CHIMERA_I2n5r = 12  # SYM_I2n5r = 13

# symmetry dictionary
# FIXME: This should not be imported here and exposed as this module constants
import pwem.constants as sciSym

# import (
#    SYM_CYCLIC, SYM_DIHEDRAL_X, SYM_TETRAHEDRAL, SYM_OCTAHEDRAL, SYM_I222,
#    SYM_I222r, SYM_In25, SYM_In25r, SYM_TETRAHEDRAL_Z3,
#    SYM_I2n3, SYM_I2n3r, SYM_I2n5, SYM_I2n5r)

CHIMERA_TO_SCIPION[CHIMERA_CYCLIC] = sciSym.SYM_CYCLIC
CHIMERA_TO_SCIPION[CHIMERA_DIHEDRAL_X] = sciSym.SYM_DIHEDRAL_X
CHIMERA_TO_SCIPION[CHIMERA_TETRAHEDRAL] = sciSym.SYM_TETRAHEDRAL
CHIMERA_TO_SCIPION[CHIMERA_TETRAHEDRALZ3] = sciSym.SYM_TETRAHEDRAL_Z3
CHIMERA_TO_SCIPION[CHIMERA_OCTAHEDRAL] = sciSym.SYM_OCTAHEDRAL
CHIMERA_TO_SCIPION[CHIMERA_I222] = sciSym.SYM_I222
CHIMERA_TO_SCIPION[CHIMERA_I222r] = sciSym.SYM_I222r
CHIMERA_TO_SCIPION[CHIMERA_In25] = sciSym.SYM_In25
CHIMERA_TO_SCIPION[CHIMERA_In25r] = sciSym.SYM_In25r
CHIMERA_TO_SCIPION[CHIMERA_I2n3] = sciSym.SYM_I2n3
CHIMERA_TO_SCIPION[CHIMERA_I2n3r] = sciSym.SYM_I2n3r
CHIMERA_TO_SCIPION[CHIMERA_I2n5] = sciSym.SYM_I2n5
CHIMERA_TO_SCIPION[CHIMERA_I2n5r] = sciSym.SYM_I2n5r

CHIMERA_SYM_NAME = dict()
CHIMERA_SYM_NAME[CHIMERA_CYCLIC] = 'Cn'
CHIMERA_SYM_NAME[CHIMERA_DIHEDRAL_X] = 'Dn'
CHIMERA_SYM_NAME[CHIMERA_TETRAHEDRAL] = 'T222'
CHIMERA_SYM_NAME[CHIMERA_TETRAHEDRALZ3] = 'TZ3'
CHIMERA_SYM_NAME[CHIMERA_OCTAHEDRAL] = 'O'
CHIMERA_SYM_NAME[CHIMERA_I222] = 'I222'
CHIMERA_SYM_NAME[CHIMERA_I222r] = 'I222r'
CHIMERA_SYM_NAME[CHIMERA_In25] = 'In25'
CHIMERA_SYM_NAME[CHIMERA_In25r] = 'In25r'
CHIMERA_SYM_NAME[CHIMERA_I2n3] = 'I2n3'
CHIMERA_SYM_NAME[CHIMERA_I2n3r] = 'I2n3r'
CHIMERA_SYM_NAME[CHIMERA_I2n5] = 'I2n5'
CHIMERA_SYM_NAME[CHIMERA_I2n5r] = 'I2n5r'

CHIMERA_CONFIG_FILE = "chimera.ini"
