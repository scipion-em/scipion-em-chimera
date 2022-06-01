# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Roberto Marabini
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
from pwem.objects.data import EMFile
import json
import math
import numpy as np

class PAE(EMFile):
    """ In Alphafold two confidence measures has been defined.
-        The first is pLDDT (predicted lDDT-Cα), a per-residue measure of local confidence on a scale from 0 - 100. 
-        The second metric is PAE (Predicted Aligned Error), which reports AlphaFold’s expected position error 
-        at residue x, when the predicted and true structures are aligned on residue y. 
-    
-        PAE is stored in a json file with the follwoing format:
-        
-        For a protein of length num_res, the JSON has the following structure of arrays format:
-        [
-            {
-                "residue1": [ 1, 1, 1, 1, 1, ...], # Length: num_res^2.
-                "residue2": [ 1, 2, 3, 4, 5, ...], # Length: num_res^2.
-                "distance": [0.2, 1.2, 3.7, 6.6, 8.7, ...], # Length: num_res^2.
-                "max_predicted_aligned_error": 31.75
-            }
-        ]
-
-        The fields in the JSON are:
-
-            residue1: The residue on which the structure is aligned for the predicted error.
-            residue2: The residue on which the error is predicted.
-            distance: The PAE value of the residue pair.
-            max_predicted_aligned_error: A scalar that denotes the maximum possible value of PAE. The minimum PAE is 0.
-
-        """
    matrix = None  # PAE data as a 2D array of floats
    max_predicted_aligned_error = None

    def __init__(self, filename=None, pseudoatoms=False, **kwargs):
        EMFile.__init__(self, filename, **kwargs)

    def read(self, index=0):
        """read json file and assign content to matrix and max_predicted_aligned_error"""
        with open(self.getFileName(), 'r') as f:
            data = json.load(f)
            data = data[index]
            if 'max_predicted_aligned_error' in data.keys():
                self.max_predicted_aligned_error = data['max_predicted_aligned_error']
            else:
                self.max_predicted_aligned_error = None
            residue1 = data['residue1']
            residue2 = data['residue2']
            distance = data['distance']
            size = int(math.sqrt(len(residue1)))
            self.matrix = np.zeros(shape=(size, size))
            for res1, res2, dis  in zip(residue1, residue2, distance):
                self.matrix[res1-1][res2-1] = dis

    def getMatrix(self):
        if self.matrix is None:
            raise("Error, matrix is None")
        return self.matrix

    def getMax_predicted_aligned_error(self):
        if self.max_predicted_aligned_error is None:
            raise("Error, max_predicted_aligned_error is None")
        return self.max_predicted_aligned_error

