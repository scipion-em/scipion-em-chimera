import numpy as np
from pyworkflow.tests import BaseTest, setupTestOutput
from chimera.objects import PAE
import json


class TestPAE(BaseTest):
    " Test PAE (Predicted Aligned Error) object"

    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)

    def test_readPAE(self):
        #create fake json file
        size = 67
        matrix = np.zeros(shape=(size, size))

        d = {}
        d['residue1'] = []
        d['residue2'] = []
        d['distance'] = []
        l = []
        for i in range(0,size):
            for j in range(0,size):
                d['residue1'].append(i + 1)
                d['residue2'].append(j + 1)
                d['distance'].append(25.*((j+1)*(j+1) + (i+1)*(i+1)*2)/(64.*64.))
                matrix[i][j] = 25.*((j+1)*(j+1) + (i+1)*(i+1)*2)/(64.*64.)
        d['max_predicted_aligned_error']=25.*(size*size + size*size)/(64.*64.)
        textFn = self.getOutputPath('data.json')
        l.append(d)
        with open(textFn, 'w') as fp:
            json.dump(l, fp)  # save dictionary as json
        
        pae = PAE(filename=textFn)
        pae.read(index=0)
        max_predicted_aligned_error = pae.getMax_predicted_aligned_error()  
        matrix = pae.getMatrix()  
        self.assertAlmostEqual(max_predicted_aligned_error, d['max_predicted_aligned_error'])
        
        # convert matrix to array
        for res1 in d['residue1']:
            for res2 in d['residue2']:
                self.assertAlmostEqual(matrix[res1-1][res2-1], d['distance'][(res1-1) * size + res2 -1])