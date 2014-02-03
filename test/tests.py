#!/usr/bin/env python
# encoding: utf-8

from __future__ import division, print_function
import unittest
import numpy as np

import sys
sys.path.append('/home/dsuess/Documents/FHOPS/')
from hierarchy import SpectrumHierarchy
from libbath import OscillatorBath
import functions_ger as fg


class TestHierarchy(unittest.TestCase):

   def test_spectrum(self):
      psi = np.load('refspec.npy')
      dt = 628. / 9999
      corr = np.sum(psi * np.conj(psi[0][None, :]), axis=1)
      A_ref = fg.fourier(corr, dt, output_w=False, hermitian=True)

      bath = OscillatorBath(1., 1., 4.)
      h = np.array([[0., .5], [.5, 0.]], dtype=np.complex128)
      myhier = SpectrumHierarchy(628, 10000, h, bath, 4)
      A = 2*myhier.calc_trajectory()[1]
      maxdiff = np.max(np.abs(A_ref - A))
      self.assertTrue(maxdiff < 10e-2)


if __name__ == '__main__':
   unittest.main()
