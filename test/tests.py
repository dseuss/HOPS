#!/usr/bin/env python
# encoding: utf-8

from __future__ import division, print_function
import unittest
import numpy as np

import sys
sys.path.append('..')
from hierarchy import SpectrumHierarchy
from libbath import OscillatorBath
import functions_ger as fg

import matplotlib.pyplot as plt


class TestHierarchy(unittest.TestCase):

   def test_spectrum(self):
      psi = np.load('refspec.npy')
      dt = 628. / 9999
      corr = np.sum(psi * np.conj(psi[0][None, :]), axis=1)
      w, A_ref = fg.fourier(corr, dt, output_w=True, hermitian=True)

      bath = OscillatorBath(1., 1., 4.)
      h = np.array([[0., .5], [.5, 0.]], dtype=np.complex128)
      myhier = SpectrumHierarchy(628, 10000, h, bath, 4)
      A = myhier.calc_trajectory()[1]
      maxdiff = np.max(np.abs(A_ref - A))
      self.assertTrue(maxdiff < 5*10e-2)

if __name__ == '__main__':
   unittest.main()
