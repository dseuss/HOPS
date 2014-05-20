#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
from __future__ import print_function

import numpy as np
import functions_ger as fg


class OscillatorBath(object):
    class InvalidBCFError(Exception):
        pass

    def __init__(self, g, gamma, Omega, antisym=False, T=0., ncoth=1,
          alpha=None):
        """
             Sets up the parameters for a harmonic bath.

             If antisym is False, we expect the parameters to be given directly
             for a bath correlation function of the form

                α(t) = Σ_j  g_j * exp(-γ_j * |t|  -  i * Ω * t)

             In that case no expansion for non-zero temperature is possible.

             For antisym=True we expect the parameters to belong to a sum of
             anti-symmetrized Lorentzians of the form

                          gγ             1                  1
                  J(ω) = ---- * ( ---------------- - ---------------- )
                          π        (ω - Ω )^2 + γ^2   (ω + Ω )^2 + γ^2

            These can be expanded for non-zero temperature in a Padé
            decomposition. ncoth is the number of extra terms for the Temperature
            expansion. In total we get (2 + ncoth) * len(g) expansion terms.
        """
        assert np.shape(g) == np.shape(gamma)
        assert np.shape(g) == np.shape(Omega)
        assert T >= 0.

        if not antisym:
            self._g = np.ravel(g)
            self._gamma = np.ravel(gamma)
            self._Omega = np.ravel(Omega)
        elif (T == 0.):
            self._gamma = np.ravel(gamma)
            self._Omega = np.ravel(Omega)
            self._g = np.pi / (4 * self.gamma * self.Omega) * np.ravel(g)
        elif (T > 0.):
            poles, residues = fg.coth_poles_pade(ncoth)
            self._g, self._gamma, self._Omega = fg.spd2bcf(np.ravel(g),
                  np.ravel(gamma), np.ravel(Omega), T, poles, residues)
        else:
            raise self.InvalidBCFError('INVALID COMBINATION FOR BCF!')

        if alpha is None:
            self._alpha = self.makeCorrelationFunction
        else:
            self._alpha = alpha

    @property
        def g(self):
    return self._g

    @property
    def gamma(self):
        return self._gamma

    @property
    def Omega(self):
        return self._Omega

    @property
    def nrmodes(self):
        return self._g.size

    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, value):
        self._alpha = value

    def modes(self):
        for i in range(self._g.size):
            yield {'g': self._g[i],
                   'gamma': self._gamma[i],
                   'Omega': self._Omega[i]}

    ## For spectral density function
    def __str__(self):
        mode_str = ['( {:2.2e} * exp(-{:2.2e} * |t| - i * {:2.2e} * t) )'
              .format(m['g'], m['gamma'], m['Omega']) for m in self.modes()]
        return '\n + '.join(mode_str)

    def SpectralDensity(self, w):
        """
           Returns spectral density as function of w
           w -- np-array
        """
        sd = [m['g'] * m['gamma'] / np.pi / ((w - m['Omega'])**2 + m['gamma']**2)
              for m in self.modes()]
        return np.sum(np.array(sd), axis=0)

    ## For bath correlation function
    def makeCorrelationFunction(self, t):
        """
           Generates the bath correlation function (as sum of exps) at t.
        """
        bcf = [m['g'] * np.exp(-m['gamma'] * np.abs(t) - 1.j * m['Omega'] * t)
              for m in self.modes()]
        return np.sum(np.array(bcf), axis=0)

    def create_copies(self, copies):
        """@todo: Docstring for create_copies.

        :copies: @todo
        :returns: @todo

        """
        g = np.array(self._g.tolist() * copies)
        gamma = np.array(self._gamma.tolist() * copies)
        Omega = np.array(self._Omega.tolist() * copies)

        Map = []
        for i in range(copies):
            Map.extend(i * np.ones(self.nrmodes, dtype=int))

        return g, gamma, Omega, Map
