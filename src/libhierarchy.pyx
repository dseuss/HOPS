import numpy as np
cimport numpy as np
import cython

cdef extern:
   void c_init(double *tLength, int *tSteps, int *depth, int *modes,
               int *hs_dim, double complex *g, double *gamma, double *Omega,
               double complex *h, int *Lmap, bint *with_terminator,
               int *populated_modes)
   void c_run_trajectory_z0_rk4(int *hs_dim, int *tSteps, complex *psi0,
                                double complex *psi)


cdef class FHierarchy(object):

   """Docstring for FHierarchy. """

   cdef int _tSteps
   cdef int _dim

   def __init__(self, double tLength, int tSteps, int depth, g, gamma, Omega,
                h, Lmap, bint with_terminator=True, populated_modes=None):
      """@todo: to be defined1. """
      cdef double complex[:] gc = np.array(g, dtype=np.complex128, order='F')
      cdef double[:] gammac = np.array(gamma, dtype=np.float64, order='F')
      cdef double[:] Omegac = np.array(Omega, dtype=np.float64, order='F')
      cdef double complex[:, :] hc = np.array(h, dtype=np.complex128, order='F')
      cdef int[:] Lmapc = np.array(Lmap, dtype=np.int32, order='F')

      self._dim = h.shape[0]
      self._tSteps = tSteps
      cdef int populated_modesc
      cdef int modes = g.size
      assert(modes == gamma.size)
      assert(modes == Omega.size)
      assert(modes == Lmap.size)
      assert(hc.shape[0] == hc.shape[1])

      if populated_modes is None:
         populated_modesc = modes
      else:
         populated_modesc = populated_modes

      c_init(&tLength, &tSteps, &depth, &modes, &self._dim, &gc[0], &gammac[0],
             &Omegac[0], &hc[0, 0], &Lmapc[0], &with_terminator,
             &populated_modesc)

   @cython.boundscheck(False)
   @cython.wraparound(False)
   def run_trajectory_z0(self, psi0):
      """@todo: Docstring for run_trajectory_z0.

      :psi0: @todo
      :returns: @todo

      """
      cdef double complex[:] psi0c = np.array(psi0, dtype=np.complex128,
                                                 order='F')
      assert(psi0c.size == self._dim)

      cdef np.ndarray[double complex, mode='fortran', ndim=2] psi = \
            np.empty([self._tSteps, self._dim], dtype=np.complex128, order='F')
      c_run_trajectory_z0_rk4(&self._dim, &self._tSteps, &psi0c[0], &psi[0, 0])
      return psi
