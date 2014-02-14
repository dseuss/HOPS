import numpy as np
cimport numpy as np
cimport openmp
import cython

cdef extern:
   void c_init(double *tLength, int *tSteps, int *depth, int *modes,
               int *hs_dim, double complex *g, double *gamma, double *Omega,
               double complex *h, int *Lmap, bint *with_terminator,
               int *populated_modes)
   void c_run_trajectory_z0_rk4(int *hs_dim, int *tSteps, complex *psi0,
                                double complex *psi)
   void c_run_trajectory_z0_zvode(int *hs_dim, int *tSteps, complex *psi0,
                                  double complex *psi, int *method,
                                  double *rtol, double *atol)
   void c_trajectory_step_z0(int *NEQ, double *T, double complex *psi,
                                 double complex *psi_dot)
   void c_get_size(int *size)
   void c_free()


cdef class FHierarchy(object):

   """Docstring for FHierarchy. """

   cdef int _tSteps
   cdef double _tLength
   cdef int _dim
   cdef int _size

   @property
   def size(self):
      return self._size

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
      self._tLength = tLength
      cdef int populated_modesc
      cdef int modes = gc.size
      assert(modes == gammac.size)
      assert(modes == Omegac.size)
      assert(modes == Lmapc.size)
      assert(hc.shape[0] == hc.shape[1])

      if populated_modes is None:
         populated_modesc = modes
      else:
         populated_modesc = populated_modes

      c_init(&tLength, &tSteps, &depth, &modes, &self._dim, &gc[0], &gammac[0],
             &Omegac[0], &hc[0, 0], &Lmapc[0], &with_terminator,
             &populated_modesc)
      c_get_size(&self._size)

   def __dealloc__(self):
      c_free()

   @cython.boundscheck(False)
   @cython.wraparound(False)
   def run_trajectory_z0(self, psi0, integrator='rk4', **kwargs):
      """@todo: Docstring for run_trajectory_z0.

      :psi0: @todo
      :returns: @todo

      """
      cdef double complex[:] psi0c = np.array(psi0, dtype=np.complex128,
                                                 order='F')
      method_tag = kwargs.get('method', 'adams')
      method_flags = {'adams': 10, 'bcf': 22}
      cdef int method = method_flags[method_tag]
      cdef double rtol = kwargs.get('rtol', 10e-6)
      cdef double atol = kwargs.get('atol', 10e-6)

      assert(psi0c.size == self._dim)
      cdef np.ndarray[double complex, mode='fortran', ndim=2] psi = \
            np.empty([self._tSteps, self._dim], dtype=np.complex128, order='F')

      if integrator == 'zvode':
         c_run_trajectory_z0_zvode(&self._dim, &self._tSteps, &psi0c[0],
                                 &psi[0, 0], &method, &rtol, &atol)
      else:
         c_run_trajectory_z0_rk4(&self._dim, &self._tSteps, &psi0c[0],
                                 &psi[0, 0])
      return psi

   def trajectory_step_z0(self, double t, double complex[:] psi):
      """@todo: Docstring for trajectoy_step_z0.

      :t: @todo
      :psi: @todo
      :returns: @todo

      """
      cdef double complex[:] psi_dot = \
            np.empty([self._size], dtype=np.complex128, order='F')
      c_trajectory_step_z0(&self._size, &t, &psi[0], &psi_dot[0])
      return psi_dot

###############################################################################
#                    Wrapper for further helper functions                     #
###############################################################################
def set_omp_threads(int omp_threads):
   if omp_threads < 0:
      openmp.omp_set_num_threads(openmp.omp_get_max_threads())
   else:
      openmp.omp_set_num_threads(omp_threads)
