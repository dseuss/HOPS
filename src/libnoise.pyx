import numpy as np
cimport numpy as np
import cython

cdef extern:
   void c_test(double *dt, int *tSteps, int *modes, double complex *g,
               double *gamma, double *Omega, int *realizations,
               double complex *EZ, double complex *EZZ, double complex *EZccZ)


def test(double dt, int tSteps, double complex[:] g, double[:] gamma,
         double[:] Omega, int realizations):
   cdef np.ndarray[double complex, mode='fortran'] EZ = \
         np.empty([tSteps], dtype=np.complex128, order='F')
   cdef np.ndarray[double complex, mode='fortran'] EZZ = \
         np.empty([tSteps], dtype=np.complex128, order='F')
   cdef np.ndarray[double complex, mode='fortran'] EZccZ = \
         np.empty([tSteps], dtype=np.complex128, order='F')
   cdef int modes = len(g)

   c_test(&dt, &tSteps, &modes, &g[0], &gamma[0], &Omega[0], &realizations,
          &EZ[0], &EZZ[0], &EZccZ[0])
   t = np.arange(tSteps)[:, None] * dt

   gp = np.array(g)
   gammap = np.array(gamma)
   Omegap = np.array(Omega)
   alpha = np.sum(gp[None, :] * np.exp((-gammap - 1.j*Omegap)[None, :] * t),
                  axis=1)

   return EZ, EZZ, EZccZ, alpha
