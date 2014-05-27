import numpy as np
cimport numpy as np
import cython

cdef extern:
   void c_test(double *dt, int *tSteps, int *modes, double complex *g,
               double *gamma, double *Omega, int *realizations,
               double complex *EZ, double complex *EZZ, double complex *EZccZ)
   void c_init_random_seed()


def test(double dt, int tSteps, double complex[:] g, double[:] gamma,
         double[:] Omega, int realizations):
    """
    Test the noise module by calculation EZ, EZZ, and EZccZ assuming a
    bath correlation function of the form

        alpha(t) = sum_j g_j * exp(-gamma_j * |t| - ii * Omega_j * t)

    :dt: Size of one time interval
    :tSteps: Number of time steps
    :g[modes]: Array of coupling strengths in bcf
    :gamma[modes]: Array of dampings in bcf
    :Omega[modes]: Array of Frequencies in bcf
    :realizations: Number of realizations averaged over
    :returns: EZ_t (mean), E(Z_t * Z_0) (covariance), E(Z_t cc(Z_0) (complex
              conjugate covariance), alpha (correlation function as it should
              be)

    """

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


def init_random_seed():
    """ Initializes the random seed with a random seed. """
    # TODO Make random seed an (optional) argument
    c_init_random_seed()
