import numpy as np
cimport numpy as np
cimport openmp
import cython
from cpython cimport bool
from warnings import warn


###############################################################################
#                       Wrapper for Hierachy Integrator                       #
###############################################################################
cdef extern:
    void c_init(double *tLength, int *tSteps, int *depth, int *modes,
               int *hs_dim, double complex *g, double *gamma, double *Omega,
               double complex *h, int *Lmap, bint *with_terminator,
               int *populated_modes)
    void c_run_trajectory_rk4(int *hs_dim, int *tSteps, complex *psi0,
                             double complex *psi, bint *normalized)
    void c_run_trajectory_z0_rk4(int *hs_dim, int *tSteps, complex *psi0,
                                double complex *psi)
    void c_run_trajectory_z0_zvode(int *hs_dim, int *tSteps, complex *psi0,
                                  double complex *psi, int *method,
                                  double *rtol, double *atol)
    void c_trajectory_step_z0(int *NEQ, double *T, double complex *psi,
                                 double complex *psi_dot)
    void c_get_size(int *size)
    void c_free()


cdef class _HierarchyIntegrator(object):

    """ Interface to the Fortran hierarchy integrator.
    IMPORTANT: Do not create an instance of this class manually, use the
               fIntegrator provided. Since the underlying Fortran object is
               a module, at most one instance at every time is allowed.
    """

    cdef int _tSteps
    cdef double _tLength
    cdef int _dim
    cdef int _size

    cdef bool _initialized

    def __init__(self):
        """
        The underlying Fortran integrator is not initialized by default, call
        the update function before calculations.
        """
        self._initialized = False

    def update(self, double tLength, int tSteps, int depth, g, gamma, Omega,
               h, Lmap, bint with_terminator=1, int populated_modes=0):
        """ Initializes the underlying Fortran integrator.
        This overwrites prior set parameters -- integrator is not reentrant.
        The bath correlation function is parametrized according to

            alpha(t) = sum_i g_i * exp(-gamma_i * |t| - ii * Omega_i * t).

        :tLength: Propagation time
        :tSteps: Number of time steps
        :depth: Depth of the hierarchy
        :g[modes]: Array of coupling strengths in bcf
        :gamma[modes]: Array of dampings in bcf
        :Omega[modes]: Array of Frequencies in bcf
        :h[dim, dim]: System hamiltonian
        :Lmap[modes]: Maps the individual modes to the corresponding coupling
                      operator, i.e. i-th mode couples with L_Lmap[i],
                      where L_j = |pi_j><pi_j|.
        :with_terminator(true): Use the standard terminator?
        :populated_modes(0): Discard all auxiliary states with excitations
                             on more than populated_modes modes. If this is
                             smaller than 1, no states are discarded.

        """
        if self._initialized:
            warn("WARNING: Overwriting old hierarchy-integrator. Libhierachy is not reentrant.")
            c_free()

        # Convert bcf parameters to numpy arrays, so that we can also pass lists
        cdef double complex[:] gc = np.array(g, dtype=np.complex128, order='F')
        cdef double[:] gammac = np.array(gamma, dtype=np.float64, order='F')
        cdef double[:] Omegac = np.array(Omega, dtype=np.float64, order='F')
        cdef double complex[:, :] hc = np.array(h, dtype=np.complex128,
                                                order='F')
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

        if populated_modes < 1:
            populated_modesc = modes
        else:
            populated_modesc = populated_modes

        c_init(&tLength, &tSteps, &depth, &modes, &self._dim, &gc[0],
               &gammac[0], &Omegac[0], &hc[0, 0], &Lmapc[0], &with_terminator,
               &populated_modesc)
        c_get_size(&self._size)
        self._initialized = True

    def __dealloc__(self):
        if self._initialized:
            c_free()

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def run_trajectory(self, psi0, bint normalized=0):
        """ Calculates a normalizable quantum trajectory.
        Initialize the integrator first using update(). The solutions are not
        normalized, but allow for a normalization when taking the average
        (non-normalized, but Girsanov-shifted trajectories).

        :psi0[dim]: Initial state
        :normalized(False): Return normalized trajectory
        :returns: psi[tSteps, dim], quantum trajectory

        """
        cdef double complex[:] psi0c = np.array(psi0, dtype=np.complex128,
                                                order='F')
        cdef np.ndarray[double complex, mode='fortran', ndim=2] psi = \
                np.empty([self._tSteps, self._dim], dtype=np.complex128,
                         order='F')
        c_run_trajectory_rk4(&self._dim, &self._tSteps, &psi0c[0], &psi[0, 0],
                             &normalized)
        return psi

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def run_trajectory_z0(self, psi0, integrator='rk4', **kwargs):
        """ Calculates the solution psi(Z=0) of the linear hierarchy

        :psi0[dim]: Initial state
        :integrator: Determines the integrator used, possible values: rk4, zvode
        :returns: psi[tSteps, dim], quantum trajectory

        Additional keyword arguments for integrator == zvode:
        :method: Determines the algorithm, possible values: adams, bcf
        :rtol: Relative tolerance
        :atol: Absolute tolerance

        """
        cdef double complex[:] psi0c = np.array(psi0, dtype=np.complex128,
                                                order='F')
        # Determine the arguments for zvode
        method_tag = kwargs.get('method', 'adams')
        method_flags = {'adams': 10, 'bcf': 22}
        cdef int method = method_flags[method_tag]
        cdef double rtol = kwargs.get('rtol', 10e-6)
        cdef double atol = kwargs.get('atol', 10e-6)

        assert(psi0c.size == self._dim)
        cdef np.ndarray[double complex, mode='fortran', ndim=2] psi = \
                np.empty([self._tSteps, self._dim], dtype=np.complex128,
                         order='F')

        if integrator == 'zvode':
            c_run_trajectory_z0_zvode(&self._dim, &self._tSteps, &psi0c[0],
             &psi[0, 0], &method, &rtol, &atol)
        else:
            c_run_trajectory_z0_rk4(&self._dim, &self._tSteps, &psi0c[0],
             &psi[0, 0])

        return psi

    def trajectory_step_z0(self, double t, double complex[:] psi):
        """ Integration step for run_trajectory_z0.
        For use with external integrator, provides the internal integration step
        used in the Fortran integrator.

        :t: Time (irrelevant since problem is time-independent)
        :psi: Full hierarchy state
        :returns: psi_dot, the lhs of the evolution equation

        """
        cdef double complex[:] psi_dot = \
                np.empty([self._size], dtype=np.complex128, order='F')
        c_trajectory_step_z0(&self._size, &t, &psi[0], &psi_dot[0])
        return psi_dot


fIntegrator = _HierarchyIntegrator()

###############################################################################
#                    Wrapper for further helper functions                     #
###############################################################################
def set_omp_threads(int omp_threads):
    """ Set the number of threads used by OpenMP

    :omp_threads: Number of threads, if omp_threads < 1 the maximum number of
                  threads is used.

    """
    if omp_threads < 1:
        openmp.omp_set_num_threads(openmp.omp_get_max_threads())
    else:
        openmp.omp_set_num_threads(omp_threads)


cdef extern:
    void c_init_random_seed()
def init_random_seed():
    """ Sets the random seed to a random value. """
    c_init_random_seed()

