#!/usr/bin/env python
# encoding: utf-8

from __future__ import division, print_function
import sys
import numpy as np
import h5py
from time import strftime, time
import multiprocessing as mp

import functions_ger as fg
from libbath import OscillatorBath
from libhierarchy import FHierarchy, set_omp_threads

from mpi4py import MPI
comm = MPI.COMM_WORLD
mpi_size = comm.Get_size()
mpi_rank = comm.Get_rank()

OUTPUT_DELAY = 10


###############################################################################
#                               class Hierarchy                               #
###############################################################################
class Hierarchy(object):

   """Docstring for Hierarchy. """

   def __init__(self, tLength, tSteps, hamiltonian, bath, depth, filename):
      """Initialize Hierarchy

      :tLength:
      :tSteps:
      :hamiltonian: @todo
      :bath: @todo
      :depth: @todo

      """
      self._hamiltonian = np.array(hamiltonian)
      self._dimension = self._hamiltonian.shape[0]
      self._bath = bath
      self._depth = depth
      self._tLength = tLength
      self._tSteps = tSteps
      self._filename = filename

      assert(self._hamiltonian.shape[0] == self._hamiltonian.shape[1])
      assert(depth >= 0)

      self._integrator = None
      self.update_integrator()

   def update_integrator(self):
      """@todo: Docstring for update_integrator.
      :returns: @todo

      """
      if self._integrator is not None:
         del self._integrator

      g = np.ones([self._dimension])[:, None] * self._bath.g[None, :]
      gamma = np.ones([self._dimension])[:, None] * self._bath.gamma[None, :]
      Omega = np.ones([self._dimension])[:, None] * self._bath.Omega[None, :]
      Lmap = np.arange(1, self._dimension + 1)[:, None] * \
            np.ones([self._bath.nrmodes], dtype=int)[None, :]
      set_omp_threads(-1)
      # FIXME Switch on terminator
      self._integrator = FHierarchy(self._tLength, self._tSteps, self._depth,
                                    g.flatten(), gamma.flatten(),
                                    Omega.flatten(), self._hamiltonian,
                                    Lmap.flatten(), with_terminator=False)

   def _set_attrs(self, ds):
      """@todo: Docstring for _set_attrs.

      :ds: @todo
      :returns: @todo

      """
      ds.attrs.create('Hamiltonian', self._hamiltonian)
      ds.attrs.create('g', self._bath.g)
      ds.attrs.create('gamma', self._bath.gamma)
      ds.attrs.create('Omega', self._bath.Omega)
      ds.attrs.create('Depth', self._depth)
      ds.attrs.create('tLength', self._tLength)
      ds.attrs.create('Monomers', self._dimension)

   def _get_attrs(self, ds):
      """@todo: Docstring for _get_attrs.

      :ds: @todo
      :returns: @todo

      """
      self._hamiltonian = ds.attrs['Hamiltonian']
      self._depth = ds.attrs['Depth']
      self._tLength = ds.attrs['tLength']
      self._bath = OscillatorBath(ds.attrs['g'], ds.attrs['gamma'],
                                  ds.attrs['Omega'])
      self._dimension = ds.attrs['Monomers']


###############################################################################
#                           class SpectrumHierarchy                           #
###############################################################################
class SpectrumHierarchy(Hierarchy):

   """Docstring for SpectrumHierarchy. """

   def __init__(self, tLength, tSteps, hamiltonian, bath, depth,
                filename='spectra.h5'):
      """@todo: to be defined1.

      :tLength: @todo
      :tSteps: @todo
      :hamiltonian: @todo
      :bath: @todo
      :depth
      :filename: @todo

      """
      Hierarchy.__init__(self, tLength, tSteps, hamiltonian, bath, depth,
                         filename)
      self._psi = None

   def calc_trajectory(self, psi0=None, label=None):
      """Calculate linear trajectory

      :tLength: @todo
      :tSteps: @todo
      :psi0: @todo
      :returns: @todo

      """
      # TODO Differntiate between noise/no noise
      if psi0 is None:
         psi_init = np.ones(self._dimension) / np.sqrt(self._dimension)
      else:
         psi_init = psi0

      self._psi = self._integrator.run_trajectory_z0(psi_init)

      if label is not None:
         self.to_file(label=label)

      return self.get()

   def get(self, wmin=None, wmax=None, sigma_w=0.):
      """@todo: Docstring for get.

      :wmin: @todo
      :wmax: @todo
      :sigma_w: @todo
      :returns: @todo

      """
      assert(self._psi is not None)

      dt = self._tLength / (self._psi.shape[0] - 1)
      corr = np.sum(self._psi * np.conj(self._psi[0]), axis=1)
      w, A = fg.fourier(corr, dt, output_w=True, hermitian=True,
                        sigma_w=sigma_w)

      ## Cut out values which are out of bound
      sel = np.ones(w.size, dtype=bool)
      if wmin is not None:
         sel *= w >= wmin
      if wmax is not None:
         sel *= w <= wmax

      return w[sel], np.real(A[sel])

   def to_file(self, label=None):
      """@todo: Docstring for to_file.

      :label: @todo
      :returns: @todo

      """
      assert(self._psi is not None)

      if label is None:
         lab = strftime('%Y-%m-%d_%H:%M')
      else:
         lab = label

      with h5py.File(self._filename, 'a') as ofile:
         if lab in ofile:
            del ofile[lab]
         ds = ofile.create_dataset(lab, data=self._psi)
         self._set_attrs(ds)
         ds.attrs.create('type', 'spectrum_fixed')

   def from_file(self, filename, label):
      """@todo: Docstring for from_file.

      :filename: @todo
      :label: @todo
      :returns: @todo

      """
      self._filename = filename
      with h5py.File(self._filename, 'r') as ofile:
         ds = ofile[label]
         assert(ds.attrs['type'] == 'spectrum_fixed')
         self._get_attrs(ds)
         self._psi = ds.value


###############################################################################
#                           class TransferHierarchy                           #
###############################################################################
# class TransferHierarchy(Hierarchy):

#    """Docstring for TransferHierarchy. """

#    def __init__(self, tLength, tSteps, hamiltonian, bath, depth,
#                 filename='transfer.h5'):
#       """@todo: to be defined1.

#       :tLength: @todo
#       :tSteps: @todo
#       :hamiltonian: @todo
#       :bath: @todo
#       :depth
#       :filename: @todo

#       """
#       Hierarchy.__init__(self, tLength, tSteps, hamiltonian, bath, depth,
#                          filename)
#       self._rho = None

#    def _get_rho_partial(self, realizations, psi0, label, id, omp_threads,
#                         result):
#       """Calculate the reduced density operator for a given number of
#       realizations in a single thread. THIS IS WHERE THE MAIN WORK IS DONE!

#       :realizations: @todo
#       :psi0: @todo
#       :label: @todo
#       :returns: @todo

#       """
#       print('Thread nr. {} [{}] with {} realizations'.format(id, mpi_rank,
#                                                              realizations))
#       rho = np.zeros((self._tSteps, self._dimension, self._dimension),
#                      dtype=complex)
#       lh.set_omp_threads(omp_threads)
#       t_output = time()
#       for i in xrange(realizations):
#          # FIXME Better status update
#          if time() - t_output > OUTPUT_DELAY:
#             print('{} [{}] -- {}/{}'.format(id, mpi_rank, i, realizations))
#             t_output = time()
#          psi = self._integrator.run_trajectory(psi0)
#          norm = np.sum(np.conj(psi) * psi, axis=1)[:, None, None]
#          rho += np.conj(psi[:, None, :]) * psi[:, :, None] / norm
#       result.put(rho)

#    def _get_rho_total(self, realizations, psi0, num_threads, omp_threads,
#                       label):
#       """Calculate the reduced density operator on the local machine by
#       distributing work over multiple threads for calculating independent
#       realizations.

#       :realizations: @todo
#       :psi0: @todo
#       :num_threads: @todo
#       :label: @todo
#       :returns: @todo

#       """
#       traj_id = np.array_split(np.arange(realizations), num_threads)
#       results = mp.Queue()
#       processes = [mp.Process(target=self._get_rho_partial,
#                               args=(traj_id[i].size, psi0, label, i,
#                                     omp_threads, results))
#                    for i in range(num_threads)]

#       for p in processes:
#          p.start()
#       rho = [results.get() for p in processes]
#       for p in processes:
#          p.join()

#       return np.sum(rho, axis=0)

#    # FIXME omp_threads isnt doing anything!
#    def calc(self, realizations, psi0, num_threads=-1, omp_threads=1,
#             label=None, normalize=True):
#       """Calculate reduced density operator for given number of realizations.
#       This function manly distributes and gathers all calculation over multiple
#       MPI-nodes.

#       :realizations:
#       :psi0:
#       :num_threads: Number of threads used on the local machine for calculating
#                     independent run_trajectories.
#       :omp_threads: Number of threads used for parallel calculation of matrix-
#                     multiplications.
#       :label:
#       :normalize: Normalize the final result to Tr ρ = 1

#       """
#       realizations_local = int(realizations/mpi_size)
#       if mpi_rank < realizations - realizations_local*mpi_size:
#          realizations_local += 1

#       rho = self._get_rho_total(realizations_local, psi0, num_threads,
#                                 omp_threads, label)

#       rho = comm.gather(rho, root=0)
#       if mpi_rank == 0:
#          # TODO Add save to file
#          t = np.linspace(0, self._tLength, self._tSteps)
#          self._rho = np.sum(rho, axis=0)
#          if label is not None:
#             self.to_file(label)
#          return t, self._get(normalize)
#       else:
#          return None, None

#    def _get(self, normalize=True):
#       """Returns the calculated density operator.

#       :normalize: @todo
#       :returns: @todo

#       """
#       assert(self._rho is not None)
#       if normalize:
#          return self._rho /\
#                np.trace(self._rho, axis1=1, axis2=2)[:, None, None]
#       else:
#          return self._rho

#    def to_file(self, label=None):
#       """@todo: Docstring for to_file.

#       :label: @todo
#       :returns: @todo

#       """
#       if mpi_rank != 0:
#          return

#       assert(self._rho is not None)

#       if label is None:
#          lab = strftime('%Y-%m-%d_%H:%M')
#       else:
#          lab = label

#       with h5py.File(self._filename, 'a') as ofile:
#          if lab in ofile:
#             del ofile[lab]
#          ds = ofile.create_dataset(lab, data=self._rho)
#          self._set_attrs(ds)
#          ds.attrs.create('type', 'transfer_fixed')

#    def from_file(self, filename, label):
#       """@todo: Docstring for from_file.

#       :filename: @todo
#       :label: @todo
#       :returns: @todo

#       """
#       self._filename = filename
#       with h5py.File(self._filename, 'r') as ofile:
#          ds = ofile[label]
#          assert(ds.attrs['type'] == 'transfer_fixed')
#          self._get_attrs(ds)
#          self._rho = ds.value


# ###############################################################################
# #                              Helper Functions                               #
# ###############################################################################
def LoadSpectrumFromFile(filename, label):
   """@todo: Docstring for LoadSpectrumFromFile.

   :filename: @todo
   :label: @todo
   :returns: @todo

   """
   hier = SpectrumHierarchy(0, 0, [[0]], OscillatorBath(0, 0, 0), 0)
   hier.from_file(filename, label)
   return hier


def StandardHamiltonian(Monomers, site_energies, couplings):
   """Create a simple aggregate Hamiltonian with site_energies in diagonal
   and couplings in off-diagonals with size Monomers x Monomers

   :Monomers: @todo
   :site_energies: @todo
   :couplings: @todo
   :returns: @todo

   """
   H = site_energies * np.identity(Monomers, dtype=complex)
   H += couplings * (np.identity(Monomers + 1, dtype=complex)[1:, :-1]
                    + np.identity(Monomers + 1, dtype=complex)[:-1, 1:])
   return H
