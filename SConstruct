#!/usr/bin/env python
# encoding: utf-8

from __future__ import division, print_function
import os
from platform import architecture
import cyscons

def moveModFiles(target=None, source=None, env=None):
   """ Moves .mod files to the same dir as the source files are in. """
   targetdir = target[0].dir
   for t in target:
      if t.name[-4:] == '.mod':
         os.rename(t.name, os.path.join(str(targetdir), t.name))


def sortObjectFiles(scfile):
   """ Splits scfile into .mod and .o files. The former are moved to the same
       dir as the source files for gfortran (since Syntastic needs them). The
       latter are returned as list.
   """
   allobjs = env.SConscript(scfile, exports=['env'])
   objs = filter(lambda o: str(o)[-4:] != '.mod', allobjs)
   if FCC == 'gfortran':
      env.AddPostAction(objs, moveModFiles)
   return objs

# TODO This has terrible form!!!!!

mode = ARGUMENTS.get('mode', 'debug')
if not (mode in ['debug', 'release']):
   print('ERROR: expected "debug" or "release", found ' + mode)
   Exit(1)

FCC = ARGUMENTS.get('FCC', 'gfortran')
if not (FCC in ['ifort', 'gfortran']):
   print('ERROR: expected "ifort" or "gfortran", found ' + FCC)
   Exit(1)

FLAGS = {'release': '-O3'.split(),
         'debug': '-g -O0 -fPIC'.split()
         }
LDFLAGS = {'ifort': '-openmp', 'gfortran': '-fopenmp'}[FCC].split()

MKLROOT = os.environ['MKLROOT']
ARCHBITS = architecture()[0][:2]
GFORTRAN = {'gfortran': '-m{}'.format(ARCHBITS), 'ifort': ''}[FCC]
ARCHDIR = {'32': 'ia32', '64': 'intel64'}[ARCHBITS]
IDENTIFIER = {'ifort32': 'intel', 'ifort64': 'intel_lp64',
              'gfortran32': 'gf', 'gfortran64': 'gf_lp64'}[FCC + ARCHBITS]
THREADEDSWITCH = 'sequential'
# THREADEDSWITCH = 'gnu_thread'
# THREADEDSWITCH = 'intel_thread'

MKLLD = ' -Wl,--start-group {MKLROOT}/lib/{ARCHDIR}/libmkl_{IDENTIFIER}.a {MKLROOT}/lib/{ARCHDIR}/libmkl_core.a {MKLROOT}/lib/{ARCHDIR}/libmkl_{THREADEDSWITCH}.a -Wl,--end-group -lpthread -lm'.format(**locals()).split()
MKLFLAGS = '{GFORTRAN} -I{MKLROOT}/include'.format(**locals()).split()

env = Environment(ENV=os.environ,
                  tools=['default', FCC],
                  CCFLAGS=' -fPIC -O3 -I/usr/include/python2.7',
                  F90FLAGS=FLAGS[mode] + MKLFLAGS,
                  _LIBFLAGS=MKLLD + LDFLAGS,
                  F90PATH='#/src #src/zvode /usr/include'.split(),
                  STATIC_AND_SHARED_OBJECTS_ARE_THE_SAME=1
                  )
cyscons.generate(env)


if FCC == 'ifort':
   env.Append(F90PATH='#build',
              FORTRANMODDIR='#build')

objs = sortObjectFiles(r'#src/SConscript')
testobjs = sortObjectFiles(r'#test/SConscript')

runtest = env.Program('runtest', objs + testobjs,
                      LIBS=['fftw3'],
                      LIBPATH=['/usr/lib/i386-linux-gnu/'])
libhierarchy = env.SharedLibrary('libhierarchy.so',
                                 ['src/libhierarchy.pyx'] + objs)
libnoise = env.SharedLibrary('libnoise.so',
                                 ['src/libnoise.pyx'] + objs)
main = env.Program('main', objs + ['src/main.f90'])
Default(libhierarchy)
