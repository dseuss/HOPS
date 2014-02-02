#!/usr/bin/env python
# encoding: utf-8

from __future__ import division, print_function
import os
import platform

def moveModFiles(target=None, source=None, env=None):
    import glob, os, os.path
    targetdir = target[0].dir
    for t in target:
        if t.name[-4:] == '.mod':
            os.rename(t.name,os.path.join(str(targetdir),t.name))

def sortObjectFiles(scfile):
   allobjs = env.SConscript(scfile,exports=['env'])
   objs = filter(lambda o: str(o)[-4:] != '.mod', allobjs)
   if FCC == 'gfortran':
      env.AddPostAction(objs, moveModFiles)
   return objs

MKLROOT = os.environ['MKLROOT']
MKLFLAGS = {
'gfortran32': '-m32 -I{MKLROOT}/include'.format(**{ 'MKLROOT': MKLROOT }).split(),
'gfortran64': '-m64 -I{MKLROOT}/include'.format(**{ 'MKLROOT': MKLROOT }).split(),
'ifort32': '-I{MKLROOT}/include'.format(**{ 'MKLROOT': MKLROOT }).split(),
'ifort64': '-I{MKLROOT}/include'.format(**{ 'MKLROOT': MKLROOT }).split()
            }
MKLLD = {
'gfortran32': ' -Wl,--start-group {MKLROOT}/lib/ia32/libmkl_gf.a {MKLROOT}/lib/ia32/libmkl_core.a {MKLROOT}/lib/ia32/libmkl_sequential.a -Wl,--end-group -lpthread -lm'.format(**{ 'MKLROOT': MKLROOT }).split(),
'ifort32': ' -Wl,--start-group {MKLROOT}/lib/ia32/libmkl_intel.a {MKLROOT}/lib/ia32/libmkl_core.a {MKLROOT}/lib/ia32/libmkl_sequential.a -Wl,--end-group -lpthread -lm'.format(**{ 'MKLROOT': MKLROOT }).split(),
'gfortran64': ' -Wl,--start-group {MKLROOT}/lib/intel64/libmkl_gf_lp64.a {MKLROOT}/lib/intel64/libmkl_core.a {MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm'.format(**{ 'MKLROOT': MKLROOT }).split(),
'ifort64': ' -Wl,--start-group {MKLROOT}/lib/intel64/libmkl_intel_lp64.a {MKLROOT}/lib/intel64/libmkl_core.a {MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm'.format(**{ 'MKLROOT': MKLROOT }).split()
         }


mode = ARGUMENTS.get('mode', 'release')
if not (mode in ['debug', 'release']):
   print('ERROR: expected "debug" or "release", found ' + mode)
   Exit(1)

FCC = ARGUMENTS.get('FCC', 'gfortran')
SETTINGS = FCC + {'i686': '32', 'x86_64': '64'}.get(platform.uname()[-2], '32')
if not (FCC in ['ifort', 'gfortran']):
   print('ERROR: expected "ifort" or "gfortran", found ' + FCC)
   Exit(1)

FLAGS = {'release': '-O3'.split(),
         'debug': '-g -O0'.split()
         }
F2PYFLAGS = {'ifort': '-lifcore -lifport'.split(),
             'gfortran': '-lgfortran'.split()}

env = Environment(ENV=os.environ,
                  tools=['default', FCC],
                  F90FLAGS=FLAGS[mode] + MKLFLAGS[SETTINGS],
                  _LIBFLAGS=MKLLD[SETTINGS],
                  F90PATH='#/src /usr/include'.split(),
                  )

if FCC == 'ifort':
   env.Append(F90PATH='#build',
              FORTRANMODDIR='#build')

objs = sortObjectFiles(r'#src/SConscript')
testobjs = sortObjectFiles(r'#test/SConscript')

runtest = env.Program('runtest', objs + testobjs,
                      LIBS=['fftw3'],
                      LIBPATH=['/usr/lib/i386-linux-gnu/'])

# f2py = env.Clone()
# f2py.Append(F90FLAGS=env['F90FLAGS'] + F2PYFLAGS[FCC])
# libhierarchy = f2py.Program('libhierarchy.so', objs,
                            # LINKCOM='f2py -c --fcompiler=gnu95 src/libhierarchy.pyf $SOURCES --f90flags="$F90FLAGS" -I#/src')
# Default(libhierarchy)
