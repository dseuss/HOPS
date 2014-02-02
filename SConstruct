#!/usr/bin/env python
# encoding: utf-8

from __future__ import division, print_function
import os

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

mode = ARGUMENTS.get('mode', 'release')
if not (mode in ['debug', 'release']):
   print('ERROR: expected "debug" or "release", found ' + mode)
   Exit(1)

FCC = 'ifort'
FLAGS = {'release': '-O3'.split(),
         'debug': '-g -O0'.split()
         }

env = Environment(ENV=os.environ,
                  tools=['default', FCC],
                  F90FLAGS=FLAGS[mode],
                  F90PATH='#/src /usr/include'.split(),
                  )

if FCC == 'ifort':
   env.Append(F90PATH='#build',
              FORTRANMODDIR='#build')

objs = sortObjectFiles(r'#src/SConscript')
testobjs = sortObjectFiles(r'#test/SConscript')

env.Program('runtest', objs + testobjs,
            LIBS='fftw3', LIBPATH=['/usr/lib/i386-linux-gnu/'])
