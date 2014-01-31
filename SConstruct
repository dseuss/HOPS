import os
from conf import *

mode = ARGUMENTS.get('mode', 'debug')
if not (mode in ['debug', 'release']):
   print('Error: expected "debug" or "release", found: ' + mode)
   Exit(1)

env = Environment(ENV=os.environ,
                  tools=['default', FCC],
                  F90FLAGS=FFLAGS[mode],
                  LINKFLAGS=LDFLAGS[mode],
                  FORTRANMODDIR=FORTRANMODDIR)

env.SConsignFile()
#SConscript('src/SConscript', exports=['env'])
SConscript('test/SConscript', exports=['env'])
