FCC = 'ifort'

FFLAGS = {'debug': ['-g', '-O0'],
          'release': []
          }

LDFLAGS = {'debug': [],
           'release': []
           }

BUILDDIR = '/tmp/FHOPSBUILD/'
FORTRANMODDIR = BUILDDIR + 'mods'

MODULES = ['hstructent.f90',
           'hstructlist.f90',
           'hstructtab.f90']
