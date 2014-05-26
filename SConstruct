# TODO Add options to link statically
import os
import SCons.SConf
import cyscons
from platform import architecture

# Select MKL directory name based on architecture
MKLLIBDIR = {'32bit': 'ia32', '64bit': 'intel64'}
# FLAGS for OpenMP that differ from the usual -fopenmp
OPENMPFLAG = {'ifort': ['-openmp']}
# FLAGS for Intel MKL depending on architecture and compiler
MKLFLAGS = {'32bit': {'gfortran': ['-m32'],
                     'ifort': []},
           '64bit': {'gfortran': ['-m64'],
                     'ifort': []}
           }
MKLLD = {'32bit': {'gfortran': ['mkl_gf', 'mkl_core', 'mkl_gnu_thread', 'dl', 'm'],
                   'ifort': ['mkl_intel', 'mlk_core', 'mkl_intel_thread', 'pthread', 'm']},
         '64bit': {'gfortran': ['mkl_gf_lp64', 'mkl_core', 'mkl_gnu_thread', 'dl', 'pthread', 'm'],
                   'ifort': ['mkl_intel_lp64', 'mkl_core', 'mkl_intel_thread', 'pthread', 'm']}
         }

# Configuration file
SCONS_CONFIG_FILE = 'scons.conf'

# First open the config file
opts = Variables(SCONS_CONFIG_FILE)

# Set all options + default values
opts.AddVariables(('FORTRAN', 'Fortran compiler to use', None),
                  ('FORTRANFLAGS', 'Flags for the Fortran compiler', None),
                  ('CC', 'C compiler to use', None),
                  ('LINKFLAGS', 'Flags for the linker', None),
                  ('MKLROOT', 'Intel MKL Root Directory', None),
                  ('MKLLD', 'Explicit link line for MKL', None),
                  ('MKLFLAGS', 'Explicit compile options for MKL', None),
                  ('ARCH', 'Target architecture', architecture()[0]),
                  )

env = Environment(ENV=os.environ,
                  tools=['default', 'gcc', 'gfortran'],
                  FORTRANFILESUFFIXES=['.f', '.f90', '.F90'],
                  STATIC_AND_SHARED_OBJECTS_ARE_THE_SAME=1,
                  FORTRANMODDIR='${TARGET.dir}',
                  FORTRANMODDIRPREFIX='-J ')

opts.Update(env)

# Custom Stuff for Fortran ####################################################
def CheckFC(context):
    context.Message('Checking whether the Fortran compiler works...')

    fortran_test_file = """
    program test
    print *, 'Fortran compiler operational.'
    end program test"""

    result = context.TryCompile(fortran_test_file, '.f90')
    context.Result(result)
    return result


def CheckFortranInclude(context, filename):
    context.Message('Checking for Fortran header file {}...'.format(filename))

    fortran_test_file = """
    program test
    use, intrinsic :: iso_c_binding
    include '{}'
    end program test""".format(filename)

    result = context.TryCompile(fortran_test_file, '.F90')
    context.Result(result)
    return result


def sortObjectFiles(env, scfile):
    """ Splits scfile into .mod and .o files. The former are moved to the same
        dir as the source files for gfortran (since Syntastic needs them). The
        latter are returned as list.
    """
    allobjs = env.SConscript(scfile, exports=['env'])
    objs = filter(lambda o: str(o)[-4:] != '.mod', allobjs)
    return objs

###############################################################################


def ConfigureCompilers(conf):
    compiler = conf.env['FORTRAN']

    if not all((conf.CheckCC(), conf.CheckFC())):
        print 'ERROR: Compiler configuration failed, exiting.'
        Exit(1)
    print 'Using Fortran compiler: ' + conf.env['FORTRAN']
    print 'Using C compiler: ' + conf.env['CC']
    conf.env.Append(FORTRANFLAGS=[OPENMPFLAG.get(compiler, '-fopenmp')],
                    _LIBFLAGS=[OPENMPFLAG.get(compiler, '-fopenmp')])


def ConfigureCython(conf):
    CYHEADER = ('pyconfig.h', 'Python.h', 'numpy/arrayobject.h', 'numpy/ufuncobject.h')
    if not all((conf.CheckHeader(header) for header in CYHEADER)):
        print 'ERROR: CYTHON configuration failed, exiting.'
        # Exit(1)


def ConfigureFFTW(conf):
    if not all((conf.CheckLib('fftw3'), conf.CheckFortranInclude('fftw3.f03'))):
        print 'ERROR: FFTW3 configuration failed, exiting.'
        Exit(1)


def ConfigureMKL(conf):
    arch = conf.env['ARCH']
    compiler = conf.env['FORTRAN']

    if (not conf.env.has_key('MKLROOT')) and (os.environ.has_key('MKLROOT')):
        conf.env['MKLROOT'] = os.environ['MKLROOT']

    print 'Checking MKL root directory...',
    print 'yes' if conf.env.has_key('MKLROOT') else 'no'

    # Use the compiler options provided -- or use sane defaults
    if conf.env.has_key('MKLFLAGS'):
        print 'Setting MKL compiler options manually: ' + conf.env['MKLFLAGS']
        conf.env.Append(FORTRANFLAGS=str(conf.env['MKLFLAGS']).split())
    else:
        if conf.env.has_key('MKLROOT'):
            conf.env.Append(FORTRANPATH=['$MKLROOT/include'])
        conf.env.Append(FORTRANFLAGS=MKLFLAGS[arch][compiler])

    ok = conf.CheckFortranInclude('mkl_spblas.fi')

    # Use the linker options provided -- or use sane defaults
    if conf.env.has_key('MKLLD'):
        print 'Setting MKL linker options manually: ' + conf.env['MKLLD']
        conf.env.Append(_LIBFLAGS=str(conf.env['MKLLD']).split())
    else:
        if conf.env.has_key('MKLROOT'):
            conf.env.Append(LIBPATH=['$MKLROOT/lib/{}'.format(MKLLIBDIR[arch])])
        ok = ok and all((conf.CheckLib(lib) for lib in MKLLD[arch][compiler]))

    if not ok:
        print 'ERROR: MKL configuration failed, exiting.'
        Exit(1)


if not any((env.GetOption(arg) for arg in ('clean', 'help'))):
    ## TODO Make this more configurable (put in config file)
    env.Append(CCFLAGS='-pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -fPIC'.split(),
               CPPPATH=['/usr/include/python2.7',
                        '/usr/local/lib64/python2.7/site-packages/numpy/core/include',
                        '/usr/lib/python2.7/dist-packages/numpy/core/include/'],
               FORTRANPATH=['/usr/include'],
               LIBPATH=['/usr/lib/x86_64-linux-gnu/'],
               FORTRANFLAGS=['-fPIC'])

    cyscons.generate(env)

    # Custom configuration
    conf = env.Configure(custom_tests={'CheckFC': CheckFC,
                                       'CheckFortranInclude': CheckFortranInclude})
    ConfigureCompilers(conf)
    ConfigureCython(conf)
    ConfigureFFTW(conf)
    ConfigureMKL(conf)
    env = conf.Finish()

    env.Append(FORTRANPATH=['#/src', '#src/zvode'])

    # FIXME MODDIR and MODDIRPREFIX for other compilers
    if env['FORTRAN'] == 'ifort':
        print 'Remark: Fortran mod files are moved to ./build'
        env['FORTRANMODDIRPREFIX'] = '-module '
        env['FORTRANMODDIR'] = '#build'

        # This doesnt work with FORTRANPATH
        env.Append(FORTRANFLAGS=['-diag-disable', '7713'],
                   F90PATH=['#build'])


objs = sortObjectFiles(env, r'#src/SConscript')
testobjs = sortObjectFiles(env, r'#test/SConscript')
Help(opts.GenerateHelpText(env))

runtest = env.Program('runtest', objs + testobjs)
libhierarchy = env.SharedLibrary('libhierarchy.so', ['src/libhierarchy.pyx'] + objs)
libnoise = env.SharedLibrary('libnoise.so', ['src/libnoise.pyx'] + objs)
main = env.Program('main', objs + ['src/main.f90'])
Default(runtest)
