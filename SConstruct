#!python

import os
import platform

# Configure build environment.
env = None
if (platform.system() == 'Windows'):

    #
    # Check for user defined compiler.
    # i.e. > scons.bat -Q compiler=mingw
    #
    # The compiler needs to be in the users path.
    #
    compiler = ARGUMENTS.get('compiler', '')
    print("Using compiler", compiler)
    if (len(compiler) > 0):
        env = DefaultEnvironment(tools = [compiler],
                                 ENV = {'PATH' : os.environ['PATH'],
                                        'TMP' : os.environ['TMP'],
                                        'TEMP' : os.environ['TEMP']})
        
# Use the current environment if nothing was specified.
if env is None:
    env = Environment(ENV = os.environ)


# C compiler flags.
#
# FIXME: Visual C flags?
if (env['CC'] == "gcc"):
    env.Append(CCFLAGS = ['-O3','-Wall'],
               LINKFLAGS = ['-Wl,-z,defs'])


# Library names and paths.
fftw_lib = 'fftw3'
fftw_lib_path = None

# OS-X apparently does not have and doesn't need the rt library.
if (platform.system() == "Darwin"):
    fftw_lib='libfftw3'
    fftw_lib_path = '/usr/local/lib'
    l1h_libs = ['lapack', 'm']
    env.Append(CCFLAGS='-I/usr/local/include')
    env.Append(LDFLAGS='-L/usr/local/include')
else:
    l1h_libs = ['lapack', 'rt', 'm']
    
lapack_lib_path = None

#
# Windows specific settings library setting. Basically we are trying
# to figure out if FFTW and LAPACK exist in the build environment or
# if we should use the versions included in this package.
#
if (platform.system() == 'Windows'):
    fftw_lib = 'fftw3-3'
    l1h_libs = ['lapack', 'm']
    conf = Configure(env)
    if not conf.CheckLib(fftw_lib):
        print("FFTW3 library not found, using storm-analysis version.")
        fftw_lib_path = ['#/storm_analysis/c_libraries/']
    if not conf.CheckLib('lapack'):
        print("LAPACK library not found, using storm-analysis version.")
        lapack_lib_path = ['#/storm_analysis/c_libraries/']


# storm_analysis/dbscan
Default(env.SharedLibrary('./storm_analysis/c_libraries/dbscan',
	                  ['./storm_analysis/dbscan/kdtree.c',
                           './storm_analysis/dbscan/dbscan.c']))


# storm_analysis/fista
if fftw_lib_path is not None:
    Default(env.SharedLibrary('./storm_analysis/c_libraries/fista_decon_utilities',
                              ['./storm_analysis/fista/fista_decon_utilities.c'],
                              CPPPATH = fftw_lib_path))
    Default(env.SharedLibrary('./storm_analysis/c_libraries/fista_fft',
                              ['./storm_analysis/fista/fista_fft.c'],
                              LIBS = [fftw_lib, 'm'], LIBPATH = fftw_lib_path, CPPPATH = fftw_lib_path))
else:
    Default(env.SharedLibrary('./storm_analysis/c_libraries/fista_decon_utilities',
                              ['./storm_analysis/fista/fista_decon_utilities.c']))
    Default(env.SharedLibrary('./storm_analysis/c_libraries/fista_fft',
                              ['./storm_analysis/fista/fista_fft.c'],
                              LIBS = [fftw_lib, 'm']))


# storm_analysis/L1H
Default(env.SharedObject(source = './storm_analysis/L1H/homotopy_common.c',
                         target = './storm_analysis/c_libraries/homotopy_common.o'))

Default(env.SharedObject(source = './storm_analysis/L1H/homotopy_imagea.c',
                         target = './storm_analysis/c_libraries/homotopy_imagea.o'))

Default(env.SharedObject(source = './storm_analysis/L1H/homotopy_imagea_common.c',
                         target = './storm_analysis/c_libraries/homotopy_imagea_common.o'))

Default(env.SharedObject(source = './storm_analysis/L1H/homotopy_sse.c',
                         target = './storm_analysis/c_libraries/homotopy_sse.o'))

Default(env.SharedObject(source = './storm_analysis/L1H/homotopy_storm.c',
                         target = './storm_analysis/c_libraries/homotopy_storm.o'))

if lapack_lib_path is not None:
    Default(env.SharedLibrary('./storm_analysis/c_libraries/homotopy_general',
                              ['./storm_analysis/L1H/homotopy_general.c',
                               './storm_analysis/c_libraries/homotopy_common.o'],
                              LIBS = l1h_libs, LIBPATH = lapack_lib_path))
else:
    Default(env.SharedLibrary('./storm_analysis/c_libraries/homotopy_general',
                              ['./storm_analysis/L1H/homotopy_general.c',
                               './storm_analysis/c_libraries/homotopy_common.o'],
                              LIBS = l1h_libs))

if lapack_lib_path is not None:
    Default(env.SharedLibrary('./storm_analysis/c_libraries/homotopy_ia_sse',
                              ['./storm_analysis/c_libraries/homotopy_imagea.o',
                               './storm_analysis/c_libraries/homotopy_sse.o',
                               './storm_analysis/c_libraries/homotopy_imagea_common.o',
                               './storm_analysis/c_libraries/homotopy_common.o'],
                              LIBS = l1h_libs, LIBPATH = lapack_lib_path))
else:
    Default(env.SharedLibrary('./storm_analysis/c_libraries/homotopy_ia_sse',
                              ['./storm_analysis/c_libraries/homotopy_imagea.o',
                               './storm_analysis/c_libraries/homotopy_sse.o',
                               './storm_analysis/c_libraries/homotopy_imagea_common.o',
                               './storm_analysis/c_libraries/homotopy_common.o'],
                              LIBS = l1h_libs))

if lapack_lib_path is not None:
    Default(env.SharedLibrary('./storm_analysis/c_libraries/homotopy_ia_storm',
                              ['./storm_analysis/c_libraries/homotopy_imagea.o',
                               './storm_analysis/c_libraries/homotopy_storm.o',
                               './storm_analysis/c_libraries/homotopy_imagea_common.o',
                               './storm_analysis/c_libraries/homotopy_common.o'],
                              LIBS = l1h_libs, LIBPATH = lapack_lib_path))
else:
    Default(env.SharedLibrary('./storm_analysis/c_libraries/homotopy_ia_storm',
                              ['./storm_analysis/c_libraries/homotopy_imagea.o',
                               './storm_analysis/c_libraries/homotopy_storm.o',
                               './storm_analysis/c_libraries/homotopy_imagea_common.o',
                               './storm_analysis/c_libraries/homotopy_common.o'],
                              LIBS = l1h_libs))

if lapack_lib_path is not None:
    Default(env.SharedLibrary('./storm_analysis/c_libraries/homotopy_sse',
                              ['./storm_analysis/c_libraries/homotopy_sse.o',
                               './storm_analysis/c_libraries/homotopy_common.o'],
                              LIBS = l1h_libs, LIBPATH = lapack_lib_path))
else:
    Default(env.SharedLibrary('./storm_analysis/c_libraries/homotopy_sse',
                              ['./storm_analysis/c_libraries/homotopy_sse.o',
                               './storm_analysis/c_libraries/homotopy_common.o'],
                              LIBS = l1h_libs))

if lapack_lib_path is not None:
    Default(env.SharedLibrary('./storm_analysis/c_libraries/homotopy_storm',
                              ['./storm_analysis/c_libraries/homotopy_storm.o',
                               './storm_analysis/c_libraries/homotopy_common.o'],
                              LIBS = l1h_libs, LIBPATH = lapack_lib_path))
else:
    Default(env.SharedLibrary('./storm_analysis/c_libraries/homotopy_storm',
                              ['./storm_analysis/c_libraries/homotopy_storm.o',
                               './storm_analysis/c_libraries/homotopy_common.o'],
                              LIBS = l1h_libs))


# storm_analysis/rolling_ball_bgr
Default(env.SharedLibrary('./storm_analysis/c_libraries/rolling_ball_lib',
	                  ['./storm_analysis/rolling_ball_bgr/rolling_ball_lib.c']))


# storm_analysis/frc
Default(env.SharedLibrary('./storm_analysis/c_libraries/frc',
	                  ['./storm_analysis/frc/frc.c'],
                          LIBS = ['m']))


# storm_analysis/multi_plane
#if lapack_lib_path is not None:
#    Default(env.SharedLibrary('./storm_analysis/c_libraries/mp_fit',
#                              ['./storm_analysis/multi_plane/mp_fit.c',
#                               './storm_analysis/c_libraries/cubic_fit.o',
#                               './storm_analysis/c_libraries/cubic_spline.o',
#                               './storm_analysis/c_libraries/multi_fit.o'],
#                              LIBS = ['lapack'], LIBPATH = lapack_lib_path))
#else:
#    Default(env.SharedLibrary('./storm_analysis/c_libraries/mp_fit',
#                              ['./storm_analysis/multi_plane/mp_fit.c',
#                               './storm_analysis/c_libraries/cubic_fit.o',
#                               './storm_analysis/c_libraries/cubic_spline.o',
#                               './storm_analysis/c_libraries/multi_fit.o'],
#                              LIBS = ['lapack']))
#
#Default(env.SharedLibrary('./storm_analysis/c_libraries/mp_utilities',
#	                  ['./storm_analysis/multi_plane/mp_utilities.c']))


# storm_analysis/sa_library
Default(env.SharedObject(source = './storm_analysis/sa_library/multi_fit.c',
                         target = './storm_analysis/c_libraries/multi_fit.o'))

if lapack_lib_path is not None:
    Default(env.SharedLibrary('./storm_analysis/c_libraries/dao_fit',
                              ['./storm_analysis/sa_library/dao_fit.c',
                               './storm_analysis/sa_library/multi_fit.c'],
                              LIBS = ['lapack', 'm'], LIBPATH = lapack_lib_path))
else:
    Default(env.SharedLibrary('./storm_analysis/c_libraries/dao_fit',
                              ['./storm_analysis/sa_library/dao_fit.c',
                               './storm_analysis/sa_library/multi_fit.c'],
                              LIBS = ['lapack', 'm']))

Default(env.SharedLibrary('./storm_analysis/c_libraries/affine_transform',
	                 ['./storm_analysis/sa_library/affine_transform.c']))

Default(env.SharedLibrary('./storm_analysis/c_libraries/grid',
	                 ['./storm_analysis/sa_library/grid.c']))

Default(env.SharedLibrary('./storm_analysis/c_libraries/ia_utilities',
	                  ['./storm_analysis/sa_library/ia_utilities.c'],
                          LIBS = ['m']))

if fftw_lib_path is not None:
    Default(env.SharedLibrary('./storm_analysis/c_libraries/matched_filter',
                              ['./storm_analysis/sa_library/matched_filter.c'],
                              LIBS = [fftw_lib], LIBPATH = fftw_lib_path, CPPPATH = fftw_lib_path))
else:
    Default(env.SharedLibrary('./storm_analysis/c_libraries/matched_filter',
                              ['./storm_analysis/sa_library/matched_filter.c'],
                              LIBS = [fftw_lib]))


# storm_analysis/sa_utilities
Default(env.SharedLibrary('./storm_analysis/c_libraries/apply-drift-correction',
	                  ['./storm_analysis/sa_utilities/apply-drift-correction.c']))

Default(env.SharedLibrary('./storm_analysis/c_libraries/avemlist',
	                  ['./storm_analysis/sa_utilities/avemlist.c'],
                          LIBS = ['m']))

Default(env.SharedLibrary('./storm_analysis/c_libraries/fitz',
	                  ['./storm_analysis/sa_utilities/fitz.c'],
                          LIBS = ['m']))

Default(env.SharedLibrary('./storm_analysis/c_libraries/tracker',
	                  ['./storm_analysis/sa_utilities/tracker.c'],
                          LIBS = ['m']))


# storm_analysis/simulator
Default(env.SharedLibrary('./storm_analysis/c_libraries/draw_gaussians',
	                  ['./storm_analysis/simulator/draw_gaussians.c'],
                          LIBS = ['m']))

Default(env.SharedLibrary('./storm_analysis/c_libraries/zernike',
	                  ['./storm_analysis/simulator/zernike.c'],
                          LIBS = ['m']))


# storm_analysis/spliner
Default(env.SharedObject(source = './storm_analysis/spliner/cubic_fit.c',
                         target = './storm_analysis/c_libraries/cubic_fit.o'))

Default(env.SharedObject(source = './storm_analysis/spliner/cubic_spline.c',
                         target = './storm_analysis/c_libraries/cubic_spline.o'))

Default(env.SharedLibrary('./storm_analysis/c_libraries/cubic_spline',
	                  ['./storm_analysis/c_libraries/cubic_spline.o']))
                           
if lapack_lib_path is not None:
    Default(env.SharedLibrary('./storm_analysis/c_libraries/cubic_fit',
                              ['./storm_analysis/spliner/cubic_fit.c',
                               './storm_analysis/c_libraries/cubic_spline.o',
                               './storm_analysis/c_libraries/multi_fit.o'],
                              LIBS = ['lapack', 'm'], LIBPATH = lapack_lib_path))
else:
    Default(env.SharedLibrary('./storm_analysis/c_libraries/cubic_fit',
                              ['./storm_analysis/spliner/cubic_fit.c',
                               './storm_analysis/c_libraries/cubic_spline.o',
                               './storm_analysis/c_libraries/multi_fit.o'],
                              LIBS = ['lapack', 'm']))

