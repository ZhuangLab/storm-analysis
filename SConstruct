#!python

import os
import platform

# Configure build environment.
env = Environment()
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
        env = Environment(tools = [compiler])
        env.Append(PATH = os.environ['PATH'])
#                                                     'TMP' : os.environ['TMP'],
#                                                     'TEMP' : os.environ['TEMP']})

# C compiler flags.
env.Append(CCFLAGS = ['-O3'])

    
# Windows specific
if (platform.system() == 'Windows'):
    fftw_lib = 'fftw3-3'
    env.Append(CPPPATH = ['./storm_analysis/c_libraries/'])
    env.Append(LIBPATH = ['./storm_analysis/c_libraries/'])
else:
    fftw_lib = 'fftw3'


# storm_analysis/dbscan
Default(env.SharedLibrary('./storm_analysis/c_libraries/dbscan',
	                  ['./storm_analysis/dbscan/kdtree.c',
                           './storm_analysis/dbscan/dbscan.c']))


# storm_analysis/fista
Default(env.SharedLibrary('./storm_analysis/c_libraries/fista_decon_utilities',
	                  ['./storm_analysis/fista/fista_decon_utilities.c']))

Default(env.SharedLibrary('./storm_analysis/c_libraries/fista_fft',
	                  ['./storm_analysis/fista/fista_fft.c'],
                          LIBS = [fftw_lib]))


# storm_analysis/L1H
Default(env.SharedLibrary('./storm_analysis/c_libraries/homotopy_general',
	                  ['./storm_analysis/L1H/homotopy_general.c',
                           './storm_analysis/L1H/homotopy_common.c'],
                          LIBS = ['lapack']))

Default(env.SharedLibrary('./storm_analysis/c_libraries/homotopy_ia_sse',
	                  ['./storm_analysis/L1H/homotopy_imagea.c',
                           './storm_analysis/L1H/homotopy_sse.c',
                           './storm_analysis/L1H/homotopy_imagea_common.c',
                           './storm_analysis/L1H/homotopy_common.c'],
                          LIBS = ['lapack']))

Default(env.SharedLibrary('./storm_analysis/c_libraries/homotopy_ia_storm',
	                  ['./storm_analysis/L1H/homotopy_imagea.c',
                           './storm_analysis/L1H/homotopy_storm.c',
                           './storm_analysis/L1H/homotopy_imagea_common.c',
                           './storm_analysis/L1H/homotopy_common.c'],
                          LIBS = ['lapack']))

Default(env.SharedLibrary('./storm_analysis/c_libraries/homotopy_sse',
	                  ['./storm_analysis/L1H/homotopy_sse.c',
                           './storm_analysis/L1H/homotopy_common.c'],
                          LIBS = ['lapack']))

Default(env.SharedLibrary('./storm_analysis/c_libraries/homotopy_storm',
	                  ['./storm_analysis/L1H/homotopy_storm.c',
                           './storm_analysis/L1H/homotopy_common.c'],
                          LIBS = ['lapack']))


# storm_analysis/rolling_ball_bgr
Default(env.SharedLibrary('./storm_analysis/c_libraries/rolling_ball_lib',
	                  ['./storm_analysis/rolling_ball_bgr/rolling_ball_lib.c']))


# storm_analysis/frc
Default(env.SharedLibrary('./storm_analysis/c_libraries/frc',
	                  ['./storm_analysis/frc/frc.c']))


# storm_analysis/sa_library
Default(env.SharedLibrary('./storm_analysis/c_libraries/dao_fit',
	                  ['./storm_analysis/sa_library/dao_fit.c',
                           './storm_analysis/sa_library/multi_fit.c'],
                          LIBS = ['lapack']))

Default(env.SharedLibrary('./storm_analysis/c_libraries/grid',
	                 ['./storm_analysis/sa_library/grid.c']))

Default(env.SharedLibrary('./storm_analysis/c_libraries/ia_utilities',
	                  ['./storm_analysis/sa_library/ia_utilities.c']))

Default(env.SharedLibrary('./storm_analysis/c_libraries/matched_filter',
	                  ['./storm_analysis/sa_library/matched_filter.c'],
                          LIBS = [fftw_lib]))


# storm_analysis/sa_utilities
Default(env.SharedLibrary('./storm_analysis/c_libraries/apply-drift-correction',
	                  ['./storm_analysis/sa_utilities/apply-drift-correction.c']))

Default(env.SharedLibrary('./storm_analysis/c_libraries/avemlist',
	                  ['./storm_analysis/sa_utilities/avemlist.c']))

Default(env.SharedLibrary('./storm_analysis/c_libraries/fitz',
	                  ['./storm_analysis/sa_utilities/fitz.c']))

Default(env.SharedLibrary('./storm_analysis/c_libraries/tracker',
	                  ['./storm_analysis/sa_utilities/tracker.c']))


# storm_analysis/simulator
Default(env.SharedLibrary('./storm_analysis/c_libraries/draw_gaussians',
	                 ['./storm_analysis/simulator/draw_gaussians.c']))

Default(env.SharedLibrary('./storm_analysis/c_libraries/zernike',
	                 ['./storm_analysis/simulator/zernike.c']))



# storm_analysis/sCMOS
Default(env.SharedLibrary('./storm_analysis/c_libraries/scmos_utilities',
	                  ['./storm_analysis/sCMOS/scmos_utilities.c']))


# storm_analysis/spliner
Default(env.SharedLibrary('./storm_analysis/c_libraries/cubic_spline',
	                  ['./storm_analysis/spliner/cubic_spline.c']))
                           
Default(env.SharedLibrary('./storm_analysis/c_libraries/cubic_fit',
	                  ['./storm_analysis/spliner/cubic_fit.c',
                           './storm_analysis/spliner/cubic_spline.c',
                           './storm_analysis/sa_library/multi_fit.c'],
                          LIBS = ['lapack']))
