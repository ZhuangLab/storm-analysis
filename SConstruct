#!python

import platform

# Configure build environment.
env = Environment()

# C compiler flags.
env.Append(CCFLAGS = ['-O3'])

# Windows specific
if (platform.system() == 'Windows'):
    env.Append(CPPPATH = ['./storm_analysis/c_libraries/'])
    env.Append(LIBPATH = ['./storm_analysis/c_libraries/'])


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
                          LIBS = ['fftw3']))


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

