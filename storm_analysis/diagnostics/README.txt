
Scripts for testing and evaluating localization finding and fitting performance.

The basic layout is that for each type of analysis we have:

1) settings.py - Settings for the simulation.

2) configure.py - Create the files, XML, etc. necessary for testing.

3) make_data.py - Create simulate SMLM movies.

4) analyze_data.py - Analyze the simulated dataset(s).

5) collate.py - Measure how well the analysis performed and create a summary of
                the results.

All of these should be run in a working directory.

Typically you would run (2) once to set everything up. Then you'd repeat (3) - (5),
changing (1) and (3) as needed, to test analysis performance with different types of
simulated data.


Running without a display:

Some of the analysis these call will open a plot window and wait for you to close
it, which is not what you want when running several diagnostics in a row. Setting

  STORM_ANALYSIS_HEADLESS=1

in the environment selects a non-interactive matplotlib backend, so pyplot.show()
does nothing anywhere in the package. Figures that are written to disk, such as the
mapping images from the multicolor diagnostic, are still written.

Note that multiplane/configure.py needs --psf-model, one of psf_fft, pupilfn or
spline. It is required and has no default, so configure.py exits without it.

That is the only required argument in any of these. fista_decon, spliner,
spliner_2d and multiplane also accept an optional --no-splines. The rest run
with no arguments.


(Linux) C profiling tools:
1. http://valgrind.org/docs/manual/cl-manual.html
2. http://kcachegrind.sourceforge.net/html/Home.html

Briefly:
$ valgrind --tool=callgrind python xyzzy.py
$ KCachegrind

Note that running in valgrind will take 5-10x longer than normal.


Python profiling tools:
https://docs.python.org/3.5/library/profile.html

$ python -m cProfile -o prof.prof xyzzy.py
$ snakeviz prof.prof
