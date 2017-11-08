
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
