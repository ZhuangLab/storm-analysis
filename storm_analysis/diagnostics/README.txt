
Scripts for testing and evaluating localization finding and fitting performance.

The basic layout is that for each type of analysis we have:

1) configure.py - Create the files, XML, etc. necessary for testing.

2) make_data.py - Create simulate SMLM movies.

3) analyze_data.py - Analyze the simulated dataset(s).

4) collate.py - Measure how well the analysis performed and create a summary of
                the results.

All of these should be run in a working directory.

Typically you would run (1) once to set everything up. Then you'd repeat (2) - (4),
changing (2) as needed, to test analysis performance with different types of
simulated data.
