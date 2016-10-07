#!/bin/bash

gcc -fPIC -g -c -Wall -O3 multi_fit.c
gcc -fPIC -g -c -Wall -O3 dao_fit.c
gcc -shared -Wl,-soname,dao_fit.so.1 -o dao_fit.so.1.0.1 dao_fit.o multi_fit.o -llapack
ln -s dao_fit.so.1.0.1 dao_fit.so

gcc -fPIC -g -c -Wall -O3 ia_utilities.c
gcc -shared -Wl,-soname,ia_utilities.so.1 -o ia_utilities.so.1.0.1 ia_utilities.o
ln -s ia_utilities.so.1.0.1 ia_utilities.so

gcc -fPIC -g -c -Wall -O3 grid.c
gcc -shared -Wl,-soname,grid.so.1 -o grid.so.1.0.1 grid.o
ln -s grid.so.1.0.1 grid.so

gcc -fPIC -g -c -Wall matched_filter.c
gcc -shared -Wl,-soname,matched_filter.so.1 -o matched_filter.so.1.0.1 matched_filter.o -lc -lfftw3
ln -s matched_filter.so.1.0.1 matched_filter.so
