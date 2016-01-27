gcc -fPIC -g -c -Wall fista_decon_utilities.c
gcc -shared -Wl,-soname,fista_decon_utilities.so.1 -o fista_decon_utilities.so.1.0.1 fista_decon_utilities.o -lc
ln -s fista_decon_utilities.so.1.0.1 fista_decon_utilities.so
