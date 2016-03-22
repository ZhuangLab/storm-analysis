gcc -fPIC -g -c -Wall fista_decon_utilities.c
gcc -shared -Wl,-soname,fista_decon_utilities.so.1 -o fista_decon_utilities.so.1.0.1 fista_decon_utilities.o -lc
ln -s fista_decon_utilities.so.1.0.1 fista_decon_utilities.so

gcc -fPIC -g -c -Wall -O3 fista_fft.c
gcc -shared -Wl,-soname,fista_fft.so.1 -o fista_fft.so.1.0.1 fista_fft.o -lc -lfftw3
ln -s fista_fft.so.1.0.1 fista_fft.so
