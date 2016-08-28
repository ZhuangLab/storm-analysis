gcc -fPIC -g -c -Wall -O3 draw_gaussians.c
gcc -shared -Wl,-soname,draw_gaussians.so.1 -o draw_gaussians.so.1.0.1 draw_gaussians.o -lc
ln -s draw_gaussians.so.1.0.1 draw_gaussians.so

gcc -fPIC -g -c -Wall -O3 zernike.c
gcc -shared -Wl,-soname,zernike.so.1 -o zernike.so.1.0.1 zernike.o -lc
ln -s zernike.so.1.0.1 zernike.so
