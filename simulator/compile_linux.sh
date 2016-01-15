gcc -fPIC -g -c -Wall draw_gaussians.c
gcc -shared -Wl,-soname,draw_gaussians.so.1 -o draw_gaussians.so.1.0.1 draw_gaussians.o -lc
ln -s draw_gaussians.so.1.0.1 draw_gaussians.so
