gcc -fPIC -g -c -Wall -O3 rolling_ball_lib.c
gcc -shared -Wl,-soname,rolling_ball_lib.so.1 -o rolling_ball_lib.so.1.0.1 rolling_ball_lib.o -lc
ln -s rolling_ball_lib.so.1.0.1 rolling_ball_lib.so
