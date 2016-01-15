C:\MinGW64\bin\x86_64-w64-mingw32-gcc -c -O3 cubic_spline.c
C:\MinGW64\bin\x86_64-w64-mingw32-gcc -c -O3 multi_fit_core.c
C:\MinGW64\bin\x86_64-w64-mingw32-gcc -c -O3 cubic_fit.c
C:\MinGW64\bin\x86_64-w64-mingw32-gcc -shared -o cubic_spline.dll cubic_spline.o
C:\MinGW64\bin\x86_64-w64-mingw32-gcc -shared -o cubic_fit.dll cubic_fit.o multi_fit_core.o cubic_spline.o -llapack -Lc:\Users\Hazen\lib

