C:\MinGW64\bin\gcc -c -O3 cubic_spline.c
C:\MinGW64\bin\gcc -c -O3 multi_fit_core.c
C:\MinGW64\bin\gcc -c -O3 cubic_fit.c
C:\MinGW64\bin\gcc -shared -o cubic_spline.dll cubic_spline.o
C:\MinGW64\bin\gcc -shared -o cubic_fit.dll cubic_fit.o multi_fit_core.o cubic_spline.o -llapack -Lc:\Users\Hazen\storm-analysis\windows_dll
