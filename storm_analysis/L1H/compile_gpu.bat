c:\MinGW64\bin\x86_64-w64-mingw32-gcc.exe -c -msse2 homotopy_gpu.c -O3 -I"c:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v5.0\include
c:\MinGW64\bin\x86_64-w64-mingw32-gcc.exe -shared -o homotopy_gpu.dll homotopy_gpu.o -lopencl -L"c:\Users\hazen\hbabcock\C\opencl" -llapack -Lc:\Users\Hazen\lib
