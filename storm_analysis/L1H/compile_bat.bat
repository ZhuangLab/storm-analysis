c:\MinGW64\bin\x86_64-w64-mingw32-gcc.exe -c homotopy_common.c -O3
c:\MinGW64\bin\x86_64-w64-mingw32-gcc.exe -c homotopy_general.c -O3
c:\MinGW64\bin\x86_64-w64-mingw32-gcc.exe -c homotopy_storm.c -O3
c:\MinGW64\bin\x86_64-w64-mingw32-gcc.exe -c homotopy_sse.c -O3
c:\MinGW64\bin\x86_64-w64-mingw32-gcc.exe -c homotopy_imagea_common.c -O3
c:\MinGW64\bin\x86_64-w64-mingw32-gcc.exe -c homotopy_imagea.c -O3

c:\MinGW64\bin\x86_64-w64-mingw32-gcc.exe -shared -o homotopy_general.dll homotopy_general.o homotopy_common.o -llapack -Lc:\Users\Hazen\lib
c:\MinGW64\bin\x86_64-w64-mingw32-gcc.exe -shared -o homotopy_storm.dll homotopy_storm.o homotopy_common.o -llapack -Lc:\Users\Hazen\lib
c:\MinGW64\bin\x86_64-w64-mingw32-gcc.exe -shared -o homotopy_sse.dll homotopy_sse.o homotopy_common.o -llapack -Lc:\Users\Hazen\lib

c:\MinGW64\bin\x86_64-w64-mingw32-gcc.exe -shared -o homotopy_ia_storm.dll homotopy_imagea.o homotopy_storm.o homotopy_imagea_common.o homotopy_common.o -llapack -Lc:\Users\Hazen\lib
c:\MinGW64\bin\x86_64-w64-mingw32-gcc.exe -shared -o homotopy_ia_sse.dll homotopy_imagea.o homotopy_sse.o homotopy_imagea_common.o homotopy_common.o -llapack -Lc:\Users\Hazen\lib

