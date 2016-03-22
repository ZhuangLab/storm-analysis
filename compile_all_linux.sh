#!/bin/bash

echo "compiling sa_library"
cd sa_library
sh compile_linux.sh

echo "compiling sa_utilities"
cd ../sa_utilities
sh compile_linux.sh

echo "compiling frc"
cd ../frc
sh compile_linux.sh

echo "compiling sCMOS"
cd ../sCMOS
sh compile_linux.sh

echo "compiling L1H"
cd ../L1H
sh compile_linux.sh

echo "compiling fista"
cd ../fista
sh compile_linux.sh

echo "compiling simulator"
cd ../simulator
sh compile_linux.sh

