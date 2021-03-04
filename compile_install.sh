set -x
cd ./Release/ 
#rm CMakeCache.txt
cmake .. -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON 
cmake --build ../Release  --target help 
cmake --build ../Release  --target all -j 30 
cmake --install . --prefix ../build/
