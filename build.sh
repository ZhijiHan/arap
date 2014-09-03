rm -rf build
mkdir build
cd build
cmake -DCMAKE_BULID_TYPE=Release ../
make
cd ../
./build/demo_bin
