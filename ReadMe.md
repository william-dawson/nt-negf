# Interoperability Prototype

This project attempts to use NTPoly as glue code for integrating libNEGF.

## Building

We have some dependencies which we can install in the external directory.
```
git submodule update --init --recursive
```
Now compile NTPoly.
```
cd external/NTPoly/Build
cmake .. -DCMAKE_INSTALL_PREFIX=../../
make -j
make install
cd ../../
```
And libNEGF
```
mkdir libnegf/build
cd libnegf/build 
cmake .. -DCMAKE_INSTALL_PREFIX=../../
make -j
make install
cd ../../
```
Now we can compile the prototype.
```
make
```

## Data Generation

We need some metadata file for understanding the indexing. 
```
pip install PyBigDFT
# install BigDFT and put it on your path
python generate.py example/sw.xyz example/log-sw.yaml 12 24 36
```
In this case, the first number 12 is the last atom in the contact region,
24 the last atom in the left region, and 36 the last atom in the right region.
The logfile log-sw.yaml comes from a BigDFT calculation.
