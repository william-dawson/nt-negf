# Interoperability Prototype

This project attempts to use NTPoly as glue code for integrating libNEGF.

## Building

First we need to download the NTPoly source code and compile it.
```
git clone git@github.com:william-dawson/NTPoly.git 
cd NTPoly/Build
cmake .. -DCMAKE_INSTALL_PREFIX=../../
make -j
```
And libNEGF
```
git clone git@github.com:libnegf/libnegf.git
cd libnegf 
cmake .. -DCMAKE_INSTALL_PREFIX=../../
make -j
```
Now we can compile the prototype.
```
make
make test
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
