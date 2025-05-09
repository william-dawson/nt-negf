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

And example python script is included 