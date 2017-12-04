## Parallel DWT
[![build status](https://gitlab.com/gnthibault/PDWT/badges/NewArch/build.svg)](https://gitlab.com/gnthibault/PDWT/commits/NewArch)

PDWT is a parallel implementation of the Discrete Wavelet Transform (DWT).
This implementation in CUDA targets Nvidia GPUs.

PDWT primarily aims at being fast, simple and versatile for an easy integration in a bigger project.
For example, the easy interface and thresholding functions make it interesting for sparse regularization of inverse problems.


## Features

* 1D and 2D transform and their inverses, multi-levels, arbitrary sizes
* Support of batched 1D transform
* Separable and non-separable transforms
* DWT, SWT, and CWT in separable mode
* 72 available separable wavelets
* Custom wavelets can be defined
* Thresholding and norms utilities
* Random shift utility for translation-invariant denoising
* Simple interface (see examples)

All the transforms are computed with the **periodic boundary extension** (the dimensions are halved at each scale, except for SWT).

## Current limitations

* 3D is not handled at the moment.
* Only the periodic boundary extension is implemented.
* The parallel part is implemented in CUDA, so only Nvidia GPUs can be used.


## Installation

# Cloning this repository

The header only pybind11 module is included as a submodule of this repository. This requires some attention when cloning this project. There are two options:

*   The simplest option is:

    ```bash
    git clone --recursive https://github.com/gnthibault/PDWT.git
    ```

    This will download the submodule up to the version that is used in this project. To update to the latest commit of the submodule itself:

    ```bash
    git submodule update --remote
    ```

*   One could also directly download the submodule from the source:

    ```bash
    git clone https://github.com/gnthibault/PDWT.git
    cd ./binding/python/pybind11/
    git submodule init
    git submodule update
    ```

### Dependencies

If you want to enjoy some developer features of this project, please consider installing the following packages (example for ubuntu):  
sudo apt-get install doxygen graphviz cppcheck

If you want to enjoy mpi related features, please consider installing the following packages (example for ubuntu):  
sudo apt-get install libopenmpi-dev libopenmpi-dev

If you want to enjoy cuda related features, you need the [NVIDIA CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit), and of course a NVIDIA GPU.

### Compilation

How to build ?
cd PDWT
mkdir build; cd build  
cmake -DINSTALLDIR=/path-to-install-dir -DUSE_CUDA=ON -DUSE_NVTX=OFF -DCMAKE_BUILD_TYPE=Release -DWITH_PYTHON=ON -DPYBIND11_PYTHON_VERSION=3.6 -DTESTING_ENABLED=ON ..  
make -j8 install

## How to test
In the build directory, do:  
make test  

## How to generate doxygen documentation
In the build directory, do:  
make doc  
The documentation should be build in the "doc" directory

## How to perform static code analysis with cppcheck
In the build directory, do:  
make cppcheck  
The cppcheck report can be found in the directory cppcheckdir-report

## Getting started

### Running the example

To run the test, you need a raw image in 32 bits floating point precision format.
As PDWT was primarily written for data crunching, the I/O part is not addressed: the input and output of PDWT are float (or double, if compiled accordingly with `make libpdwtd.so`) arrays.

If you have python and scipy installed, you can generate an image input file with

```bash
cd test
python generate_image.py [Nr] [Nc]
```
where Nr, Nc are optional arguments which are the number of rows/columns of the generated image (default is 512).

You can then run an example with

```bash
make demo
./demo
```

and tune the wavelet, number of levels, etc. in the prompt.


### Calling PDWT

A typical usage would be the following :

```C
#include "wt.h"

// ...
// float* img = ...
int Nr = 1080; // number of rows of the image
int Nc = 1280; // number of columns of the image
int nlevels = 3;

// Compute the wavelets coefficients with the "Daubechies 7" wavelet
Wavelets W(img, Nr, Nc, "db7", nlevels);
W.forward();

// Do some thresholding on the wavelets coefficients
float norm1 = W1.norm1();
printf("Before threshold : L1 = %e\n", norm1);
W.soft_threshold(10.0);
norm1 = W1.norm1();
printf("After threshold : L1 = %e\n", norm1);

// Inverse the DWT and retrieve the image from the GPU
W.inverse();
W.get_image(img);
```


