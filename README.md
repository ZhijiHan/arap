# As-Rigid-As-Possible implementation

## Overview

This small project implements the algorithm described in the paper [As-Rigid-As-Possible surface modeling](http://sites.fas.harvard.edu/~cs277/papers/sorkine_asrigid.pdf). It uses [libigl](http://igl.ethz.ch/projects/libigl/), [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) and [CHOLMOD](http://www.cise.ufl.edu/research/sparse/cholmod/)(http://eigen.tuxfamily.org/index.php?title=Main_Page) for UI, matrix manipulation and optimization.

## How to run the demo

Demos can be found in [this youtube link](https://www.youtube.com/watch?v=uK000vVJCes). Moreover, libigl's implementation can be found [here](https://www.youtube.com/watch?v=-BQWvz1zyOI).

Unfortunately you have to configure libigl and Eigen first. Put them in any folder you want, and change any relevant path in all the cmake files in libigl and this project. Moreover, make sure you can use all the GL libs correctly (glew, glfw, etc).

The mesh files used in this demo is from libigl's tutorials. Please copy them into some places you like, and change file paths in build.sh accordingly.

Now you can run ./build.sh in the root directory. If everything goes well, you can see a decimated knight moving on the screen after pressing space.

## About the code

The project contains 3 files: main.cc, arapsolver.h and arapsolver.cc.

### main.cc

libigl contains its own implementation for ARAP algorithm. By default their implementation is used as the benchmark to test whether our code is correct. You can set USE_IGL_AS_BENCHMARK to 0 if you don't want to use it.

If USE_IGL_AS_BENCHMARK is enabled, ***pre_draw*** function will compare libigl's solution and our own implementation. If the error goes beyond our pre defined threshold, the program will complain.

### arapsolver.h

This file declares an ArapSolver class. The two most important methods are ***Precompute*** and ***Solve***. In ***Precompute***, the cotangent weights for all the triangles are computed, and the Laplace-Beltrami operator are computed and decomposed. ***Solve*** takes a couple of fixed vertices, and optimizes alternatively between rotations and vertices.

### arapsolver.cc

This file defines all the methods in ArapSolver. Please refer to the comments for details.