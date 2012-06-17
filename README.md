attsim
===========

attsim is a super basic rigid-body attitude dynamics simulator, using
C and GSL.  

Installation
-------------

Install the GNU Scientific Library and the BLAS and ATLAS packages for
linear algebra

    sudo apt-get install libgsl0-dev libatlas3gf-base libblas3gf libatlas-dev libblas-dev libatlas-base-dev

Install GNUPlot to generate nice PDF plots

    gnuplot-x11

Clone this repository

    git clone git://github.com/peddie/attsim.git
    cd attsim

Clone and build Greg Horn's mathlib library for spatial rotations

    git clone git://github.com/ghorn/mathlib.git
    pushd mathlib && make && popd

Build the simulator

    make

Run
-----------

    ./sim <seconds>

Analyze
-----------

    gnuplot results.plot
    <pdf viewer> *.pdf

Extend
-----------

    emacs dynamics.c
    <hack>
    make && ./sim 222 && gnuplot results.plot

Features
-----------

- Simulates attitude dynamics of a body subjected to arbitrary
  body-frame or inertial-frame torques.

- Currently allows for feedback control via measurement gyros and
  magnetometer and using magnetic torquers

- Fast (approximately 7000-10000 times realtime)

Bugs
-----------

- Not sure yet whether it does torque-free precession correctly.

- Not sure yet whether the bdot control is reasonable (just made some
  stuff up).

If you find anything please let us know!

Coming Soon
-----------

- OpenGL output!

- Unit tests!

