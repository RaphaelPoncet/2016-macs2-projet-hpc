# Explict PDE engine for wave propagation

## Summary
This code solves the acoustic wave equation for heterogeneous
media (e.g. variable characteristic wave speed):

It solves the above equations on cartesian grid, using a classical
second order centered finite differences in time and space.

It is open source (Apache 2 license).

## Features

* flexible input: in the parameter file, you can input initial conditions and
  characteristic velocity using either binary files or mathematic formulas,
  without modifying the code.
* flexible variables output: you can output the grid variables either
  in VTK or simple gridded binary and ascii format. The VTK files can be
  opened in [paraview](http://www.paraview.org) or any visualisation
  software supporting the VTK file formats. The gridded binary and ascii
  format can -- for example -- be read and processed using python. You will
  find two python scripts (`plot_output.py` and `read_receivers.py`) in the directory.
* periodical boundary conditions. Optionally, the user can specify
  imperfect absorbing boundary conditions using the sponge layer technique (Cerjean et al 1985).

## Install, build and execution

### Dependencies

The code is written in C/C++, and has a few optional python scripts
for post-processing results. It depends on [plog](https://github.com/SergiusTheBest/plog) for logging,
[muparser](http://beltoforion.de/article.php?a=muparser)
for parsing mathematical formulas, and [picojson](https://github.com/kazuho/picojson) (for parsing the
JSON parameter file). However the source of these libraries are bundled
in the git repository, so you don't need to install anything.

### Generating the executable

The main steps are:

* get the source code on your machine

        git clone https://github.com/RaphaelPoncet/2016-macs2-projet-hpc

* go to the source directory

        cd 2016-macs2-projet-hpc

* compile

        make wave.exe

* (optional) if the compilation fails with a message involving the muparser
  library, compile the library by hand:

        cd ./external/muparser-2.2.5/
        ./configure --enable-shared = no
        make
        cd ../../

### Execution

The executable is named `wave.exe`. To execute the code, you need to provide a *parameter file*:

    ./wave.exe param_file.json

This file is in [json](http://www.json.org) format. The repository
comes with 4 examples of configuration files:

* convergence.json : propagation on a 1D plane wave in a homognenous
  medium. This file demonstrates the simplest way to study the convergence of the numerical scheme.

* homogeneous.json : propagation of a gaussian bump in an homogeneous medium.

    ![Pressure at final time](./images/radial_homogeneous.png)

* marmousi.json : propagation of a gaussian bump in the original
  Marmousi model that mimics a geologically realistic velocity model
  of the earth subsurface

* marmousi2.json: propagation of a gaussian bump in the Marmousi 2
  velocity model (an extension of the original Marmousi model)

## How to use the parameter file

