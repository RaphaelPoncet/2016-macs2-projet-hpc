# Explicit PDE engine for wave propagation

## Summary
This code solves the acoustic wave equation for heterogeneous
media (e.g. variable characteristic wave speed):

![Acoustic wave equation](./images/wave-eq.png)

on a cartesian grid, using a classical
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

* `convergence.json` : propagation on a 1D plane wave in a homognenous
  medium. This file demonstrates the simplest way to study the convergence of the numerical scheme.

* `homogeneous.json` : propagation of a gaussian bump in an homogeneous medium.

    ![Pressure at final time](./images/radial_homogeneous.png)

* `marmousi.json` : propagation of a gaussian bump in the original
  Marmousi model that mimics a geologically realistic velocity model
  of the earth subsurface

* `marmousi2.json` : propagation of a gaussian bump in the Marmousi 2
  velocity model (an extension of the original Marmousi model)

    ![Marmousi 2 pressure snapshot](./images/marmousi2.png)


    The above image represent a snapshot of pressure waves propagating
    through the model. Below, we repesent a *shot gather*, e.g. the
    recording of pressure time series on an array of receivers.

    ![Marmousi 2 shot gather](./images/marmousi2_shot.png)

## How to use the parameter file

### Looking at an example

Let us look in detail at the `convergence.json` parameter file.


    "parameters" : ["x0=0.5*(xmin+xmax)", 
                    "z0=0.5*(zmin+zmax)", 
                    "lambda=100.0", 
                    "V0=0.1"],

The *parameters* field define a few parameters that can be used
subsequently in the remainder of the file.

    "init" : {"velocity" : {"formula" : "V0"},
              "pressure_0" : {"formula" : "exp(-lambda*((z-z0)^2))"},
              "pressure_1" : {"formula" : "pressure_0"}},

The *init* field defines initialization of grid variables. The grid
variables are `pressure_0` and `pressure_1` (the acoustic pressure at
2 successive time steps) , `velocity` (the velocity in the grid, which
is a coefficient of the equation), `laplace_p` (storing the pressure
derivatives) and `pressure_ref` (helpful for storing analytical
formulas). In this case, all variables are defined using mathematical
formulas, but one can also initialize them from a file (see
`marmousi.json` or `marmousi2.json`).

* for the velocity, we define it as a constant, `V0`, which has been
  defined in the `parameters` section.

* the `pressure_0` variable is defined as a gaussian, centered on `z0`.

* the `pressure_1` variable is defined as equal to `pressure_0`

        "grid" : {"nx" : 100,
                  "ny" : 1,
                  "nz" : 200,
                  "xmin" : 0.0,
                  "xmax" : 10.0,
                  "ymin" : 0.0,
                  "ymax" : 0.0,
                  "zmin" : 0.0,
                  "zmax" : 10.0},

The `grid` section defines the geometry of the computational grid. We
simply define the dimensions of the grid `nx`, `ny` and `nz`, and its
extents in the 3 dimensions `xmin`, `xmax` and so on. *Important*: for
a 2D grid, `ny` *must be set to 1* (e.g. in 2D, the dimensions are x and
z, not x and y).

    "timeloop" : {"dt": ".99*CFL",
                  "tfinal": 10.0},

The `timeloop` section relates to the main loop parameters. First, we
must specify the timestep `dt`. It is the user responsibility to use a
timestep small enough so that the scheme is stable. Fortunately, there
is a special variable, `CFL`, which gives the Courant-Friedrichs-Levy
stability condition. Any number smaller than `CFL` should result in a
stable scheme. Here, we take the timestep to be 99% of the CFL number.
Then, we must specify the maximum simulation time. Alternatively,
*instead of* `tfinal`, one can specify `niter`, the maximum number of
iterations. *Exactly one* of `tfinal`, `niter` must be set.

    "output" : [ {"type" : "EvalVariable", 
                  "rhythm" : "end", 
                  "name": "pressure_ref",
                  "formula" : "0.5*(exp(-lambda*((z-z0-t*V0)^2)) + exp(-lambda*((z-z0+t*V0)^2)))"},
                 {"type" : "CheckVariables",
                  "rhythm" : 100},
                 {"type" : "OutputVariables", 
                  "rhythm" : "end", 
                  "format": "gridded",
                  "file" : "./output/out_variables_%i.dat"},
                 {"type" : "OutputVariables", 
                  "rhythm" : 10, 
                  "format": "VTK",
                  "file" : "./output/out_variables_%i.vtr"},
                 {"type" : "OutputReceivers", 
                  "rhythm" : 10, 
                  "format": "gridded",
                  "iz" : 200,
                  "file" : "./output/receivers.dat"},
                 {"type" : "OutputNorm",
                  "rhythm" : "end",
                  "name" : "pressure_0 - pressure_ref",
                  "file" : "output/norm.txt"}
                 ]

In that parameter file, the output section is very big, and consists
in 6 different outputs. Let us break it down.

     {"type" : "EvalVariable", 
      "rhythm" : "end", 
      "name": "pressure_ref",
      "formula" : "0.5*(exp(-lambda*((z-z0-t*V0)^2)) + exp(-lambda*((z-z0+t*V0)^2)))"}
