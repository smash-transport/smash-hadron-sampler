# SMASH Hadron Sampler

This sampler is meant to be applied in hybrid models, simulating heavy-ion collisions in the high temperature or baryon-density region.
More precisely it provides an interface between the macroscopic hydrodynamic evolution of the fireball and the hadronic afterburner.
During the hydrodynamic evolution a hypersurface of constant energy density (the switching energy density) is created.
Each element on the hypersurface needs then be transformed into a list of particles with properties loosely provided by the macroscopic properties of the hypersurface elements.
This process of particlization is performed by means of the hadron sampler provided within this project.
It is designed to couple the [3+1D viscous hydrodynamic code vhlle](https://github.com/yukarpenko/vhlle) to the [hadronic transport model SMASH](https://smash-transport.github.io).
For details about the sampling algorithm, please consult [I. Karpenko et al., Phys. Rev. C 91 (2015) 6, 064901](https://inspirehep.net/literature/1343339).

> [!NOTE]
> Please cite the following two papers when using the SMASH-hadron-sampler:
> - [I. Karpenko et al., Phys. Rev. C 91 (2015) 6, 064901](https://inspirehep.net/literature/1343339)
> - [A. Schäfer et al., arXiv:2112.08724](https://arxiv.org/abs/2112.08724)


## Prerequisites
- [cmake](https://cmake.org) version &ge; 3.16 or higher
- [SMASH](https://github.com/smash-transport/smash) version 3.2, as well as prerequisites therein
- [ROOT](https://root.cern.ch) version &ge; 6.06

> [!IMPORTANT]
> Only tagged versions are guaranteed to be compatible with SMASH.

The minimum requirement regarding the [SMASH transport approach](https://github.com/smash-transport/smash) is that it can be used as a library.
Detailed install instructions how to compile SMASH are provided in SMASH'S [README.md](https://github.com/smash-transport/smash/blob/main/README.md) and [INSTALL.md](https://github.com/smash-transport/smash/blob/main/INSTALL.md) files.
If the output of the sampler is expected to be fed into SMASH for particle propagation after the sampling (Afterburner), then SMASH needs to be compiled entirely &ndash; not just as a library.


## Install instructions
Before getting started with the sampler installation, please check that all prerequisites are satisfied and needed software ready to be used.

To compile the project, first set the environment variable to the SMASH directory:

    export SMASH_DIR=[...]/smash

Copy the cmake files to the sampler directory:

    cd [...]/smash-hadron-sampler
    cp -r $SMASH_DIR/cmake ./

Execute the following commands to build the project:

    mkdir build
    cd build
    cmake -DPythia_CONFIG_EXECUTABLE=[...]/pythia83XX/bin/pythia8-config ..
    make

where `[...]/pythia83XX` is the path to the Pythia directory.
The `XX` needs to be exchanged to match the Pythia version that is used by the compiled SMASH library.


In continuation, the executable `sampler` is created.


## Execute the sampler
To run the sampler, execute the following command in the `build` directory:

    ./sampler --config <file>

where `<file>` needs to be the path of the configuration file.
In this config file the location of the freeze-out hypersurface file, the path to the output directory, and all other necessary parameters can be specified.

There are additional command line parameters with which the freeze-out hypersurface file, the output directory, and a number prefix for parallel runs can be specified.
All possible command line parameters are:

    -h, --help                  Print overview of command line options.
    -c, --config <file>         Mandatory parameter to specify the config file.
    -n, --num <integer>         Optional number to create a random seed, useful to run several
                                instances of the sampler in parallel. The number is also the
                                ROOT output filename.
    -o, --output <directory>    Optional parameter to overwrite the output directory given
                                by "output_dir" in the config file.
    -s, --surface <file>        Optional parameter to overwrite the freeze-out hypersurface
                                file given by "surface_file" in the config file.
    -q, --quiet                 Suppress the disclaimer and config parameters print-out.
    --version                   Print version of the sampler executable.

Example usage:

    ./sampler -c /path/to/config-example --num 1 --output /path/to/output/dir -s /path/to/freezeout.dat


## Config file

The following lists **all possible config parameters**:

Mandatory parameters:

    surface_file                  Path to the freeze-out hypersurface file that gets sampled.
    output_dir                    Path to the output directory.
    number_of_events              Number of events that are sampled.
    ecrit                         Critical energy density at which the hydro stopped in a particular
                                  cell and the freeze-out hypersurface was constructed.


Optional parameters:

    bulk                          Enables bulk viscosity if set to 1.   Default is 0 (false).
    shear                         Enables shear viscosity if set to 1.  Default is 0 (false).
    cs2                           Speed of sound squared.               Default is 0.15.
    ratio_pressure_energydensity  Pressure divided by energy density.   Default is 0.15.
    create_root_output            Enables ROOT output if set to 1.      Default is 0 (false).
    hydro_coordinate_system       Coordinate system in which the position of the freeze-out
                                  hypersurface elements from the hydro evolution are provided
                                  (by space-time four-vectors). Possible coordinate systems
                                  are 'tau-eta' (default) or 'cartesian'.


> [!NOTE]
> Lines in the config file can be commented out by using an exclamation point in the beginning of a line (e.g., `! This is a comment`).

> [!TIP]
> The repository provides an example config file named `config-example`.

In case this example config file is used, the `surface_file` parameter has to be set to the location of the freeze-out file and the `output_dir` parameter to the desired output path, either in the config file itself (instead of _/path/to/freezeout/file_ and _/output/path_) or via the command line parameters mentioned above!
