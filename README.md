# Hadron Sampler

This sampler is meant to be applied in hybrid models, simulating heavy-ion collisions in the high temperature or baryon-density region. More precisely it provides an interface between the macroscopic hydrodynamic evolution of the fireball and the hadronic afterburner. During the hydrodynamic evolution a hypersurface of constant energy density (the switching energy density) is created. Each element on the hypersurface needs then be transformed into a list of particles with properties loosely provided by the macroscopic properties of the hypersurface elements. This process of particlization is performed by means of the hadron sampler provided within this project. It is designed to couple the [3+1D viscous hydrodynamic code vhlle](https://github.com/yukarpenko/vhlle) to the [hadronic transport model SMASH](https://smash-transport.github.io). For details about the sampling algorithm, please consult [I. Karpenko et al., Phys.Rev.C 91 (2015) 6, 064901](https://inspirehep.net/literature/1343339).

When using the smash-hadron-sampler, please cite:
- [I. Karpenko et al., Phys.Rev.C 91 (2015) 6, 064901](https://inspirehep.net/literature/1343339)
- [A. Schäfer et al., arXiv:2112.08724](https://arxiv.org/abs/2112.08724).


### Prerequisites:
- [cmake](https://cmake.org) version &ge; 3.15.4
- [SMASH](https://github.com/smash-transport/smash) version 3.2, as well as prerequisites therein
- [ROOT](https://root.cern.ch) version &ge; 6.06

Please note that only tagged versions are guaranteed to be compatible with SMASH.


### Install instructions:
It is expected that the output of this sampler is used in combination with the SMASH transport model. We therefore assume SMASH was already compiled and is available to be used as an external library. All necessary prerequisites are also assumed to already be installed.
If not, install instructions can be found [here](https://github.com/smash-transport/smash/blob/main/README.md).

To compile the project, first set the environment variable to the smash directory:

    export SMASH_DIR=[...]/smash

Copy the cmake files to the sampler directory:

    cd [...]/smash-hadron-sampler
    cp -r $SMASH_DIR/cmake ./

Execute the following commands to build the project:

    mkdir build
    cd build
    cmake .. -DPythia_CONFIG_EXECUTABLE=[...]/pythia8312/bin/pythia8-config
    make
where `[...]/pythia8312` is the path to the pythia directory to which SMASH is also coupled.

In continuation, the executable `sampler` is created.


### Execute the sampler
To run the sampler, execute the following command:

    ./sampler events NUM PATH_TO_CONFIG_FILE

where `NUM` is a random number set by the user. It can be useful to run several instances of the sampler in parallel. `PATH_TO_CONFIG_FILE` provides the path to the configuration file. Therein, the location of the freezeout hypersurface file, the path to the output directory, and all other necessary parameters can be specified.


### Config file
The repository provides a config file `config-example`.  
:warning: Attention: In case this config file is used, the `surface` parameter has to be set to the location of the freezeout file (instead of _/path/to/freezeout/file_) and the `spectra_dir` parameter to the desired output path (instead of _/output/path_)!  

The following lists **all possible config parameters** (to read some explanations entirely scroll to the right):

Mandatory parameters:
```
surface                       Path to the freezeout hypersurface file that gets sampled.
spectra_dir                   Path to the output directory.
number_of_events              Number of events that are sampled.
ecrit                         Critical energy density at which the hydro stopped in a particular cell and the freezeout hypersurface was constructed.
```

Optional parameters:
```
bulk                          Enables bulk viscosity if set to 1.   Default is 0 (false).
shear                         Enables shear viscosity if set to 1.  Default is 0 (false).
cs2                           Velocity of sound squared.            Default is 0.15.
ratio_pressure_energydensity  Pressure divided by energy density.   Default is 0.15.
createRootOutput              Enables ROOT output if set to 1.      Default is 0 (false).
```
