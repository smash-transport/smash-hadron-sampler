# Hadron Sampler

This sampler is meant to be applied in hybrid models, simulating heavy-ion collisions in the high temperature or baryon-density region. More precisely it provides an interface between the macroscopic hydrodynamic evolution of the fireball and the hadronic afterburner. During the hydrodynamic evolution a hypersurface of constant energy density (the switching energy density) is created. Each element on the hypersurface needs then be transformed into a list of particles with properties loosely provided by the macroscopic properties of the hypersurface elements. This process of particlization is performed by means of the hadron sampler provided within this project. It is designed to couple the [3+1D viscous hydrodynamic code vhlle](https://github.com/yukarpenko/vhlle) to the [hadronic transport model SMASH](https://smash-transport.github.io). For details about the sampling algorithm, please consult [I. Karpenko et al., Phys.Rev.C 91 (2015) 6, 064901](https://inspirehep.net/literature/1343339).

When using the smash-hadron-sampler, please cite:
- [I. Karpenko et al., Phys.Rev.C 91 (2015) 6, 064901](https://inspirehep.net/literature/1343339)
- [A. Schäfer el al., arXiv:2112.08724](https://arxiv.org/abs/2112.08724).

### Prerequisites:
- [cmake](https://cmake.org) version &ge; 3.15.4
- [SMASH](https://github.com/smash-transport/smash) version 3.0, as well as prerequisites therein
- [ROOT](https://root.cern.ch) version &ge; 6.06

Please note that only tagged versions are guaranteed to be compatible with SMASH.

### Install instructions:
It is expected that the output of this sampler is used in combination with the SMASH transport model. We therefore assume SMASH was already compiled and is available to be used as an external library. All necessary prerequisites are also assumed to already be installed.
If not, install instructions can be found [here](https://github.com/smash-transport/smash-devel/blob/master/README.md).

To compile the project, first set the environment variable to the smash directory:

    export SMASH_DIR=[...]/smash

Copy the cmake files to the sampler directory:

    cd [...]/smash-hadron-sampler
    cp -r $SMASH_DIR/cmake ./

Execute the following commands to build the project:

    mkdir build
    cd build
    cmake .. -DPythia_CONFIG_EXECUTABLE=[...]/pythia8309/bin/pythia8-config
    make
where `[...]/pythia8309` is the path to the pythia directory to which also SMASH is coupled.

In continuation, the executable `sampler` is created.


### Execute the sampler
To run the sampler, execute the following command:

    ./sampler events NUM PATH_TO_CONFIG_FILE

where `NUM` is a random number set by the user. It can be useful to run several instances of the sampler in parallel. `PATH_TO_CONFIG_FILE` provides the path to the configuration file. Therein the location of the freezeout hypersurface file, the path to the output directory and all other necessary parameters are set.
