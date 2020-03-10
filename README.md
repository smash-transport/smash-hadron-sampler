# Hadron Sampler

This sampler is meant to be applied in hybrid models, simulating heavy-ion collisions in the high temperature or baryon-density region. More precisely it provides an interface between the macroscopic hydrodynamic evolution of the fireball and the hadronic afterburner. During the hydrodynamic evolution a hypersurface of constant energy density (the switching energy density) is created. Each element on the hypersurface needs then be transformed into a list of particles with properties loosely provided by the macroscopic properties of the hypersurface elements. This process of particlization is performed by means of the hadron sampler provided within this project. It is designed to couple the [3+1D viscous hydrodynamic code vhlle](https://github.com/yukarpenko/vhlle) to the [hadronic transport model SMASH](https://smash-transport.github.io).

#### Install instructions:
It is expected that the output of this sampler is used in combination with the SMASH transport model. We therefore assume SMASH was already compiled and is available to be used as an external library. All necessary prerequisites are also assumed to already be installed.
If not, install instructions can be found [here](https://github.com/smash-transport/smash-devel/blob/master/README.md).

To compile the project, first set the environment variable to the smash directory:

    export SMASH_DIR=[...]/smash

Copy the cmake files to the sampler directory:

    cd [...]/hadronSampler
    cp -r $SMASH_DIR/cmake ./

Execute the following commands to build the project:

    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=[...]/eigen3 -DPythia_CONFIG_EXECUTABLE=[...]/pythia8235/bin/pythia8-config
    make

In continuation, the executable `sampler` is created.


#### Execute the sampler
To run the sampler, execute the following command:

    ./sampler events NUM PATH_TO_CONFIG_FILE

where `NUM` is a random number set by the user. It can be useful to run several instances of the sampler in parallel. `PATH_TO_CONFIG_FILE` provides the path to the configuration file. Therein the location of the freezeout hypersurface file, the path to the output directory and all other necessary parameters are set.