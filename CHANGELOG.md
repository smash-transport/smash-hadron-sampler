# Changelog


All notable changes to this project will be documented in this file.
The format is inspired by [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
This project follows the versions of [SMASH](https://github.com/smash-transport/smash), that means that each new version of SMASH triggers a new version of this project, in order to assure compatibility.
For intermediate releases of the SMASH-hadron-sampler (not induced by a new SMASH release) the third digit of the versioning number will be bumped, which is therefore not solely reserved for bug fixes.
A [release checklist](#release-checklist) is provided at the end of the file.
The described procedure should be followed for every new release.

This changelog is in place since version SMASH-hadron-sampler-3.0.

The major categories to group changes in this log are:

* :left_right_arrow: for all, in particular breaking, changes, fixes and additions to the in- and output files.
* :heavy_plus_sign: for new features.
* :recycle: for changes in existing functionality.
* :sos: for any bug fixes.
* :heavy_minus_sign: for now removed features.

Also possible, but for this project less relevant, is `Deprecated` for soon-to-be removed features.


## Unreleased

* :left_right_arrow: Renamed config key `surface`, used to specify the freezeout surface file, to `surface_file` to clarify key.
* :left_right_arrow: Renamed config key `createRootOutput`, used to enable ROOT output, to `create_root_output` to unify naming convention.
* :left_right_arrow: Renamed config key `spectra_dir`, used to specify the output directory, to `output_dir`.
* :left_right_arrow: Changed the command line parameters used to run the sampler to clarify their usage.
* :heavy_plus_sign: Added optional command line parameters to overwrite the `surface` and `output_dir` keys in the config file.
* :sos: Fix problems with thread safety from ROOT objects
* :left_right_arrow: ROOT output is now disabled by default, it can be enabled in the config by setting createRootOutput parameter to 1


## SMASH-hadron-sampler-3.1.1
Date: 2025-02-17

* :heavy_plus_sign: Fix that introduces the command line option `./sampler --version` to get the version of the sampler executable

[Link to diff from previous version](https://github.com/smash-transport/smash-hadron-sampler/compare/SMASH-hadron-sampler-3.1...SMASH-hadron-sampler-3.1.1)


## SMASH-hadron-sampler-3.1
Date: 2023-02-29

* :heavy_plus_sign: Preliminary support for delta f corrections for bulk
* :sos: Adjusted compiler flags for Apple machines

[Link to diff from previous version](https://github.com/smash-transport/smash-hadron-sampler/compare/SMASH-hadron-sampler-3.0...SMASH-hadron-sampler-3.1)


## SMASH-hadron-sampler-3.0
Date: 2023-04-28

* :sos: The particle buffer was increased in order to sample from a bigger hypersurface
* :recycle: ⚠️ The `master` branch has been renamed to `main`
* :recycle: The hadron sampler makes now use of the C++17 standard

[Link to diff from previous version](https://github.com/smash-transport/smash-hadron-sampler/compare/SMASH-hadron-sampler-1.1...SMASH-hadron-sampler-3.0)


## SMASH-hadron-sampler-1.1
Date: 2022-05-17


## Release checklist

1. Create a `release` branch from `develop` and switch to it.
2. Adapt the _CMakeLists.txt_ file:
  * bump the hard coded version macro (`VERSION` argument of `project(...)` command);
  * comment out `set(SAMPLER_VERSION "${SAMPLER_VERSION}-next")` with a `#` to mark the codebase as clean.
3. Adjust the _CHANGELOG.md_ file:
  * add `## SMASH-hadron-sampler-X.Y.Z` below the `## Unreleased` line (leaving one empty line in-between);
  * prepare a `Date:` line below the release candidate section, which will be completed on the day of the public release, as well as the _Link to diff from previous version_ at the end of the upcoming release block (this will become a meaningful link after the public release).
4. If the new release is due to a new SMASH release, bump the mentioned SMASH version in the _README.md_ file section `Prerequisites`.
5. Close the `release` branch in the git-flow sense:
  * merge it into the `main` branch;
  * switch to `main` and tag the last commit;
  * switch to `develop` and merge the `release` branch back into it.
6. Publish the new release by pushing the changes and the new tag on `main`.
7. On branch `develop`:
  * uncomment `set(SAMPLER_VERSION "${SAMPLER_VERSION}-next")` to mark the codebase as dirty (remove the `#`);
  * push the changes on `develop`.
8. Remove the `release` branch.

This checklist is inspired by the [release checklist of the SMASH-vHLLE-Hybrid](https://smash-transport.github.io/smash-vhlle-hybrid/latest/developer/release_procedure/#release-checklist).
