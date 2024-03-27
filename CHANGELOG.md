# Changelog


All notable changes to this project will be documented in this file. The format is inspired by [Keep a Changelog](https://keepachangelog.com/en/1.0.0/). This project follows the versions of [SMASH](https://github.com/smash-transport/smash), that means that each new version of SMASH triggers a new version of this project, in order to assure compatibility.
This changelog is in place since version SMASH-hadron-sampler-3.0.

The major categories to group changes in this log are:

* :left_right_arrow: for all, in particular breaking, changes, fixes and additions to the in- and output files.
* :heavy_plus_sign: for new features.
* :recycle: for changes in existing functionality.
* :sos: for any bug fixes.
* :heavy_minus_sign: for now removed features.

Also possible, but for this project less relevant, is `Deprecated` for soon-to-be removed features.


## Unreleased

* :left_right_arrow: ROOT output is now disabled by default, it can be enabled in the config by setting createRootOutput parameter to 1

## SMASH-hadron-sampler-3.1
Date: 2023-02-29

* :heavy_plus_sign: Preliminary support for delta f corrections for bulk
* :sos: Adjusted compiler flags for Apple machines

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-hadron-sampler-3.0...SMASH-hadron-sampler-3.1)

## SMASH-hadron-sampler-3.0
Date: 2023-04-28

* :sos: The particle buffer was increased in order to sample from a bigger hypersurface
* :recycle: ⚠️ The `master` branch has been renamed to `main`
* :recycle: The hadron sampler makes now use of the C++17 standard

[Link to diff from previous version](https://github.com/smash-transport/smash/compare/SMASH-hadron-sampler-1.1...SMASH-hadron-sampler-3.0)

## SMASH-hadron-sampler-1.1
Date: 2022-05-17
