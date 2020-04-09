## Resubmission
This is a resubmission. In this version:

* Fixed a minor bug introduced by a new feature in 1.0.1
* Changed version number from 1.0.1 to 1.0.2


## Test environments
* local MacOS 10.15.4, R 3.6.3
* Ubuntu 16.04.6 LTS (on travis-ci), R 3.6.2
* win-builder (devel, release, old-release)


## R CMD check results
There were no ERRORs or WARNINGs.

Rhub and win-builder have a note for the days since last update: 0. This is a 
fix for a small bug identified after the previous submission yesterday.

There is a WARN for v1.0.0 of the package on r-patched-osx-x86_64. The warning 
appears to be related to an issue with a missing driver in the pandoc 
installation. There was no error for this check when v1.0.1 was submitted, and
the vignette causing the error hasn't been changed.


## Downstream dependencies
saotd: 1 Note: found 826 marked UTF-8 strings -- This is not unexpected for a 
  package that contains datasets for natural language processing.


