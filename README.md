# Long-slit-reduction-with-PyRAF
Long slit reduction with PyRAF

This is a tutorial example for the introductory course "Long-slit reduction with PyRAF" given at the ["School on Long-Slit Spectroscopy"](https://sites.google.com/site/schoolonspectroscopy/home) at University of Rwanda from August 26th to September 1st of 2018.

Author: Rubén García-Benito (RGB)

The FITS files used are not included in this repository. 

The functions and scripts are meant to be adapted (rewritten) to your paticular circumstances.

# Dependences

The tutorial require these packages to be installed in order to run:

+ numpy
+ astropy
+ IRAF/PyRAF

For IRAF/PyRAF I recommend the [astroconda](https://astroconda.readthedocs.io/en/latest/) package (with [astroconda with PyRAF
](https://astroconda.readthedocs.io/en/latest/installation.html#legacy-software-stack-with-iraf))

# Content

* "Long slit PyRAF reduction I - Rwanda 2018 - RGB.pdf" and "Long slit PyRAF reduction II - Rwanda 2018 - RGB.pdf" are the slides of the reduction tutorial given and the [School](https://sites.google.com/site/schoolonspectroscopy/home).

* "relspyraf.py": container of the functions needed for reduction ("Reduce Long Slit Pyraf" module). Put it in the same directory of "long_slit_reduction_with_PyRAF.py" or in some of your "PYTHONPATH" directories

* "long_slit_reduction_with_PyRAF.py": logical sequence of reduction steps 

* "data" directory: additional files used in the tutorial for the reduction examples (except the FITS files)
   * "nightlog_2005-07-08.txt": log file of the observations used in the tutorial
   * "lapalma_extinction.dat": Extinction file of [Roque de los Muchachos Observatory](http://www.ing.iac.es/Astronomy/observing/conditions/wlext.html)
   * "standards.men": config file for IRAF "standard" & "sensfunc"
   * "mbd25d4655.dat" and "mbd33d2642.dat": Standard stars calibration files (from ESO: ftp://ftp.eso.org/pub/stecf/standards/okestan/)


# Instructions

Put "relspyraf.py" in the same directory of "long_slit_reduction_with_PyRAF.py" or in some of your "PYTHONPATH" directories.

Browse "long_slit_reduction_with_PyRAF.py" to be familiar with the reduction process and change accordingly to your needs. You can also modify "relspyraf.py" to include some more functionalities or new functions. 

This is an example of a standard (easy) reduction. There might be additional steps required by your observations not included in this tutorial (e.g. illumination correction)
