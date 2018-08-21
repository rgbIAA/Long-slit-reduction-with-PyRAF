# This is a tutorial example for the introductory course "Long-slit reduction with PyRAF"
# given at the "School on Long-Slit Spectroscopy" at University of Rwanda from August 26th 
# to September 1st of 2018:
#
# https://sites.google.com/site/schoolonspectroscopy/home
#
# Author: Ruben Garcia-Benito (RGB)
#
# This are the long-slit reducing steps for the tutorial
# Requieres the "relspyraf.py" file (module)

from collections import OrderedDict
from astropy.table import Table
from relspyraf import *

# List of tasks to use to apply unlearn ---------------------------------
ltask = ['ccdproc','zerocombine','imstatistics','flatcombine','response',
	'identify','reidentify','fitcoords','transform',
	'imcombine','lacos_spec','background','apall',
	'standard','sensfunc','calibrate','setairmass']

# Observation information for the tutorial at University of Rwanda
# 2007/07/08 - ISIS@WHT
# Dispersion [A/pix] = 0.45 (B) 
# Spatial Resolution =  0.2 (arcsec/pixel)
# Slit width = 1.''
# GAIN_R2 = 1.17 [e-/ADU]
# RDNOISE_R2 = 4.0 [e-]
# ARC LAMPS in Blue ARM: CuNe+CUAR

# Parameters ------------------------------------------------------------
# Change your parameters according to your observations
rdnoise = 4.0    
gain    = 1.17

# Bias & Trim from header
biassec = '[316:1201,4105:4190]'
trimsec = '[1:1201,50:4000]'

# Files and directories for IRAF
database = "database"             # Database for wavelength calibration and apertures
coor     = 'linelists$cuar.dat'   # List of arc lines

# Measure the FWHM in a reduced arc with splot 'k' (file.fits[1][1000,*])
# For IDENTIFY
fwidth = 3.4

# Flux calibration
observatory = 'lapalma'
extinction  = 'lapalma_extinction.dat'

# Names files --------------------------------------------------
bias  = 'bias.fits'   # Master BIAS
flat  = 'flat.fits'   # Master FLAT

# Name lists ---------------------------------------------------
fbias = 'bias.list'   # Bias list
fflat = 'flat.list'   # Flat list

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# ------------------------------------------------------------------------
# --------------------------- REDUCCION ----------------------------------
# ------------------------------------------------------------------------

# Read LOG table. Prepare an organized table with a column with 'object' 
# with unique description for bias, flat, etc

tlog = Table.read('nightlog_2005-07-08.txt', format='ascii')

# Unlearn list of tasks
unlearn_list(ltask)

# ---------------------------- BIAS ------------------------------------

# Some fits files contain their data in the second HDU. You need to add "[1]" to 
# the name of the files for IRAF to read it

#idbias = (tlog['object'] == 'bias')
#wfile(fbias, tlog[idbias]['file'], suffix='[1]')

# Check if there is something wrong with an image
#imstatistics('@'+fbias, fields='image,mean,midpt,stddev,min,max')

#delfiles(bias)
# Creates an image centered in 0 of pixel to pixel variation (overscan removed)
#zerocombine('@'+fbias, output=bias, combine='average', ccdtype="", rdnoise=rdnoise, gain=gain)


# ---------------------------- FLATS ------------------------------------

#idflat = (tlog['object'] == 'flat')
#wfile(fflat, tlog[idflat]['file'], suffix='[1]')

#imstatistics('@'+fflat, fields='image,mean,midpt,stddev,min,max')

## FLATCOMBINE, CCDPROC y RESPONSE de FLATS ------------

##flatccdres(fflat, lflats=tlog[idflat]['file'], suffix='[1]', rdnoise=rdnoise, gain=gain, 
##           biassec=biassec, trimsec=trimsec, zero=bias)
#flatccdres(fflat, rdnoise=rdnoise, gain=gain, biassec=biassec, trimsec=trimsec, zero=bias)


# ------------------ Apply CCDPROC to all IMAGES: OBJETs, STANDARDs, ARCs --------

# Select all images except bias and flats
#idim = (tlog['object'] != 'bias') & (tlog['object'] != 'flat')# & (tlog['object'] != 'arc')

# Apply ccdproc
#ccdp(tlog['file'][idim], flat=flat, zero=bias, trimsec=trimsec, biassec=biassec, suffix='[1]')

# ------------ Combine science and standard star images to remove cosmic rays --------------

idobj  = (tlog['object'] == 'GALAXY')
imgal  = 'galaxy.fits'

idstd1 = (tlog['object'] == 'bd33')
std1   = 'bd33.fits'

idstd2 = (tlog['object'] == 'bd25')
std2   = 'bd25.fits'

#for selection, name in zip([idobj, idstd1, idstd2], [imgal, std1, std2]):
#    combine(tlog['file'][selection], name, prefix='c', rdnoise=rdnoise, gain=gain)

# Remove cosmic rays if there is only one image per object/standard star ------

#rmcosmic(tlog['file'][idobj], fmt_in='c%s', fmt_out='clean%s', rdnoise=rdnoise, gain=gain)


# ---------- WAVELENGTH CALIBRATION: IDENTIFY, REIDENTIFY & FITCOORD ------------

idarc  = (tlog['object'] == 'arc')
arcref = 'c%s[1]' % tlog['file'][idarc][0]

# Find solution in one arc file
#wlsolution(arcref, fwidth=fwidth)

# --------- Apply wavelength solution to science images: TRANSFORM -------------
#idobj = (tlog['object'] == 'GALAXY')
#lobj = ['c%s' % fits for fits in tlog['file'][idobj]]
lobj     = [imgal, std1, std2]
fitnames = arcref

#transf(lobj, fitnames)


# --------------------------- SETAIRMASS ---------------------------------------

#setair(lobj, observatory=observatory)


# ---------------------- STANDARD STAR EXTRACTION ------------------------------

lstd = OrderedDict()
lstd['lbd25.fits'] = dict(line=617, lower=-50, upper=50, b_sample='-300:-200,200:300')
lstd['lbd33.fits'] = dict(line=642, lower=-50, upper=50, b_sample='-300:-200,200:300')

#for std in lstd:
#    extract_aperture(std, std=True, rdnoise=rdnoise, gain=gain, **lstd[std])


# ---------------------------- SENSITIVITY FUNCTION ----------------------------

dsens = OrderedDict()
dsens['aplbd25.0001.fits'] = dict(star_names='mbd25d4655', sensitivity='sensf_bd25') 
dsens['aplbd33.0001.fits'] = dict(star_names='mbd33d2642', sensitivity='sensf_bd33') 
#dsens['aplbd25.0001.fits'] = dict(star_names='bd253941', sensitivity='sensf_bd25')
#dsens['aplbd33.0001.fits'] = dict(star_names='bd332642', sensitivity='sensf_bd33')

##for std in dsens:
##    sensf(std, **dsens[std])

lstd       = dsens.keys()
star_names = [dsens[std]['star_names'] for std in dsens]

# caldir="" if using external file 'mbd25d4655' and 'mbd33d2642'
#sensf(lstd, star_names=star_names, caldir="", observatory=observatory, extinction=extinction)

## In case you have different sensitivity functions and want to combine them
## iraf.stsdas.toolbox.imgtools(_doprint=0)
##imcalc('sens1.fits,sens2.fits', 'sensfunc.fits', '(im1+im2)/2.')

# ------------------------- Remove SKY (background) ------------------------

sample = '200:400,900:1100'
# Output of "back" adds a "b" to the input file by default (use variable "outfits" otherwise)
#back('lgalaxy.fits', sample=sample)


# ---------------------------- Flux calibration ----------------------------

## Output of "calflux" adds a "cf" to the input file by default (use variable "outfits" otherwise)
#calflux('blgalaxy.fits', observatory=observatory, extinction=extinction)


# ------------------------- 1D OBJECT EXTRACTION ---------------------------

galfits = 'cfblgalaxy.fits'

lgal = OrderedDict()
# Define your initial appertures (you can change them interactively)
lgal['galaxy_knotb.fits'] = dict(line=3670, lower=-24, upper=21)

#for gal in lgal:
#    extract_aperture(galfits, outfits=gal, find='no', rdnoise=rdnoise, gain=gain, **lgal[gal])


# ---------------------------------------------------------------
# Final unlearn
unlearn_list(ltask)
