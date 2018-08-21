# This is a tutorial example for the introductory course "Long-slit reduction with PyRAF"
# given at the "School on Long-Slit Spectroscopy" at University of Rwanda from August 26th 
# to September 1st of 2018:
#
# https://sites.google.com/site/schoolonspectroscopy/home
#
# Author: Ruben Garcia-Benito (RGB)

# %%%%%%%%%%%%%%%% IMPORTS %%%%%%%%%%%%%%%%%%%
# PYRAF import --------------------
from pyraf import iraf
iraf.noao.imred(_doprint=0)
iraf.noao.twodspec(_doprint=0)
iraf.noao.onedspec(_doprint=0)
iraf.noao.imred.ccdred(_doprint=0)
iraf.noao.twodspec.longslit(_doprint=0)
iraf.noao.twodspec.apextract(_doprint=0)
iraf.stsdas.toolbox(_doprint=0)
iraf.stsdas.toolbox.imgtools(_doprint=0)
from iraf import ccdproc, zerocombine, imstatistics
from iraf import flatcombine, response, hedit, unlearn
from iraf import identify, reidentify, fitcoords, transform
from iraf import imcombine, lacos_spec, background, apall
from iraf import standard, sensfunc, calibrate
from iraf import setairmass, imcalc
# -------------------------------

from collections import OrderedDict
from astropy.table import Table
from glob import glob
import numpy as np
import astropy
import os
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# List of tasks to use to apply unlearn ---------------------------------
ltask = ['ccdproc','zerocombine','imstatistics','flatcombine','response',
	'identify','reidentify','fitcoords','transform',
	'imcombine','lacos_spec','background','apall',
	'standard','sensfunc','calibrate','setairmass','imcalc']

# Observation information for the tutorial at University of Rwanda
# 2007/07/08 - ISIS@WHT
# Dispersion [A/pix] = 0.45 (B) 
# Spatial Resolution =  0.2 (arcsec/pixel)
# Slit width = 1.''
# GAIN_R2 = 1.17 [e-/ADU]
# RDNOISE_R2 = 4.0 [e-]
# ARC LAMPS in Blue ARM: CuNe+CUAR

# Parameters ------------------------------------------------------------
# Change your parameters according to your observations. "rdnoise" and "gain" are 
# the only global parameters in this tutorial
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

# Names files --------------------------------------------------
bias  = 'bias.fits'   # Master BIAS
flat  = 'flat.fits'   # Master FLAT

# Name lists ---------------------------------------------------
fbias = 'bias.list'   # Bias list
fflat = 'flat.list'   # Flat list

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ---------------------- FUNCTIONS -----------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Check if object is list --------------------------------------
def cklist(obj):
    if not isinstance(obj, (list, tuple, np.ndarray, astropy.table.table.Column)):
        obj = [obj]
    return obj

# Function to write file from a list/array ---------------------
def wfile(nfile, values, suffix=None, prefix=None):
    suffix = '' if suffix is None else suffix
    prefix = '' if prefix is None else prefix
    values = cklist(values)
    values = ['%s%s%s' % (prefix, item, suffix) for item in values]
    np.savetxt(nfile, values, fmt='%s')

# Delete list of files -----------------------------------------
def delfiles(files, directory=None):
    files = cklist(files)
    for nfile in files:
        if directory is not None:
            nfile = os.path.join(directory, nfile)
        if nfile is not None and os.path.exists(nfile):
            os.remove(nfile)

# Change name of file with a given format ----------------------
def fmt_file(fmt, name):
    dirname  = os.path.dirname(name)
    rname    = os.path.basename(name)
    new_name = os.path.join(dirname, fmt % rname)
    return new_name

# Unlearn list of functions ------------------------------------
def unlearn_list(ltask):
    ltask = cklist(ltask)
    for task in ltask: 
        unlearn(task)

# Master flat --------------------------------------------------
def flatccdres(flat_files, flat='flat.fits', nflat='nflat.fits', cflat='cflat.fits', iccd='no', iresponse='yes', 
	combine='average', reject='avsigclip', ccdtype="", process='no', fixpix='no', overscan='yes', 
	zerocor='yes', trim='yes', darkcor='no', flatcor='no', illumco='no', fringec='no', function='legendre', 
	corder=5, rorder=18, readaxi='column', biassec=biassec, zero=bias, trimsec=trimsec, low_reject=3,
	high_reject=3, scale='mode', lflats=None, suffix=None):

    unlearn('ccdproc')
    unlearn('flatcombine')

    delfiles(flat)
    if lflats is not None:
        wfile(flat_files, lflats, suffix=suffix)

    flatcombine('@'+flat_files, output=flat, combine=combine, reject=reject, ccdtype=ccdtype,
	scale=scale, process=process, rdnoise=rdnoise, gain=gain)

    delfiles(cflat)
    ccdproc(images=flat, output=cflat, ccdtype=ccdtype, fixpix=fixpix, overscan=overscan, trim=trim, zerocor=zerocor,
	darkcor=darkcor, flatcor=flatcor, illumco=illumco, fringec=fringec, zero=zero, function=function, 
	order=corder, interactive=iccd, readaxi=readaxi, trimsec=trimsec, biassec=biassec)

    # Dispersion goes up-down
    # hedit(cfout,'DISPAXIS', 2, add='yes', verify='no')
    delfiles(nflat)
    response(cflat, cflat, nflat, interac=iresponse, order=rorder, low_reject=low_reject, high_reject=high_reject)

    unlearn('ccdproc')
    unlearn('flatcombine')

# ccdproc function  --------------------------------------------
def ccdp(files, suffix=None, input_files='tmp_input.list', output_files='tmp_output.list', **kwargs):

    # Dictionary with default parameters for ccdproc
    dpar = dict(ccdtype="", fixpix='no', overscan='yes', trim='yes', zerocor='yes', flatcor='no', darkcor='no',
	        illumco='no', fringec='no', zero="", illum="", function='legendre', order=4, 
                interactive='no', readaxi='column', trimsec=None, biassec=None, flat="")

    # Update parameters dictionary for ccdproc if new are provided (kwargs)
    dpar.update(kwargs)

    # Make "files" a list, in case there is only one file
    files = cklist(files)

    # Write the list of input files to a temporary file
    wfile(input_files, files, suffix=suffix)

    # Create a list of names with the output name after ccdproc is applied
    new_files_ccdproc = ['c%s' % item for item in files]
    # Write the new files list to a file
    wfile(output_files, new_files_ccdproc)

    # Delete, if exits, the new files to be created
    delfiles(new_files_ccdproc)

    # Apply ccdproc 
    unlearn('ccdproc')
    ccdproc(images='@'+input_files, output='@'+output_files, **dpar)
    unlearn('ccdproc')

    # Delete the temporary lists
    delfiles([input_files, output_files])

# Identify (Wavelength calibration solution) ---------------------
def wlsolution(arc, arc_ref=None, iden=True, reiden=True, section="middle column", coor='linelists$cuar.dat',
               database='database', fwidth=3.4, function='legendre', order=4, nsum=10, overrid='yes', refit='no',
               trace='yes', nlost=20, step=10, xorder=3, yorder=4, rinterac='no', finterac='yes', 
               fitc=True, suffix=None):
    # twodspec.longslit
    # Measure the FWHM (fwidth) in a reduced arc with splot 'k' (file.fits[1][1000,*])

    lw_tasks = ['identify', 'reidentify', 'fitcoords']
    unlearn_list(lw_tasks)

    if suffix is not None:
        arc = '%s%s' % (arc, suffix)

    if arc_ref is None: 
        arc_ref = arc

    arc_ref = arc_ref.split('[')[0]

    if iden:
        # 'm': select lines | 'l': searching in a list | 'f': final fit
        identify(arc, section=section, database=database, fwidth=fwidth, coordli=coor, 
                 function=function, order=order, nsum=nsum)

    if reiden:
        reidentify(images=arc, reference=arc_ref, interac=rinterac, section=section, overrid=overrid, 
                  refit=refit, trace=trace, nlost=nlost, step=step, nsum=nsum, coordlist=coor)

    if fitc:
        arc = '.'.join(arc.split('.')[0:-1])
        fitcoords(arc, fitname="", interac=finterac, xorder=xorder, yorder=yorder, function=function)

    unlearn_list(lw_tasks)

# Apply Wavelength solution ------------------------------------
def transf(lfits, fitnames, interp='linear', database='database', flux='yes', y1='INDEF', suffix=None,
           input_files='tmp_input.list', output_files='tmp_output.list', dispaxis=None, **kwargs):

    fitnames = '.'.join(fitnames.split('.')[0:-1])

    if isinstance(lfits, str) and input_files is None:
        input_files = lfits
        lfits = np.loadtxt(input_files, dtype=np.str, unpack=True)
    else:
        lfits = cklist(lfits)
        wfile(input_files, lfits, suffix=suffix)

    wfits = ['l%s' % fits.split('[')[0] for fits in lfits]
    wfile(output_files, wfits)
    delfiles(wfits)

    if dispaxis is not None:
        for fits in lfits:
            hedit(fits, 'DISPAXIS', dispaxis, add='yes', verify='no') 

    unlearn('transform')
    transform(input='@'+input_files, output='@'+output_files, fitnames=fitnames, interp=interp,
              database=database, flux=flux, y1=y1, **kwargs)
    unlearn('transform')

    # for fits in wfits:
    #     hedit(fit, 'WAT2_001', 'wtype=linear label=Wavelength units=Angstroms', add='yes', verify='no') 

    delfiles([input_files, output_files])

# Combine images -----------------------------------------------
def combine(lfits, outfits, rej='crreject', combine='average', scale='exposure', weight='exposure', nlow=1,
            expname='EXPTIME', suffix=None, prefix=None, input_files='tmp_input.list', delete_list=True):

    wfile(input_files, lfits, suffix=suffix, prefix=prefix)
    delfiles(outfits)

    unlearn('combine')
    imcombine('@'+input_files, output=outfits, combine=combine, reject=rej, scale=scale,
		weight=weight, expname=expname, rdnoise=rdnoise, gain=gain, nlow=nlow)
    unlearn('combine')

    if delete_list:
        delfiles(input_files)

# Remove cosmic rays of single images --------------------------
def rmcosmic(lfits, xorder=5, yorder=5, sigfrac=2.0, niter=5, objlim=2, sigclip=3.0, fmt_in=None, fmt_out='cr%s'):

    lfits = cklist(lfits)
    fmt_in = '%s' if fmt_in is None else fmt_in

    for fits in lfits:
        fits   = fmt_in  % fits
        output = fmt_out % fits
        mask   = 'mask%s' % fits
        delfiles([mask, output])

        unlearn('lacos_spec')
        lacos_spec(fits, output, mask, gain=gain, readn=rdnoise, xorder=xorder, yorder=yorder,
                   sigclip=sigclip, sigfrac=sigfrac, objlim=objlim, niter=niter)
        unlearn('lacos_spec')

        delfiles(mask)

# Set airmass in images ----------------------------------------
def setair(files, observa='lapalma', equinox='EPOCH', ut='UTSTART', input_files='tmp_input.list'):

    files = cklist(files)
    wfile(input_files, files)

    unlearn('setairmass')
    setairmass(images='@'+input_files, observa=observa, equinox=equinox, ut=ut)
    unlearn('setairmass')

    delfiles(input_files)


# Extract aperture ---------------------------------------------
def extract_aperture(fits, line, lower=-5, upper=5, b_sample='-10:-6,6:10',
    outfits=None, fmt_out='ap%s', fmt_ap='%s.000?.fits',
    database='database/', std=False, **kwargs):

    # line, nsum: The dispersion line (line or column perpendicular to the dispersion axis) and 
    # number of adjacent lines (half before and half after unless at the end of the image) 
    # used in finding, recentering, resizing, and editing operations. A line of INDEF selects 
    # the middle of the image along the dispersion axis. A positive nsum selects a sum of 
    # lines and a negative selects a median of lines. 

    # lower, upper: Default lower and upper aperture limits relative to the aperture center. 
    # These limits are used for apertures found with apfind and when defining the first 
    # aperture in apedit. 

    # b_sample = "-10:-6,6:10": Default background sample. The sample is given by a set of colon 
    # separated ranges each separated by either whitespace or commas. The string "*" refers to 
    # all points. Note that the background coordinates are relative to the aperture center and not 
    # image pixel coordinates so the endpoints need not be integer. 

    # Default parameter dictionary
    dapall = dict(format='onedspec', nsum=18, t_function='spline3', t_order=5, recenter='no', 
                  b_naverage=-50, resize='no', gain=gain, readnoise=rdnoise, interac='yes',
                  extras='no', trace='no', references="", apertures="1", b_function='legendre')

    # Default values for standard or object
    # For standard stars: background='fit',  weights='variance', clean='yes'
    # For objects:        background='none', weights='none',     clean='no'
    if std:
        dextra = dict(background='fit', weights='variance', clean='yes')
    else:
        dextra = dict(background='none', weights='none', clean='no')
    dapall.update(dextra)

    # Overwrite default values
    dapall.update(kwargs)

    # outfits name
    if outfits is None:
        outfits = fmt_file(fmt_out, fits)

    # Delete output file
    delfiles(outfits)

    # Delete aperture files created from previous extractions
    ap_file = fmt_file(fmt_ap, '.'.join(outfits.split('.')[0:-1]))
    list_ap_files = glob(ap_file)
    delfiles(list_ap_files)

    # Delete aperture database for this file
    root_name_ap = '.'.join(outfits.split('.fit')[0:-1])
    db_files = glob(os.path.join(database, root_name_ap))
    delfiles(db_files)
    
    # When using optimum extraction sometimes the spectrum is zero
    unlearn('apall')
    apall(input=fits, output=outfits, line=line, lower=lower, upper=upper, b_sample=b_sample, **dapall)
    unlearn('apall')

# Sensitivity function -----------------------------------------
def sensf(lstd, sensitivity='sensfunc', standards='standards.dat', observatory='lapalma', graphs='srei',
          caldir="onedstds$iidscal/", extinct="lapalma_extinction.dat", istand='yes', isens='yes', 
          order=5, logfile="", star_names=None, **kwargs):
    # La Palma extinction directory
    # http://www.ing.iac.es/Astronomy/observing/conditions/wlext.html
    # Directory of the standard star file 
    # iraf/noao/lib/onedstds/iidscal/bd253941.dat
    # iraf/noao/lib/onedstds/iidscal/bd332642.dat
    # iraf/iraf/noao/lib/onedstds/oke1990/feige110.dat
    # caldir = "onedstds$oke1990/"
    # iraf/iraf/noao/lib/obsdb.dat -> observatory

    # The ones in IRAF have very poor resolution. Better to use new ones by ESO:
    # https://www.eso.org/sci/observing/tools/standards/spectra/okestandards_rev.html
    # Download them from:
    # ftp://ftp.eso.org/pub/stecf/standards/okestan/

    # The 'f' files list wavelength ( A ), flux ( ergs/cm/cm/s/A * 10**16 ) and flux ( milli-Jy ) and bin (A)
    # The file name consists of a prefix 'f' and the star name
    # The 'm' files list wavelength ( A ), AB magnitude and bin (A). The file 
    # name consists of a prefix 'm' and the star name --> TAKE 'm' FILES FOR IRAF

    # The in the current directory create a file called "standars.men" with the folling format

    # Standard stars in .:
    #
    # mbd25d4655
    # mbd33d2642

    # Set the variable caldir="" and run standard (sensf)

    sensfile = '%s.fits' % sensitivity
    delfiles([sensfile, standards])

    unlearn('standard')
    unlearn('sensfunc')

    lstd = cklist(lstd)
    star_names = cklist(star_names)

    if star_names is None:
        star_names = ['.'.join(os.path.basename(std).split('.')[0:-1]) for std in lstd]

    for std, star_name in zip(lstd, star_names):
        standard(input=std, output=standards, observatory=observatory, caldir=caldir, 
                 extinct=extinct, star_name=star_name, interact=istand)

    # ignorepars='yes' --> All the observations are combined into a single sensitivity function
    sensfunc(standards=standards, sensitivity=sensitivity, extinct=extinct, observatory=observatory,
             logfile=logfile, graphs=graphs, interactive=isens, order=order, ignoreaps='yes', **kwargs)

    unlearn('standard')
    unlearn('sensfunc')

# Remove sky (background) --------------------------------------
def back(fits, outfits=None, fmt='b%s', **kwargs):
    # sample: Lines or columns to be used in the background fits. The default "*" selects all lines or columns
    # sample = '200:400,1000:1200'
    dback = dict(sample=sample, axis=1, low_reject=3., high_reject=3., inter='yes', order=1, niterat=3)
    dback.update(kwargs)

    if outfits is None:
        outfits = fmt_file(fmt, fits)
    delfiles(outfits)

    unlearn('background')
    background(fits, outfits, **dback)
    unlearn('background')

# Apply flux calibration (sensitivity) -------------------------
def calflux(lfits, sensitivity='sensfunc.fits', extinction="lapalma_extinction.dat",
            observatory='lapalma', extinct='yes', outfits=None, fmt_out='cf%s',
            in_files='tmp_in.list', out_files='tmp_out.list'):

    lfits = cklist(lfits)

    if outfits is None:
        outfits = [fmt_file(fmt_out, fits) for fits in lfits]

    delfiles(outfits)

    wfile(in_files, lfits)
    wfile(out_files, outfits)

    unlearn('calibrate')
    # output=pref+'//@'+tlist
    calibrate(input='@'+in_files, output='@'+out_files, sensiti=sensitivity,
              extinct=extinct, extinction=extinction, observatory=observatory, 
              ignoreaps='yes')
    unlearn('calibrate')

    delfiles([in_files, out_files])

# Change dispersion correction (re-bin) ------------------------
def cdispcor(fits, outfits, ref_disp=None, flux='yes', linearize='yes', 
             key_ref='REFSPEC1', add='yes', verify='no', **kwargs):

    # To use the dispersion solution of other image, create a keyword
    # in the header "REFSPEC1" to point to the file with the dispersion
    # solution you want to use
    if ref_disp is not None:
        hedit(fits, key_ref, ref_disp, add=add, verify=verify)

    delfiles(outfits)

    unlearn('dispcor')
    dispcor(fits, outfits, linearize=linearize, flux=flux)
    unlearn('dispcor')

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

##flatccdres(fflat, lflats=tlog[idflat]['file'], suffix='[1]')
#flatccdres(fflat)


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
#    combine(tlog['file'][selection], name, prefix='c')

# Remove cosmic rays if there is only one image per object/standard star ------

#rmcosmic(tlog['file'][idobj], fmt_in='c%s', fmt_out='clean%s')


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

#setair(lobj)


# ---------------------- STANDARD STAR EXTRACTION ------------------------------

lstd = OrderedDict()
lstd['lbd25.fits'] = dict(line=617, lower=-50, upper=50, b_sample='-300:-200,200:300')
lstd['lbd33.fits'] = dict(line=642, lower=-50, upper=50, b_sample='-300:-200,200:300')

#for std in lstd:
#    extract_aperture(std, std=True, **lstd[std])


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
#sensf(lstd, star_names=star_names, caldir="")

## In case you have different sensitivity functions and want to combine them
##imcalc('sens1.fits,sens2.fits', 'sensfunc.fits', '(im1+im2)/2.')

# ------------------------- Remove SKY (background) ------------------------

sample = '200:400,900:1100'
# Output of "back" adds a "b" to the input file by default (use variable "outfits" otherwise)
#back('lgalaxy.fits', sample=sample)


# ---------------------------- Flux calibration ----------------------------

## Output of "calflux" adds a "cf" to the input file by default (use variable "outfits" otherwise)
#calflux('blgalaxy.fits')


# ------------------------- 1D OBJECT EXTRACTION ---------------------------

galfits = 'cfblgalaxy.fits'

lgal = OrderedDict()
# Define your initial appertures (you can change them interactively)
lgal['galaxy_knotb.fits'] = dict(line=3670, lower=-24, upper=21)

#for gal in lgal:
#    extract_aperture(galfits, outfits=gal, find='no', **lgal[gal])


# ---------------------------------------------------------------
# Final unlearn
unlearn_list(ltask)
