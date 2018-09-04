# This is a tutorial example for the introductory course "Long-slit reduction with PyRAF"
# given at the "School on Long-Slit Spectroscopy" at University of Rwanda from August 26th 
# to September 1st of 2018:
#
# https://sites.google.com/site/schoolonspectroscopy/home
#
# Author: Ruben Garcia-Benito (RGB)
#
# MODULE: relspyraf (Reduce Long Slit PyRAF)

# %%%%%%%%%%%%%%%% IMPORTS %%%%%%%%%%%%%%%%%%%
# PYRAF import --------------------
from pyraf import iraf
iraf.noao.imred(_doprint=0)
iraf.noao.twodspec(_doprint=0)
iraf.noao.onedspec(_doprint=0)
iraf.noao.imred.ccdred(_doprint=0)
iraf.noao.twodspec.longslit(_doprint=0)
iraf.noao.twodspec.apextract(_doprint=0)
from iraf import imcombine, background, apall
from iraf import ccdproc, zerocombine, imstatistics
from iraf import flatcombine, response, hedit, unlearn
from iraf import identify, reidentify, fitcoords, transform
from iraf import standard, sensfunc, calibrate, setairmass
# -------------------------------

from astropy.io import fits as pyfits
from glob import glob
import astropy.table
import numpy as np
import astropy
import os
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ---------------------- FUNCTIONS -----------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Check if object is list --------------------------------------
def cklist(obj):
    if not isinstance(obj, (list, tuple, np.ndarray, astropy.table.table.Column)) and obj is not None:
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

# Change name of a list/array of files with a given format -----
def fmt_files(names, fmt='%s', prefix=None, suffix=None, dirname=None):
    if prefix is not None:
        fmt = prefix + fmt
    if suffix is not None:
        fmt = fmt + suffix
    if dirname is not None:
        fmt = os.path.join(dirname, os.path.basename(fmt))
    names = cklist(names)
    new_names = [fmt % name for name in names]
    return new_names

# Unlearn list of functions ------------------------------------
def unlearn_list(ltask):
    ltask = cklist(ltask)
    for task in ltask: 
        unlearn(task)

# Extract HDU (chip) from a multi HDU fits ---------------------
def extract_hdu(input_fits, outfits, extension, primary_header=True, exclude_primary='COMMENT', exclude=None, dkeys=None):
    lhdu = pyfits.open(input_fits)
    exclude = cklist(exclude)
    exclude_primary = cklist(exclude_primary)
    if (extension + 1) > len(lhdu):
        print ('>>> Extension "%i" is larger than the number of HDUs [%i]' % (extension, len(lhdu)))
        lhdu.close()
        return 
    hdu = lhdu[extension]
    nhdu = pyfits.PrimaryHDU(data=hdu.data, header=hdu.header)
    if exclude is not None:
        for key in exclude:
            if key in nhdu.header:
                nhdu.header.remove(key)
    if primary_header:
        phdr = lhdu[0].header 
        for key in phdr:
            if (not key in nhdu.header) and (not key in exclude_primary):
                nhdu.header.set(key, value=phdr[key], comment=phdr.comments[key])
    lhdu.close()
    if isinstance(dkeys, dict):
        for key in dkeys:
            nhdu.header[key] = dkeys[key]
    nhdu.writeto(outfits, clobber=True)

# Master flat --------------------------------------------------
def flatccdres(flat_files, flat='flat.fits', nflat='nflat.fits', cflat='cflat.fits', iccd='no', iresponse='yes', 
	combine='average', reject='avsigclip', ccdtype="", process='no', fixpix='no', overscan='yes', 
	zerocor='yes', trim='yes', darkcor='no', flatcor='no', illumco='no', fringec='no', function='legendre', 
	corder=5, rorder=18, readaxi='column', biassec=None, zero=None, trimsec=None, low_reject=3,
	high_reject=3, scale='mode', rdnoise=0.0, gain=1.0, lflats=None, prefix=None, suffix=None):

    unlearn('flatcombine')
    unlearn('ccdproc')
    unlearn('response')

    delfiles(flat)
    if lflats is not None:
        wfile(flat_files, lflats, prefix=prefix, suffix=suffix)

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

    unlearn('flatcombine')
    unlearn('ccdproc')
    unlearn('response')

# ccdproc function  --------------------------------------------
def ccdp(files, prefix=None, suffix=None, input_files='tmp_input.list', output_files='tmp_output.list', fmt_out='c%s', **kwargs):

    # Dictionary with default parameters for ccdproc
    dpar = dict(ccdtype="", fixpix='no', overscan='yes', trim='yes', zerocor='yes', flatcor='no', darkcor='no',
	        illumco='no', fringec='no', zero="", illum="", function='legendre', order=4, 
                interactive='no', readaxi='column', trimsec=None, biassec=None, flat="")

    # Update parameters dictionary for ccdproc if new are provided (kwargs)
    dpar.update(kwargs)

    # Make "files" a list, in case there is only one file
    files = cklist(files)

    # Write the list of input files to a temporary file
    wfile(input_files, files, prefix=prefix, suffix=suffix)

    # Create a list of names with the output name after ccdproc is applied
    new_files_ccdproc = [fmt_out % item for item in files]
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
def transf(lfits, fitnames, interp='linear', database='database', flux='yes', y1='INDEF', prefix=None, suffix=None,
           input_files='tmp_input.list', output_files='tmp_output.list', fmt_out='l%s', dispaxis=None, **kwargs):

    fitnames = '.'.join(fitnames.split('.')[0:-1])

    if isinstance(lfits, str) and input_files is None:
        input_files = lfits
        lfits = np.loadtxt(input_files, dtype=np.str, unpack=True)
    else:
        lfits = cklist(lfits)
        wfile(input_files, lfits, prefix=prefix, suffix=suffix)

    wfits = [fmt_out % fits.split('[')[0] for fits in lfits]
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
            expname='EXPTIME', rdnoise=0.0, gain=1.0, suffix=None, prefix=None, 
            input_files='tmp_input.list', delete_list=True):

    wfile(input_files, lfits, suffix=suffix, prefix=prefix)
    delfiles(outfits)

    unlearn('imcombine')
    imcombine('@'+input_files, output=outfits, combine=combine, reject=rej, scale=scale,
		weight=weight, expname=expname, rdnoise=rdnoise, gain=gain, nlow=nlow)
    unlearn('imcombine')

    if delete_list:
        delfiles(input_files)

# Remove cosmic rays of single images --------------------------
def rmcosmic(lfits, xorder=5, yorder=5, sigfrac=2.0, niter=5, objlim=2, sigclip=3.0, 
             fmt_in=None, fmt_out='cr%s', rdnoise=0.0, gain=1.0):

    # Download "lacos_spec.cl" from: http://www.astro.yale.edu/dokkum/lacosmic/lacos_spec.cl
    # Save it into a directory ("lacos_directory" for our purposes)
    # In your IRAF "login.cl" include the following line:
    # task    lacos_spec = lacos_directory/lacos_spec.cl

    from iraf import lacos_spec

    lfits = cklist(lfits)
    fmt_in = '%s' if fmt_in is None else fmt_in

    unlearn('lacos_spec')
    for fits in lfits:
        fits   = fmt_in  % fits
        output = fmt_out % fits
        mask   = 'mask%s' % fits
        delfiles([mask, output])

        lacos_spec(fits, output, mask, gain=gain, readn=rdnoise, xorder=xorder, yorder=yorder,
                   sigclip=sigclip, sigfrac=sigfrac, objlim=objlim, niter=niter)

        delfiles(mask)
    unlearn('lacos_spec')

# Set airmass in images ----------------------------------------
def setair(files, observatory='lapalma', equinox='EPOCH', ut='UTSTART', input_files='tmp_input.list', prefix=None, suffix=None):

    files = cklist(files)
    wfile(input_files, files, prefix=prefix, suffix=suffix)

    unlearn('setairmass')
    setairmass(images='@'+input_files, observatory=observatory, equinox=equinox, ut=ut)
    unlearn('setairmass')

    delfiles(input_files)


# Extract aperture ---------------------------------------------
def extract_aperture(fits, line, lower=-5, upper=5, b_sample='-10:-6,6:10',
    outfits=None, fmt_out='ap%s', fmt_ap='%s.000?.fits', database='database/', 
    rdnoise=0.0, gain=1.0, std=False, **kwargs):

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
          caldir='onedstds$iidscal/', extinction='lapalma_extinction.dat', istand='yes', isens='yes', 
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
                 extinction=extinction, star_name=star_name, interact=istand)

    # ignorepars='yes' --> All the observations are combined into a single sensitivity function
    sensfunc(standards=standards, sensitivity=sensitivity, extinction=extinction, observatory=observatory,
             logfile=logfile, graphs=graphs, interactive=isens, order=order, ignoreaps='yes', **kwargs)

    unlearn('standard')
    unlearn('sensfunc')

# Remove sky (background) --------------------------------------
def back(fits, outfits=None, fmt='b%s', **kwargs):
    # sample: Lines or columns to be used in the background fits. The default "*" selects all lines or columns
    # sample = '200:400,1000:1200'
    dback = dict(axis=1, low_reject=3., high_reject=3., inter='yes', order=1, niterat=3)
    dback.update(kwargs)

    if outfits is None:
        outfits = fmt_file(fmt, fits)
    delfiles(outfits)

    unlearn('background')
    background(fits, outfits, **dback)
    unlearn('background')

# Apply flux calibration (sensitivity) -------------------------
def calflux(lfits, sensitivity='sensfunc.fits', extinction="lapalma_extinction.dat", observatory='lapalma', 
            extinct='yes', outfits=None, fmt_out='cf%s', prefix=None, suffix=None,
            in_files='tmp_in.list', out_files='tmp_out.list'):

    lfits = cklist(lfits)

    if outfits is None:
        outfits = [fmt_file(fmt_out, fits) for fits in lfits]

    delfiles(outfits)

    wfile(in_files, lfits, prefix=prefix, suffix=suffix)
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
