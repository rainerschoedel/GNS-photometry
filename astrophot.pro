pro astrophot, field, filter, chip_nr


; PURPOSE: Run photo astrometry on each holographic region
;            
; INPUT:  filter must be string

chip = 'chip' + strn(chip_nr)

; Basic parameters
; size of final image
xaxis = 2400L
yaxis = 1200L  
sub_size_x0 = 600
sub_size_y0 = 600

rebfac = 2
psf_fwhm = 2.0 * rebfac  ; approximate PSF FWHM
;psf_fwhm = 1.5 * rebfac  ; approximate PSF FWHM

maskrad = round(10 * rebfac)         ; circular mask for PSF

satlevel = 4.0e4   ; saturation threshold in long exposure
wt_threshold = 0.1 ; only consider pixels with a
                   ; a weight of at least wt_threshold of maximum
                   ; possible weight in an image
                   ; weight = number of contributing
                   ; exposures (ESO terminology: N)
rebfac = 2   ; rebin factor
sub_size_x0 = sub_size_x0 * rebfac
sub_size_y0 = sub_size_y0 * rebfac
border =  round(maskrad); borders of input images will be masked
                           ; to suppress edge effects (holography,StarFinder)

; Numbers of first and last sub-images (usually not to be changed, see
; runholo.pro)
i_x0 = 0
i_x1 = 6

i_y0 = 0
i_y1 = 2

; number of sub-images from jackknifing
nsub = 10

; input and output directories
basedir = '/home/data/GNS/2015/'

indir = basedir + filter + '/' + field + '/cubeims/chip' + strn(chip_nr) + '/'
photdir = basedir + filter + '/' + field + '/photo/chip' + strn(chip_nr) + '/lists/'
subcubes_dir = basedir + filter + '/' + field + '/subcubes/chip' + strn(chip_nr) + '/'
tmpdir = basedir + filter + '/' + field + '/tmp/tmp' + strn(chip_nr) + '/'

for i_x = i_x0, i_x1 do begin
  for i_y = i_y0, i_y1 do begin 

    name = strn(i_x) + '_' +  strn(i_y)

   ; Compute longexposure image to identify saturated sources
   ; ---------------------------------------------------------
   readcol, subcubes_dir + 'masklist_' + name + '.txt', mnames, FORMAT='A'
   readcol, subcubes_dir + 'list_' + name + '.txt', cnames, FORMAT='A'
   support = fltarr(sub_size_x0, sub_size_y0)
   lxp = fltarr(sub_size_x0/2, sub_size_y0/2)
   nm = n_elements(mnames)
   for j = 0, nm-1 do begin
    maskcube = readfits(subcubes_dir + mnames[j])
    support = support + total(maskcube,3)
    cube = readfits(subcubes_dir + cnames[j])
    lxp = lxp + total(cube,3)
   endfor
   accept = where(support gt 0)
   lxp[accept] = lxp[accept]/support[accept]
   ; rebin lxp
   lxp = CREBIN(lxp,sub_size_x0,sub_size_y0)
   writefits, tmpdir + 'lxp_' + name +'.fits', lxp
   writefits, tmpdir + 'support_' + name +'.fits', support

   ; Mask saturated pixels generously
   ; (so that no reference star near a saturated source is chosen)
    sat_pixels = where(lxp gt satlevel,n_saturated)
    sat_mask = replicate(1,sub_size_x0,sub_size_y0)
    if (n_saturated gt 0) then begin
      xy = array_indices(sat_mask,sat_pixels)
      xbad = xy[0,*]
      ybad = xy[1,*]
      for i_sat = 0, n_saturated-1 do begin
        sat_mask = circ_mask(sat_mask,xbad[i_sat],ybad[i_sat],maskrad,/INNER)
      endfor
    endif
    writefits, tmpdir + 'sat_mask_' + name +'.fits', sat_mask

   ; Read holography image, noise and weight map
   im = readfits(indir + 'holo_' +  name + '.fits.gz')
   noise = readfits(indir + 'holo_' + name + '_sigma.fits.gz')
   exp = readfits(indir + 'holo_' + name + '_expmap.fits.gz')

   ; Create mask (borders, low S/N regions
   mask = fltarr(sub_size_x0,sub_size_y0)
   exp = exp/max(exp)
   good = where(exp gt wt_threshold,complement=bad)
   mask[bad] = 0
   mask[good] = 1
   mask = mask * sat_mask
   mask[sub_size_x0-border-1:sub_size_x0-1,*] = 0
   mask[*,sub_size_y0-border-1:sub_size_y0-1] = 0
   mask[*,0:border-1] = 0
   mask[0:border-1,*] = 0
   writefits, tmpdir + 'mask_' + name +'.fits', mask

   ; EXTRACT PSF
   ; important: Pass mask to extractpsf.pro, but do not
   ; apply it to the image/noise map
   psf = extractpsf(im, noise, mask, maskrad, psf_fwhm, tmpdir, /UNWEIGHTED)
   writefits, indir + 'psf_'+ name + '.fits' , psf

   ; RUN STARFINDER
   ; write results to file
   estim_bg = 1
   back_box = maskrad
   min_correlation = 0.7
   correl_mag = 4
   deblend = 0
   deblost = 0
   niter = 2
   rel_thresh = 1
   GUIDE_X = ''
   GUIDE_Y = ''
   Threshold = [5.,5.]
   starfinder, im, psf, X_BAD=x_bad, Y_BAD = y_bad, $
        BACKGROUND = background, BACK_BOX = back_box, $
        threshold, REL_THRESHOLD = rel_thresh, /PRE_SMOOTH, $
        NOISE_STD = noise, min_correlation, $
        CORREL_MAG = correl_mag, INTERP_TYPE = 'I', $
        ESTIMATE_BG = estim_bg, DEBLEND = deblend, N_ITER = niter, SILENT=1, $
        GUIDE_X = guide_x, GUIDE_Y = guide_y, $
        SV_SIGMA_R = 0.0085, SV_SIGMA_A= 0.0050, $
        x, y, f, sx, sy, sf, c, STARS = stars, $
        LOGFILE = logfilename, /CUBIC;, /NO_SLANT
    writefits, indir + 'subtracted_' + name + '.fits', (im - stars), /COMPRESS
;    writefits, indir + 'bg_' + name + '.fits', background, /COMPRESS
    forprint, TEXTOUT =  photdir + 'stars_' + name + '.txt', x, y, f, sx, sy, sf, c, format='(7(A, 3X))',/NOCOMMENT

  ; RUN STARFINDER on jackknife images
    for i_s = 1, nsub do begin
;      im = readfits(indir + 'mosaic_' + name + '_s' + strn(i_s) + '.fits.gz')
;      noise = readfits(indir + 'mosaic_' + name + '_sigma_s' + strn(i_s) + '.fits.gz')
      im = readfits(indir + 'holo_' + name + '_s' + strn(i_s) + '.fits.gz')
;      noise = readfits(indir + 'holo_' + name + '_sigma_s' + strn(i_s) + '.fits.gz')
      starfinder, im, psf, X_BAD=x_bad, Y_BAD = y_bad, $
        BACKGROUND = background, BACK_BOX = back_box, $
        threshold, REL_THRESHOLD = rel_thresh, /PRE_SMOOTH, $
        NOISE_STD = noise, min_correlation, $
        CORREL_MAG = correl_mag, INTERP_TYPE = 'I', $
        ESTIMATE_BG = estim_bg, DEBLEND = deblend, N_ITER = niter, SILENT=1, $
        GUIDE_X = guide_x, GUIDE_Y = guide_y, $
        SV_SIGMA_R = 0.0085, SV_SIGMA_A= 0.0050, $
        x, y, f, sx, sy, sf, c, STARS = stars, $
        LOGFILE = logfilename, /CUBIC;, /NO_SLANT
      forprint, TEXTOUT =  photdir + 'stars_' + name + '_s' + strn(i_s) + '.txt', x, y, f, sx, sy, sf, c, format='(7(A, 3X))',/NOCOMMENT
   endfor ; loop over jacknife images
     
 endfor ; loop over i_x
endfor  ; loop over i_y

end
