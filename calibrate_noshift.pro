pro calibrate_noshift, field, band, chip


;t_exp = 0.85
t_exp = 1.26
dmax_sirius = 3.  ; max HAWK-I pixel distance between SIRIUS and HAWK-I for calibration star selection
                 ; 1 pixel = 0.053" 
dmax = 2.        ; max HAWK-I pixel distance for identification of reference stars
mag_lim_VVV = 18 ; only use stars brighter than this (J-band) from VVV

;dmax_near_stars = 10 ; eliminate nearby stars within this amount of HAWK-I pixels (10 pix = 0.53")
; Any secondary star within within delta_r of a photometric reference
; star needs to be at least delta_mag magnitudes fainter
delta_r = 10. ; ; with 0.053"  pixel scale 20 = 1" 
delta_mag = 2.5
SIGMA_CUT = 2.0 ; sigma cut to be applied in resistant_mean

magerr_si = 0.03 ; max acceptable SIRIUS magnitude uncertainty for photometric calibration stars
magerr_hawki = 0.03 ; max acceptable HAWK-I magnitude uncertainty for photometric calibration stars

mag_max = 12.5 ; Bright magnitude cut off for calibration stars (HAWK-I)
scale_holo = 0.053 ; pixel scale of holography image (" /pixel)
scale_vvv = 0.34   ; pixel scale of VVV image (" /pixel)

; Radius for creating median ZP
rmed = 1200. ; with 0.053"  pixel scale 1200 ~1 arcmin
minexp_cal = 0.3 ; valid calibration stars should be covered by at leasy this fraction of all exposures

; Define input and output paths
; and file names
; -----------------------------

indir = '/home/data/GNS/2015/' + band + '/' + strn(field) + '/cubeims/chip' + strn(chip) + '/'
indir_list = '/home/data/GNS/2015/' + band + '/' + strn(field) + '/photo/chip' + strn(chip) + '/lists/'

; output
res_dir = '/home/data/GNS/2015/' + band + '/' + strn(field) + '/photo/chip' + strn(chip) + '/lists/'
res_dir_im = '/home/data/GNS/2015/' + band + '/' + strn(field) + '/cubeims/chip' + strn(chip) + '/'

; for Ks and J bands use the images and lists aligned with the H-band
; (produced by align_fields.pro)
if (band eq 'H') then align = '' else align = '_aligned_H'
star_list = 'astrophot' + align + '_noshift.txt'

; The J-band data will be used to establish astrometry (via VVV)
indir_list_J = '/home/data/GNS/2015/H' + '/' + strn(field) + '/photo/chip' + strn(chip) + '/lists/' 
star_list_J = 'astrophot_aligned_H.txt' ; MUST BE '*_aligned_H*'!
res_dir_imJ = '/home/data/GNS/2015/H' + strn(field) + '/photo/chip' + strn(chip) + '/'


tmpdir =  '/home/data/GNS/2015/' + band + '/' + strn(field) + '/tmp/' 


; Calibration files and images are taken from VVV J-band
vvv_dir = '/home/data/VVV/Fields/J/'
ref = readfits(vvv_dir + 'Field' + field  + '.fits.gz',vvvheader)
ref_file = vvv_dir + 'Field' + field  + '_stars.txt'
;crval1 = SXPAR(vvvheader, 'CRVAL1')
;crval2 = SXPAR(vvvheader, 'CRVAL2')
;crval1 = double(265.965030272004)
;crval2 = double(-29.2927997156051)
;print, crval1, FORMAT='(D)'
;print, crval2, FORMAT='(D)'
;STOP

; Holography image is aligned just as lnx_aligned (see alignquadrants)
; To transfrom HAWK-I positions from holography image
; to VVV we merely need to scale and shift (x_off, y_off from alignquadrants)
; VVV --> HAWKI
rebfac = 2
scale = (scale_vvv/scale_holo)
RESTORE, '../scripts/xy_off_f' + strn(field) + '.SAV'


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; (1) Calibrate astrometry
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; Read images

im = readfits(indir +  strn(field) + '_chip' + strn(chip) + '_holo' + align + '.fits.gz', im_head)
noise = readfits(indir +  strn(field) + '_chip' + strn(chip) + '_noise' + align + '.fits.gz', noise_head)
exp = readfits(indir +  strn(field) + '_chip' + strn(chip) + '_exp' + align + '.fits.gz', exp_head)

sz = size(im)
sz_vvv = size(ref)
xsize_vvv = sz_vvv[1]
ysize_vvv = sz_vvv[2]
xsize_hawki = sz[1]
ysize_hawki = sz[2]
sub_size = 1200 ; size of sub-fields for check of ZP variations across the field

; Read star lists
readcol, ref_file, x_vvv, y_vvv, a_vvv, d_vvv, J_vvv, Ks_vvv, FORMAT='D,D,D,D,F,F'  ; astrometric reference file from VVV
vvv_good = where(J_vvv lt mag_lim_VVV)
x_vvv = x_vvv[vvv_good]
y_vvv = y_vvv[vvv_good]
a_vvv = a_vvv[vvv_good]
d_vvv = d_vvv[vvv_good]
J_vvv = J_vvv[vvv_good]
x_vvv_scaled = x_vvv*scale - rebfac*x_off[chip-1]
y_vvv_scaled = y_vvv*scale - rebfac*y_off[chip-1]
readcol, indir_list + star_list, x, y, f, sx, sy, sf, corr, ndet, FORMAT='D,D,D,D,D,D,F,F' ; HAWK-I list for field to be calibrated


; find common stars VVV J and HAWK-I
; create model image to check whether scaled VVV positions
; correspond to HAWK-I image
; Create map of all potential reference stars for this field
dat = ptr_new({X_size: 10, Y_size: 10, Sigma_x: 1.5, Sigma_y: 1.5, Angle: 0.0})
f_vvv = 1.0e7 * 10^(-0.4*J_vvv)
imvvv  = image_model(x_vvv_scaled,y_vvv_scaled,f_vvv,xsize_hawki,ysize_hawki,'gaussian', dat)
writefits, tmpdir +'HAWKI_in_VVV.fits', imvvv, /COMPRESS

compare_lists, x, y, x_vvv_scaled, y_vvv_scaled, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax,$
   SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2   
n_stars = n_elements(subc1)
print, 'Number of common stars: ' + strn(n_stars)
 
;STOP

; ;Now compute transformation VVV --> HAWK-I to improve identification
; of common stars
n_astro = n_elements(subc1)
print, 'Using ' + strn(n_astro) +  ' stars to compute preliminary alignment.'
polywarp, x[subc1], y[subc1], x_vvv[subc2], y_vvv[subc2], 1, Kx, Ky
print, Kx
print, Ky
xi = Kx[0,0] + Kx[0,1]*x_vvv + Kx[1,0]*y_vvv + Kx[1,1]*x_vvv*y_vvv
yi = Ky[0,0] + Ky[0,1]*x_vvv + Ky[1,0]*y_vvv + Ky[1,1]*x_vvv*y_vvv
compare_lists, x, y, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
nc = n_elements(subc1)
print, 'Iteration 0'
print, 'Found ' + strn(nc) + ' common stars.'
; Use robust statistics to trim outliers
dp = sqrt((x1c-x2c)^2 + (y1c-y2c)^2)
RESISTANT_Mean,dp,3.0,Mean_h,Sigma_h,Num_Rej,goodvec = goodvec
subc1 = subc1[goodvec]
subc2 = subc2[goodvec]
nc = n_elements(subc1)
print, 'After outlier trimming we have ' + strn(nc) + ' common stars.'
;STOP

; iterative degree 1 alignment
; ------------------------------

  lim_it = 1
  count=0
  comm=[]
  it=0
  lim_it=1 ;consider convergence when the number of common stars  'lim_it' times.
	 
;  while count lt lim_it do begin
;  it=it+1
for it = 1, 1 do begin
  degree = 1
  polywarp, x[subc1], y[subc1], x_vvv[subc2], y_vvv[subc2], 1, Kx, Ky
  print, Kx
  print, Ky
  xi = Kx[0,0] + Kx[0,1]*x_vvv + Kx[1,0]*y_vvv + Kx[1,1]*x_vvv*y_vvv
  yi = Ky[0,0] + Ky[0,1]*x_vvv + Ky[1,0]*y_vvv + Ky[1,1]*x_vvv*y_vvv
  compare_lists, x, y, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
  nc = n_elements(subc1)
  print, 'Iteration ' + strn(it)
  print, 'Found ' + strn(nc) + ' common stars.'
  ; Use robust statistics to trim outliers
  dp = sqrt((x1c-x2c)^2 + (y1c-y2c)^2)
  RESISTANT_Mean,dp,3.0,Mean_h,Sigma_h,Num_Rej,goodvec = goodvec
  subc1 = subc1[goodvec]
  subc2 = subc2[goodvec]
  nc = n_elements(subc1)
  print, 'After outlier trimming we have ' + strn(nc) + ' common stars.'
  comm=[comm,nc]
    if (n_elements(comm) gt 2) then begin;
	   if comm[-2] ge comm[-1] then begin
	   count=count+1
	   endif else begin
	   count=0
	   endelse
	endif
;  endwhile
endfor
compare_lists, x, y, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2

; Create region of VVV alignment stars on HAWK-I image
n_stars = n_elements(subc2)
openw, out, 'VVV.reg', /get_lun
printf, out, 'global color=magenta dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
printf, out, 'image'
for s = 0, n_stars-1 do begin
   printf, out, 'circle('+ strn(x2c[s])+','+strn(y2c[s])+',5)'
endfor
free_lun, out
; Create region of HAWKI alignment stars on HAWK-I image
openw, out, 'HAWKI.reg', /get_lun
printf, out, 'global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
printf, out, 'image'
for s = 0, n_stars-1 do begin
   printf, out, 'circle('+ strn(x[subc1[s]])+','+strn(y[subc1[s]])+',5)'
endfor
free_lun, out

; Plot distributions of differences of
; astrometric  positions between reference stars in HAWK-I and VVV
; -----------------------------------------------------------------

dx_tmp = (xi[subc2]-x[subc1])*scale_holo
dy_tmp = (yi[subc2]-y[subc1])*scale_holo

set_plot,'PS',/interpolate
device, XOFFSET=0, YOFFSET=0, $
      FILENAME =  res_dir + 'position_uncertainty_pixels_' + strn(chip) + '_noshift.eps', XSIZE=20., YSIZE=12., $
      /portrait, /color,BITS_PER_PIXEL=8, encapsulated=1
!P.MULTI=[0,2,1]
!P.CHARSIZE=1.2
!P.THICK=4.0
!P.CHARTHICK=4.0

binsize = 0.01
cgHistoplot, dx_tmp, /FILL, HISTDATA=h1, LOCATIONS=loc, BINSIZE=binsize, XTITLE = 'Distance in X-axis (")', YTITLE = 'Number of Stars',xrange = [-0.25,0.25], xticks = 4

binCenters1 = loc + (binsize / 2.0)
yfit = GaussFit(binCenters1, h1, coeff, NTERMS=3, sigma = error)
cgPlot, binCenters1, yfit, COLOR='dodger blue', THICK=2, /OVERPLOT
  
centerfit = String(coeff[1], FORMAT='(F0.3)')
cgText, 0.17, 0.84, /NORMAL, 'Center ' + centerfit + ' arcsec', COLOR='navy'
cgText, 0.17, 0.80, /NORMAL, 'Error ' + strn(coeff[2], FORMAT='(F0.3)') + ' arcsec', COLOR='navy'

cgHistoplot, dy_tmp, /FILL, HISTDATA=h1, LOCATIONS=loc, BINSIZE=binsize, XTITLE = 'Distance in Y-axis (")' , YTITLE = '' ,xrange = [-0.25,0.25], xticks = 4
binCenters1 = loc + (binsize / 2.0)

yfit = GaussFit(binCenters1, h1, coeff, NTERMS=3, sigma = error)
cgPlot, binCenters1, yfit, COLOR='dodger blue', THICK=2, /OVERPLOT
  
centerfit = String(coeff[1], FORMAT='(F0.3)')
cgText, 0.65, 0.84, /NORMAL, 'Center ' + centerfit + ' arcsec', COLOR='navy'
cgText, 0.65, 0.80, /NORMAL, 'Error ' + strn(coeff[2], FORMAT='(F0.3)') + ' arcsec', COLOR='navy'

device, /close

 
; Now compute and apply the astrometric solution 
; ----------------------------------------------

astro = SOLVE_ASTRO(a_vvv[subc2], d_vvv[subc2], x[subc1]+1, y[subc1]+1, success=success, distort = 'none');, reject = 2.0, NITER=3)
print, astro


putast, im_head, astro          ;, crpix, crval, ctype
putast, noise_head, astro;, crpix, crval, ctype
putast, exp_head, astro;, crpix, crval, ctype
sxaddpar, im_head, 'Band',band,'HAWK-I VLT'
sxaddpar, noise_head, 'Band',band,'HAWK-I VLT'
sxaddpar, exp_head, 'Band',band,'HAWK-I VLT'

;hhh = sxpar(im_head,'CRVAL1')
;print, hhh
;print, hhh, FORMAT='(D)'
;STOP

; Write astrometrically calibrated images to disc
;writefits, res_dir_im + strn(field) + '_chip' + strn(chip) + '_holo_cal.fits', im, im_head, /COMPRESS
;writefits, res_dir_im + strn(field) + '_chip' + strn(chip) + '_noise_cal.fits', noise, noise_head, /COMPRESS
;writefits, res_dir_im + strn(field) + '_chip' + strn(chip) + '_exp_cal.fits', exp, exp_head, /COMPRESS



; Write astrometrically calibrated list to disc
; combine statistical and PSF uncertainties
; NEED TO READ INPUT LIST AGAIN
; ================================================
XY2AD, x, y, astro, a, d 
star_list = 'calibrated_' + band + '_' + strn(field) + '_' + strn(chip) + '_noshift.txt'
forprint, TEXTOUT =  indir_list + star_list, a, d, x, y, f, sx, sy, sf, format='(8F16.6)',/NOCOMMENT

; Plot distributions of differences of CALIBRATED
; astrometric  positions between reference stars in HAWK-I and VVV
; in RA, Dec and in GALACTIC coordinates
; -----------------------------------------------------------------



; First in equatorial coordinates
; ----------------------------

dra_tmp =  cos(d_vvv[subc2] * !PI/180.) * (a_vvv[subc2] - a[subc1]) * 3600.
ddec_tmp = (d_vvv[subc2] - d[subc1]) * 3600.

set_plot,'PS',/interpolate
device, XOFFSET=0, YOFFSET=0, $
      FILENAME =  res_dir + 'dpos_equatorial_' + strn(chip) + '_noshift.eps', XSIZE=20., YSIZE=12., $
      /portrait, /color,BITS_PER_PIXEL=8, encapsulated=1
!P.MULTI=[0,2,1]
!P.CHARSIZE=1.2
!P.THICK=4.0
!P.CHARTHICK=4.0

binsize = 0.01
cgHistoplot, dra_tmp, /FILL, HISTDATA=h1, LOCATIONS=loc, BINSIZE=binsize, XTITLE = 'Difference in R.A. (")', YTITLE = 'Number of Stars',xrange = [-0.25,0.25], xticks = 4

binCenters1 = loc + (binsize / 2.0)
yfit = GaussFit(binCenters1, h1, coeff, NTERMS=3, sigma = error)
cgPlot, binCenters1, yfit, COLOR='dodger blue', THICK=2, /OVERPLOT
  
centerfit = String(coeff[1], FORMAT='(F0.3)')
cgText, 0.17, 0.84, /NORMAL, 'Center ' + centerfit + ' arcsec', COLOR='navy'
cgText, 0.17, 0.80, /NORMAL, 'Error ' + strn(coeff[2], FORMAT='(F0.3)') + ' arcsec', COLOR='navy'

cgHistoplot, ddec_tmp, /FILL, HISTDATA=h1, LOCATIONS=loc, BINSIZE=binsize, XTITLE = 'Difference in Dec. (")' , YTITLE = '' ,xrange = [-0.25,0.25], xticks = 4
binCenters1 = loc + (binsize / 2.0)

yfit = GaussFit(binCenters1, h1, coeff, NTERMS=3, sigma = error)
cgPlot, binCenters1, yfit, COLOR='dodger blue', THICK=2, /OVERPLOT
  
centerfit = String(coeff[1], FORMAT='(F0.3)')
cgText, 0.65, 0.84, /NORMAL, 'Center ' + centerfit + ' arcsec', COLOR='navy'
cgText, 0.65, 0.80, /NORMAL, 'Error ' + strn(coeff[2], FORMAT='(F0.3)') + ' arcsec', COLOR='navy'

device, /close

; Now in galactic coordinates
; ----------------------------

GLACTC, a_vvv, d_vvv, 2000.0,  glvvv, gbvvv, 1, /DEGREE
GLACTC, a, d, 2000.0, gl, gb, 1, /DEGREE

dra_tmp =  cos(gbvvv[subc2] * !PI/180.) * (glvvv[subc2] - gl[subc1]) * 3600.
ddec_tmp = (gbvvv[subc2] - gb[subc1]) * 3600.

set_plot,'PS',/interpolate
device, XOFFSET=0, YOFFSET=0, $
      FILENAME =  res_dir + 'dpos_galactic_' + strn(chip) + '_noshift.eps', XSIZE=20., YSIZE=12., $
      /portrait, /color,BITS_PER_PIXEL=8, encapsulated=1
!P.MULTI=[0,2,1]
!P.CHARSIZE=1.2
!P.THICK=4.0
!P.CHARTHICK=4.0

binsize = 0.01
cgHistoplot, dra_tmp, /FILL, HISTDATA=h1, LOCATIONS=loc, BINSIZE=binsize, XTITLE = 'Difference in gl (")', YTITLE = 'Number of Stars',xrange = [-0.25,0.25], xticks = 4

binCenters1 = loc + (binsize / 2.0)
yfit = GaussFit(binCenters1, h1, coeff, NTERMS=3, sigma = error)
cgPlot, binCenters1, yfit, COLOR='dodger blue', THICK=2, /OVERPLOT
  
centerfit = String(coeff[1], FORMAT='(F0.3)')
cgText, 0.17, 0.84, /NORMAL, 'Center ' + centerfit + ' arcsec', COLOR='navy'
cgText, 0.17, 0.80, /NORMAL, 'Error ' + strn(coeff[2], FORMAT='(F0.3)') + ' arcsec', COLOR='navy'

cgHistoplot, ddec_tmp, /FILL, HISTDATA=h1, LOCATIONS=loc, BINSIZE=binsize, XTITLE = 'Difference in gb (")' , YTITLE = '' ,xrange = [-0.25,0.25], xticks = 4
binCenters1 = loc + (binsize / 2.0)

yfit = GaussFit(binCenters1, h1, coeff, NTERMS=3, sigma = error)
cgPlot, binCenters1, yfit, COLOR='dodger blue', THICK=2, /OVERPLOT
  
centerfit = String(coeff[1], FORMAT='(F0.3)')
cgText, 0.65, 0.84, /NORMAL, 'Center ' + centerfit + ' arcsec', COLOR='navy'
cgText, 0.65, 0.80, /NORMAL, 'Error ' + strn(coeff[2], FORMAT='(F0.3)') + ' arcsec', COLOR='navy'

device, /close


end
