pro align_fields, field, chip_nr


; PURPOSE: Take the input lists from merge_sublists.pro and
; merge_subims.pro
;      
; Align J+Ks images and lists with the H-band 
;

chip = 'chip' + strn(chip_nr)
max_distance = 1 ; Correponds to dmax in search for 
              ; common stars between a subimages
              ;  for fine-alignment of subimages
              ; with rebfac = 2, max_distance = 2
              ; corresponds to half the resolution limit
 
;
; input and output directories
basedir = '/home/data/GNS/2015/'

photdirh = basedir + 'H/' + field + '/photo/chip' + strn(chip_nr) + '/lists/'
imdirh = basedir + 'H/' + field + '/cubeims/chip' + strn(chip_nr) + '/'
photdirj = basedir + 'J/' + field + '/photo/chip' + strn(chip_nr) + '/lists/'
imdirj = basedir + 'J/' + field + '/cubeims/chip' + strn(chip_nr) + '/'
photdirk = basedir + 'Ks/' + field + '/photo/chip' + strn(chip_nr) + '/lists/'
imdirk = basedir + 'Ks/' + field + '/cubeims/chip' + strn(chip_nr) + '/'


name = field + '_' + chip

readcol, photdirh + 'astrophot.txt', xh, yh, fh, sxh, syh, sfh, ch, ndet

; Determine offsets for J
; Apply offsets to lists and images
; and write them to file
; --------------------------
photdir = photdirj
imdir = imdirj
readcol, photdir + 'astrophot.txt', x, y, f, sx, sy, sf, c, ndet
compare_lists, xh, yh, x, y, x1c, y1c, x2c, y2c, $
           MAX_DISTANCE = max_distance, $
           SUBSCRIPTS_1 = subc1, SUBSCRIPTS_2 = subc2
RESISTANT_Mean, (x1c-x2c), 3.0, x_off, x_off_Sigma, Num_Rej, goodvec = goodvec
RESISTANT_Mean, (y1c-y2c), 3.0, y_off, y_off_Sigma, Num_Rej, goodvec = goodvec
x = x + x_off
y = y + y_off
print, 'Offsets in x and y of J relative to H: ' + strn(x_off) + ' +- ' + strn(x_off_Sigma) +  ', ' + strn(y_off) + ' +- ' + strn(y_off_Sigma)
forprint, TEXTOUT=photdir + 'astrophot_aligned_H.txt', x, y, f, sx, sy, sf, c, ndet, FORMAT='(8F13.3)', $
              COMMENT='#x[pix]  y[pix]  f[ADU]  dx[pix]  dy[pix]  c0[corr]  output from align_fields.pro'

im = readfits(imdir + name + '_holo.fits.gz')
noise = readfits(imdir + name + '_noise.fits.gz')
wt = readfits(imdir + name + '_exp.fits.gz')
writefits, imdir + name + '_holo_aligned_H.fits', image_shift(im,x_off,y_off), /COMPRESS
writefits, imdir + name + '_noise_aligned_H.fits', image_shift(noise,x_off,y_off), /COMPRESS
writefits, imdir + name + '_exp_aligned_H.fits', image_shift(wt,x_off,y_off), /COMPRESS

; Determine offsets for Ks
; Apply offsets to lists and images
; and write them to file
; --------------------------
photdir = photdirk
imdir = imdirk
readcol, photdir + 'astrophot.txt', x, y, f, sx, sy, sf, c, ndet
compare_lists, xh, yh, x, y, x1c, y1c, x2c, y2c, $
           MAX_DISTANCE = max_distance, $
           SUBSCRIPTS_1 = subc1, SUBSCRIPTS_2 = subc2
RESISTANT_Mean, (x1c-x2c), 3.0, x_off, x_off_Sigma, Num_Rej, goodvec = goodvec
RESISTANT_Mean, (y1c-y2c), 3.0, y_off, y_off_Sigma, Num_Rej, goodvec = goodvec
x = x + x_off
y = y + y_off
print, 'Offsets in x and y of Ks relative to H: ' + strn(x_off) + ' +- ' + strn(x_off_Sigma) +  ', ' + strn(y_off) + ' +- ' + strn(y_off_Sigma)
forprint, TEXTOUT=photdir + 'astrophot_aligned_H.txt', x, y, f, sx, sy, sf, c, ndet, FORMAT='(8F13.3)', $
              COMMENT='#x[pix]  y[pix]  f[ADU]  dx[pix]  dy[pix]  c0[corr]  output from align_fields.pro'

im = readfits(imdir + name + '_holo.fits.gz')
noise = readfits(imdir + name + '_noise.fits.gz')
wt = readfits(imdir + name + '_exp.fits.gz')
writefits, imdir + name + '_holo_aligned_H.fits', image_shift(im,x_off,y_off), /COMPRESS
writefits, imdir + name + '_noise_aligned_H.fits', image_shift(noise,x_off,y_off), /COMPRESS
writefits, imdir + name + '_exp_aligned_H.fits', image_shift(wt,x_off,y_off), /COMPRESS

END
