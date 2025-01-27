PRO COMBINE_CHIPS, field, band

; PURPOSE: Last step of astro-photometry
;          Combine lists from the four chips of each band
;
; INPUT:   Calibrated lists for each chip and band
;          created by CALIBRATE

dmax = 0.05 ; max distance for list comparison in arcseconds
chip = 1

phot_dir = '/data/GNS/2015/' + band + '/' + strn(field) + '/photo/chip' + strn(chip) + '/lists/'
readcol, phot_dir + 'stars_calibrated_' + band + '_chip' + strn(chip) + '.txt', a ,d, x, y, f, m, sx, sy, sf, sm
; convert RA and Dec to arcseconds (double)
a = 3600. * double(a)
d = 3600. * double(d)

for chip = 2, 4 do begin

  phot_dir = '/data/GNS/2015/' + band + '/' + strn(field) + '/photo/chip' + strn(chip) + '/lists/'
  readcol, phot_dir + 'stars_calibrated_' + band + '_chip' + strn(chip) + '.txt', a2, d2, x2, y2, f2, m2, sx2, sy2, sf2, sm2
  ; convert RA and Dec to arcseconds (double)
  a2 = 3600. * double(a2)
  d2 = 3600. * double(d2)

  ; compare list for this chip with current list of 
  ; all detected stars
  ; -----------------------------------------------------
  if (n_detected gt 0) then begin
    x1 = xl
    y1 = yl
  endif else begin
    x1 = [0.0]
    y1 = [0.0]
  endelse
  x2 = x
  y2 = y

  compare_lists, a1, d1, a2, d2, a1c, d1c, a2c, d2c, MAX_DISTANCE=tolerance, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
  nc = n_elements(subc2)
  if (i_x eq 0) and (i_y eq 0) then nc = 0L
  print, 'Number of common stars: ' + string(nc)
  if (sub1[0] gt (-1)) then n_new = n_elements(sub2)
  print, 'Number of new stars: ' + string(n_new)


  ; save values of new common stars in common list
  ; i.e. create additional entries (measurements) for
  ; stars already in the common list
  ; ------------------------------------------------

  for jj = 0L, (nc-1) do begin
    n_all[subc1[jj]] = n_all[subc1[jj]] + 1.0
    ind = n_all[subc1[jj]] - 1
    x_all[ind,subc1[jj]] = x[subc2[jj]]
    y_all[ind,subc1[jj]] = y[subc2[jj]]
    f_all[ind,subc1[jj]] = f[subc2[jj]]
    sx_all[ind,subc1[jj]] = sx[subc2[jj]]
    sy_all[ind,subc1[jj]] = sy[subc2[jj]]
    sf_all[ind,subc1[jj]] = sf[subc2[jj]]
    c_all[ind,subc1[jj]] = c[subc2[jj]]
  endfor

  ; save values of new stars in common list
  ; i.e. create entries for new stars
  ; -----------------------------------------

   for jj = 0L, (n_new-1) do begin
      n_all[(n_detected+jj)] = 1.0
      x_all[0,(n_detected+jj)] = x[sub2[jj]]
      y_all[0,(n_detected+jj)] = y[sub2[jj]]
      f_all[0,(n_detected+jj)] = f[sub2[jj]]
      sx_all[0,(n_detected+jj)] = sx[sub2[jj]]
      sy_all[0,(n_detected+jj)] = sy[sub2[jj]]
      sf_all[0,(n_detected+jj)] = sf[sub2[jj]]
      c_all[0,(n_detected+jj)] = c[sub2[jj]]
    endfor

    ; Now create a list with the mean positions
    ; of already detected stars
    n_detected = n_detected + n_new ; current total number of detected stars
    nl = n_all[0:n_detected-1]
    xl = fltarr(n_detected)
    yl = fltarr(n_detected)
    fl = fltarr(n_detected)
    for i_star = 0, n_detected-1 do begin    
      vals = f_all[*,i_star]
      good = where(vals gt 0)
      fl[i_star] = mean(vals[good])
      vals = x_all[*,i_star]
      xl[i_star] = mean(vals[good])
      vals = y_all[*,i_star]
      yl[i_star] = mean(vals[good])
    endfor

  endfor


endfor

END
