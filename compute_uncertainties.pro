pro compute_uncertainties, field, filter, chip_nr


; PURPOSE: Take the input lists from astrophot.pro
;          Compare deep list with jackknife lists.
;          Only acccept stars common to all.
;          Compute the jackknife uncertainties.
;          Export new lists with measurements and uncertainties.
;



chip = 'chip' + strn(chip_nr)
max_distance = 2 ; Correponds to dmax in search for 
              ; common stars between a subimmages
              ;  for fine-alignment of subimages
              ; with rebfac = 2, tolerance = 2
              ; corresponds to half the resolution limit
 
; Numbers of first and last sub-images (usually not to be changed, see
; runholo.pro)
i_x0 = 0
i_x1 = 6

i_y0 = 0
i_y1 = 2

; number of sub-images
n_jack = 10

; input and output directories
basedir = '/home/data/GNS/2015/'

photdir = basedir + filter + '/' + field + '/photo/chip' + strn(chip_nr) + '/lists/'

for i_x = i_x0, i_x1 do begin
  for i_y = i_y0, i_y1 do begin 

     name = strn(i_x) + '_' +  strn(i_y)
     readcol, photdir + 'stars_' + name + '.txt', x0, y0, f0, sx0, sy0, sf0, c0
     n_deep = n_elements(f0)
     print, 'Deep image contains ' + strn(n_deep) + ' point sources.'

     ; Find stars common to all jackknife images
     ; Compare each list with the deep list 
     ; ----------------------------------------
     for i_j = 1, n_jack do begin
       readcol, photdir + 'stars_' + name + '_s' + strn(i_j) + '.txt', x, y, f, sx, sy, sf, c, /SILENT
       compare_lists, x0, y0, x, y, x1c, y1c, x2c, y2c, $
           MAX_DISTANCE = max_distance, $
           SUBSCRIPTS_1 = subc1, SUBSCRIPTS_2 = subc2
       x0 = x0[subc1]
       y0 = y0[subc1]
       f0 = f0[subc1]
       c0 = c0[subc1]
;       print, 'Jackknife image ' + strn(i_j) + ':'
;       print, 'Common stars with deep image: ' + strn(n_elements(subc1))
     endfor
     n_detected = n_elements(subc1)
     print, 'Number of stars common to all jackknife images: ' + strn(n_detected)

     ; Gather all measurements of common stars
     ; ----------------------------------------
     x_jack = fltarr(n_detected, n_jack)
     y_jack = fltarr(n_detected, n_jack)
     f_jack = fltarr(n_detected, n_jack)
     for i_j = 1, n_jack do begin
       readcol, photdir + 'stars_' + name + '_s' + strn(i_j) + '.txt', x, y, f, sx, sy, sf, c, /SILENT
       compare_lists, x0, y0, x, y, x0c, y0c, xc, yc, $
           MAX_DISTANCE = max_distance, $
           SUBSCRIPTS_1 = subc1, SUBSCRIPTS_2 = subc2
;       print, 'Jackknife image ' + strn(i_j) + ':'
;       print, 'Common stars with deep image: ' + strn(n_elements(subc1))
       x_jack[*,i_j-1] = x[subc2]
       y_jack[*,i_j-1] = y[subc2]
       f_jack[*,i_j-1] = f[subc2]
     endfor

    ; Now compute unbiased estimators and uncertainties
    xub = fltarr(n_detected)
    yub = fltarr(n_detected)
    fub = fltarr(n_detected)
    dx = fltarr(n_detected)
    dy = fltarr(n_detected)
    df = fltarr(n_detected)
    nnjack = (n_jack-1)/float(n_jack)
    for s = 0, n_detected-1 do begin
      xub[s] = mean(x_jack[s,*])
      yub[s] = mean(y_jack[s,*])
      fub[s] = mean(f_jack[s,*])
      dx[s] = sqrt(total((x_jack[s,*]-xub[s])^2)*nnjack)
      dy[s] = sqrt(total((y_jack[s,*]-yub[s])^2)*nnjack)
      df[s] = sqrt(total((f_jack[s,*]-fub[s])^2)*nnjack)
    endfor

    forprint, TEXTOUT=photdir + 'list_' + name + '.txt', x0, y0, f0, dx, dy, df, c0, FORMAT='(7F13.3)', $
              COMMENT='#x[pix]  y[pix]  f[ADU]  dx[pix]  dy[pix]  c0[corr]  output from compute_uncertainties.pro'
    print, 'Finished filter ' + filter + ', field ' + field + ' , chip ' + chip + '.'

  endfor
endfor

END
