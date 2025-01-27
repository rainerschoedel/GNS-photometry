# Photometry pipeline for GNS-1


Astro-photometry in the new pipeline is done on the individual holographically reduced sub-images. Lists are subsequently merged and a final mosaic image is created.

1) ASTROPHOT, field_nr, filter, chip_nr


EXTRACTPSF.PRO is needed and called by ASTROPHOT.

The correlation threshold is currently set to 0.7
The detection threshold is set to two iterations at 5 sigma. 
Going deeper will slow things down and will probably not help much because we will run into the crowding limit anyway.


2) COMPUTE_UNCERTAINTIES, field_nr, filter, chip_nr

Uses output from ASTROPHOT as input. 
Eliminates all stars that are not detected in all images (deep image and jackknife images).
Computes astrometric and photometric jackknife uncertainties. 


WARNING ON TEMPORARY FILES
Temporary files will be stored in the Band/Field/tmp/tmpX
directory. Temporary files WILL GENERALLY BE REMOVED when a new script is started. 

Results and some intermediate files will be stored in the directories:
year/band/field/photo/chipX/lists


3) MERGE_SUBLISTS, field_nr, Band, chip_nr

Takes input from ALIGN_FIELDS

Merge all the sublists: Determine and apply x/y offsets for each sub-field; combine measurements and uncertainties; compute PSF uncertainties from multiple measurements; estimate PSF uncertainty for stars without multiple measurements.

OUTPUT will be written to year/band/field/photo/chipX/lists/astrophot.txt

Some plots are written to year/band/field/photo/chipX/plots (Statistical and PSF uncertainties).

4) merge_subims, field_nr, chip_nr, Band

Merges sub-images to create final mosaics.
Does to images what merge_lists.pro does to star lists

OUTPUT will be written to [year/band/field]/cubims/chip/





5) ALIGN_FIELDS, field_nr, chip_nr

This scripts aligns the lists and images produced by merge_sublists and merge-subims all lists with the H-band list and images. Only a simple shift (typically a fraction of a pixel) is computed and applied. 

output for H-band: astrophot.txt (unchanged)
output for J, Ks: astrophot_aligned_H.txt


6) calibrate, field, band, chip

Calibrates star lists  astrometrically (with VVV) and photometrically (with SIRIUS).

You will have to manually .cont the code after it has looked for alignment stars. If it does not find several hundreds to > 1000 alignment stars, then there is probably some trouble.

t_exp must be set correctly (usually 1.26s, but there are observations from the begining of GNS that have 0.85s, I think).

INPUT is output from align_fields, i.e. the files  
year/band/field/photo/chipX/lists/astrophot.txt or astrophot_aligned_H.txt
The VVV data are located in /home/data/VVV/
SIRIUS in /home/data/SIRIUS

OUTPUT 
Calibrated star lists (stars_calibrated_[band]_chip[nr].txt ) and map of calibrated stars in in/ home/data/GNS/2015/‘ + band + '/' + strn(field) + '/photo/chip' + strn(chip) + ‘/lists/calibrated_' + band + '_' + strn(field) + '_' + strn(chip) + '.txt'

Plots in /home/data/GNS/2015/‘ + band + '/' + strn(field) + '/photo/chip' + strn(chip) + ‘/lists/'
'position_uncertainty_pixels_' + strn(chip) + ‘.eps' (distributions of differences of
astrometric  positions between reference stars in HAWK-I and VVV) measured in image coordintes
'dpos_calibrated_' + strn(chip) + ‘.eps' (distributions of differences of CALIBRATED astrometric  positions between reference stars in HAWK-I and VVV) measured in arcseconds and actual coordinates after calibration
'ZP_hist_' + band + ‘.eps'  - histogram of zero points
'ZP_calib_' + band + ‘.eps' - scatter plot of zero points
'uncertainty_phot_' + strn(chip) + ‘.eps' - photometric uncertainty
'uncertainty_xpos_' + strn(chip) + ‘.eps' - astrometric uncertainty in x
'uncertainty_ypos_' + strn(chip) + ‘.eps' - astrometric uncertainty in y

Astrometrically calibrated images/uncertainty/exposure maps in /
res_dir_im = '/home/data/GNS/2015/' + band + '/' + strn(field) + '/cubeims/chip' + strn(chip) + '/'







3) Cleanup

Delete raw data:
dark
flat
sky
science

Delete intermediate data products in 
tmp and sub-directories
cleaned
aligned
ims: *cube*.fits.gz
cubes: rm cube*, rm mask*, rm lnx*
subcubes/chip?/holomask*cd 

We may want to keep the files in cubeims (holographic reconstructions of sub-fields), so that we can re-run the phot0-astrometry if necessary.

We can also keep the subcubes and masks in subcubes/chip? (~20 GB per chip).
Then we can repeat holography when necessary.
Avoid deleting any other files. In particular, do not delete the 'xy_off_f'+strn(field_nr)+'.SAV'
files created by alignquadrants (otherwise, a lot of the process needs to be repeated if it turns out to be necessary).

