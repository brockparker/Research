'''Quick script to correct Gaia GSP-Phot [M/H]'''

import gdr3apcal
import pandas as pd
# BP Necessary imports

names = ['PK-04','PK-13','PK-24','PK-26','PK-27']

targets = [ 
	[[5721.588, 4.1274, -0.2181, 1.0033, 0.4399, 0.8136, 4, 'MARCS', 249.10352805285947, -10.544322055137767]],
	[[6633.719, 3.9449, -0.2314, 1.6174, 0.7235, 1.3338, 3, 'PHOENIX', 249.23674301875465, -10.48716639328424]],
	[[5090.673, 4.3066, -0.8000, 0.4066, 0.1763, 0.3261, 6, 'MARCS', 249.46164303415972, -10.551895829331551]],
	[[6731.061, 3.8976, -0.4726, 1.4831, 0.6765, 1.2400, 3, 'MARCS', 249.50491920428755, -10.579727553285723]],
	[[6671.852, 4.1343, -0.6076, 1.3251, 0.6064, 1.1116, 3, 'A', 249.10352805285947, -10.544322055137767]]]
# BP Creating list of Gaia GSP-Phot parameters.

cols = ['teff_gspphot', 'logg_gspphot', 'mh_gspphot', 'azero_gspphot', 'ebpminrp_gspphot', 'ag_gspphot', 'mg_gspphot', 'libname_gspphot', 'ra', 'dec']
# BP Definding pandas column names.

for i in range(len(names)):
	df = pd.DataFrame(targets[i], columns=cols)
	# BP Creating data frame object with proper columns.

	calib = gdr3apcal.GaiaDR3_GSPPhot_cal()
	# BP initializing calibration object.
	metal_calib = calib.calibrateMetallicity(df)
	#teff_calib = calib.calibrateTeff(df)
	# BP Calibrating metallicity.
	
	print('{}'.format(names[i]))
	
	print('The orignial metallicity is {:.3f}.\
		  \nThe new calibrated metallicity is {:.3f}'.format(targets[i][0][2], metal_calib[0]))