import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error

########################################################################################
path = 'wiro'
path2 = '_wiro'
date = 2459368
name = 'wiro_wiro_broadband_final.pdf'
title = 'HD-189733b 2021-06-02 WIRO w/ WIRO Derived Priors'#\nRosenthal Priors and Detrending'
detrend = 0
#########################################################################################

plt.rc('axes', labelsize=12)
plt.rc('figure', titlesize=16)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)

filters = ['i','r','g','u']
colors = ['r','y','g','b']
offsets = [-0.03, -0.02, -0.01, 0]

fig, (ax1,ax2) = plt.subplots(2, 1, figsize=(8,6), gridspec_kw={'height_ratios': [2, 1]}, sharex='all')
plt.subplots_adjust(hspace=.0)

for k in range(4):
	transit_data = pd.read_csv(r'/d/tel2/brock/Data/RBO/hd189733/20210602/exofast/' + path + '_' + filters[k] + path2 + '/hd189733.mcmc.model.transit_0.txt', sep=" ", header=None)
	df_transit = pd.DataFrame(transit_data)
	time_transit = np.array(df_transit[0]) - date
	light_transit = np.array(df_transit[1])


	raw_data = pd.read_csv(r'/d/tel2/brock/Data/RBO/hd189733/20210602/exofast/n20210602.Sloan' + filters[k] + '.WIRO.dat', sep=" ", header=None)
	df_raw = pd.DataFrame(raw_data)
	time_raw = np.array(df_raw[0]) - date
	light_raw = np.array(df_raw[1])
	error_raw = np.array(df_raw[2])


	model_data = pd.read_csv(r'/d/tel2/brock/Data/RBO/hd189733/20210602/exofast/' + path + '_' + filters[k] + path2+ '/hd189733.mcmc.detrendedmodel.transit_0.planet_0.txt', sep=" ", header=None)
	df_model = pd.DataFrame(model_data)
	time_model = np.array(df_model[0]) - date
	light_model = np.array(df_model[1]) + 1


	pretty_data = pd.read_csv(r'/d/tel2/brock/Data/RBO/hd189733/20210602/exofast/' + path + '_' + filters[k] + path2+ '/hd189733.mcmc.prettymodel.transit_0.planet_0.txt', sep=" ", header=None)
	df_pretty = pd.DataFrame(pretty_data)
	time_pretty = np.array(df_pretty[0]) - date
	light_pretty = np.array(df_pretty[1]) + 1


	residual_data = pd.read_csv(r'/d/tel2/brock/Data/RBO/hd189733/20210602/exofast/' + path + '_' + filters[k] + path2+ '/hd189733.mcmc.residuals.transit_0.txt', sep=" ", header=None)
	df_residual = pd.DataFrame(residual_data)
	time_residual = np.array(df_residual[0]) - date
	residuals = np.array(df_residual[1])
	res_max = np.maximum(-1*residuals.min(), residuals.max())*1.25

	rms = mean_squared_error(residuals, np.zeros(np.size(residuals)), squared=False)*1e6

	#ax1.plot(time_model, light_model + offsets[k], color = colors[k], linestyle = '-', zorder = 0, label = filters[k] + '\': '+ str(round(rms,2)) + ' ppm')
	ax1.plot(time_pretty, light_pretty + offsets[k], color = colors[k], linestyle = '-', zorder = 0, label = filters[k] + '\': '+ str(round(rms,2)) + ' ppm')
	if detrend==1:
		ax1.plot(time_transit, light_transit + offsets[k], color = 'b', linestyle = '-', zorder = 5, label=filters[k] + '\' Detrended')
	ax1.plot(time_raw, light_raw + offsets[k], color = 'k', linestyle = '', marker='.', markersize=5 , zorder = 0)

	ax2.plot(time_residual, residuals + offsets[k]/2, color = 'k', linestyle = '', marker='.', markersize=5 , zorder = 0, label = filters[k] + '\': ' + str(round(rms,2)) + ' ppm')
	ax2.axhline(0 + offsets[k]/2, color=colors[k], linestyle='dotted', zorder=5, alpha=0.5, linewidth=2)


ax1.set_ylabel('Normalized Flux')
ax1.tick_params(axis='both', direction='in', which='both')
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles[::-1], labels[::-1], loc = 'best')
ax1.set_title(title)
ax1.set_xlim(left=time_residual.min()-0.0075,right=time_residual.max()+0.0075)
ax2.set_xlabel('BJD - ' + str(date) + ' (Days)')
ax2.set_ylabel('Residuals')
#ax2.set_ylim(top=res_max*1.5,bottom=-res_max)
#ax2.legend(loc='best', handles[::-1], labels[::-1])
ax2.tick_params(axis='both', direction='in', which='both')
fig.tight_layout()
plt.savefig('/d/tel2/brock/Data/RBO/hd189733/20210602/exofast/lightcurves/' + name)
plt.show()
