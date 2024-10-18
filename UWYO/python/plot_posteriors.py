"""
Script to plot the outputs of the pseudo-MC code posteriors for Av and Rv.
"""

import numpy as np
import matplotlib.pyplot as plt
import corner
import os
os.chdir('/d/users/brock/python/')
import extinction
# BP Necessary imports.

plt.rc('axes', labelsize=16)
plt.rc('figure', titlesize=16)
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('legend', fontsize=16)

red = extinction.reddening()

path = '/d/tel2/brock/Data/Extinction/'
# BP Defining file locations.

names = ['pk04', 'pk13', 'pk24', 'pk26', 'pk27']
#names = ['pk04']
# BP Defining target names.
ext_models = ['CCM89', 'F19', 'G23']
#ext_models = ['CCM89']
# BP Defining extinction models.

for name in names:
	for ext_model in ext_models:
		data = np.loadtxt(path + name + '_posterior_' + ext_model + '.txt')
		# BP Reading in data
		
		temps = data[:, 0]
		gravs = data[:, 1]
		metals = data[:, 2]
		chis = data[:, 3]
		avs = data[:, 4]
		rvs = data[:, 5]
		# BP Defining each of the parameters
				
		posteriors = np.vstack((temps, metals, avs, rvs)).T
		# BP Adding data into correct array shape for corner plotting.
		
		#print(np.shape(posteriors))
		
		figure = corner.corner(
			posteriors,
			labels=[
				r"$T_{{eff}}$",
				r"[Fe/H]",
				r"$A_V$",
				r"$R_V$",
				r"$\chi ^2$"
			],
			levels=[0.39, 0.86],
			quantiles=[0.16, 0.5, 0.84],
			smooth=0.75,
			bins=40,
			show_titles=True,
			title_kwargs={"fontsize": 15})
		# BP Creating corner plot.
		figure.gca().annotate(name.replace('pk', 'PK-') + ': ' + ext_model, size=16, xy=(0.085,0.875), xycoords='figure fraction',xytext=(-20,-10), textcoords='offset points', ha='center', va='center', rotation=90)
		# BP Adding target and extinction model to figures.
		
		plt.savefig('/d/users/brock/paper_figures/' + name + '_posterior_' + ext_model + '.pdf')
		plt.show()
		# BP Saving figures.