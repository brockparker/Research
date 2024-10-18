import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import CubicSpline
from scipy.ndimage import gaussian_filter

import matplotlib.cm as cm
import matplotlib.colors as colors
import extinction
from matplotlib.lines import Line2D

from dust_extinction.parameter_averages import F19, CCM89, VCG04, G23

######################################################################################################
lam1, lam2, cdelt1 = 	1300, 25000, 0.1		# wavelength range
subsample =		5
red = 			[2.3, 3.1, 5.0]
exts = 			[1.0]
colors = 		['tab:red', 'tab:green','tab:blue']
######################################################################################################

def CCM89cust(lamda, AV, RV):
	inverse = np.array(1./lamda)
	IR = np.array(np.nonzero(inverse < 1.1)).ravel()
	Optical = np.array(np.nonzero((inverse >= 1.1) & (inverse <= 3.3))).ravel()
	UV = np.array(np.nonzero((inverse > 3.3) & (inverse <= 8.))).ravel()

	if len(IR) > 0:
		a_IR = 0.574 * (inverse[IR]**1.61)
		b_IR = -0.527 * (inverse[IR]**1.61)
		lam_IR = (a_IR + b_IR/RV) * AV
	else:
		lam_IR = np.array([])

	if len(Optical) > 0:
		y = inverse[Optical]-1.82
		a_opt = 1. + 0.17699*y - 0.50447*(y**2) - 0.02427*(y**3) + 0.72085*(y**4) + 0.01979*(y**5) - 0.77530*(y**6) + 0.32999*(y**7)
		b_opt = 1.41338*y + 2.28305*(y**2) + 1.07233*(y**3) - 5.38434*(y**4) - 0.62251*(y**5) + 5.30260*(y**6) - 2.09002*(y**7)
		lam_opt=(a_opt + b_opt/RV) * AV
	else:
		lam_opt=np.array([])
	
	if len(UV) > 0:
		UV_2 = np.array(np.nonzero((inverse >= 5.9) & (inverse <= 8.))).ravel()
		UV_1 = np.array(np.nonzero((inverse < 5.9) & (inverse > 3.3))).ravel()

		Fa_2 = -0.04473*(inverse[UV_2]-5.9)**2 - 0.009779*(inverse[UV_2]-5.9)**3
		Fb_2 = 0.2130*(inverse[UV_2]-5.9)**2 + 0.1207*(inverse[UV_2]-5.9)**3
		
		Fa_1 = UV_1 - UV_1	# 0
		Fb_1 = Fa_1		# 0

		F_a = np.concatenate([Fa_2, Fa_1])
		F_b = np.concatenate([Fb_2, Fb_1])

		a_UV = 1.752 - 0.316*inverse[UV] - 0.104/((inverse[UV]-4.67)**2 + 0.341) + F_a
		b_UV = -3.090 + 1.825*inverse[UV] + 1.206/((inverse[UV]-4.62)**2 + 0.263) + F_b
   
		lam_UV = (a_UV + b_UV/RV) * AV
	else:
		lam_UV=np.array([])

	A_lam = np.concatenate([lam_UV,lam_opt,lam_IR])

	return A_lam

def VCG04(microns, Av, Rv):
	inverse = np.array(1./microns)
	IR = np.array(np.nonzero(inverse < 1.1)).ravel()
	Optical = np.array(np.nonzero((inverse >= 1.1) & (inverse <= 3.3))).ravel()
	UV = np.array(np.nonzero((inverse > 3.3) & (inverse <= 8.))).ravel()

	if len(IR) > 0:
		a_IR = 0.574 * (inverse[IR]**1.61)
		b_IR = -0.527 * (inverse[IR]**1.61)
		lam_IR = (a_IR + b_IR/Rv) * Av
	else:
		lam_IR = np.array([])

	if len(Optical) > 0:
		y = inverse[Optical]-1.82
		a_opt = 1. + 0.17699*y - 0.50447*(y**2) - 0.02427*(y**3) + 0.72085*(y**4) + 0.01979*(y**5) - 0.77530*(y**6) + 0.32999*(y**7)
		b_opt = 1.41338*y + 2.28305*(y**2) + 1.07233*(y**3) - 5.38434*(y**4) - 0.62251*(y**5) + 5.30260*(y**6) - 2.09002*(y**7)
		lam_opt=(a_opt + b_opt/Rv) * Av
	else:
		lam_opt=np.array([])
	
	if len(UV) > 0:
		UV_2 = np.array(np.nonzero((inverse >= 5.9) & (inverse <= 8.))).ravel()
		UV_1 = np.array(np.nonzero((inverse < 5.9) & (inverse > 3.3))).ravel()

		Fa_2 = -0.0077*(inverse[UV_2]-5.9)**2 - 0.0030*(inverse[UV_2]-5.9)**3
		Fb_2 = 0.2060*(inverse[UV_2]-5.9)**2 + 0.0550*(inverse[UV_2]-5.9)**3
		
		Fa_1 = UV_1 - UV_1	
		Fb_1 = Fa_1		

		F_a = np.concatenate([Fa_2, Fa_1])
		F_b = np.concatenate([Fb_2, Fb_1])

		a_UV = 1.708 - 0.215*inverse[UV] - 0.134/((inverse[UV]-4.558)**2 + 0.566) + F_a
		# a_UV = 1.808 - 0.215*inverse[UV] - 0.134/((inverse[UV]-4.558)**2 + 0.566) + F_a
		b_UV = -2.350 + 1.403*inverse[UV] + 1.103/((inverse[UV]-4.587)**2 + 0.263) + F_b
   
		lam_UV = (a_UV + b_UV/Rv) * Av
	else:
		lam_UV=np.array([])

	A_lam = np.concatenate([lam_UV,lam_opt,lam_IR])

	return A_lam


lamA1 = np.arange(lam1,lam2,cdelt1)	
finalwavs=lamA1[0::subsample]

plt.rc('axes', labelsize=16)
plt.rc('figure', titlesize=16)
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('legend', fontsize=18)

fig, ax = plt.subplots(figsize=(8,6), layout = 'tight')
ax2 = ax.twiny()
tick_loc = [0.25, 2, 5]
tick_labels = ['IR', 'Visible', 'UV']

for i in range(len(red)):

	CCM89func = CCM89(red[i])
	CCM89F = CCM89func.extinguish(1/finalwavs*1e4, Av =  1.0)
	CCM89ext = 2.5 * -np.log10(CCM89F)
	ax.plot(1/finalwavs*1e4, CCM89ext, color=colors[i],label = 'CCM89', linestyle = ':')

	F19func = F19(red[i])
	F19F = F19func.extinguish(1/finalwavs*1e4, Av =  1.0)
	F19ext = 2.5 * -np.log10(F19F)
	ax.plot(1/finalwavs*1e4, F19ext, color=colors[i],label = 'F19', linestyle = '-')

	G23func = G23(red[i])
	G23F = G23func.extinguish(1/finalwavs*1e4, Av =  1.0)
	G23ext = 2.5 * -np.log10(G23F)
	ax.plot(1/finalwavs*1e4, G23ext, color=colors[i],label = 'G23', linestyle = '-.')

	VCG04ext = VCG04(finalwavs/1e4, 1.0, red[i])
	#VCG04F = 2.5**VCG04ext
	#ax.plot(1/finalwavs*1e4, VCG04ext, color = colors[i], label = 'VGC04', linestyle='--')

#plt.text(0.25, 4.5, '$R_V$ = 2.0', c = colors[0], fontsize = 14)
ax.text(5.25, 4, '$R_V$ = 2.3', c = colors[0], fontsize = 20)
#ax.text(5.25, 4.35, 'Small Dust', c = 'k', fontsize = 16)
ax.text(5.25, 2.1, '$R_V$ = 3.1', c = colors[1], fontsize = 20)
ax.text(5.25, 1.1, '$R_V$ = 5.0', c = colors[2], fontsize = 20)
#ax.text(5.25, 0.8, 'Large Dust', c = 'k', fontsize = 16)

ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(tick_loc)
ax2.set_xticklabels(tick_labels)

ccm_patch = Line2D([0], [0], color='k',linestyle=':' , label='CCM89')
f_patch = Line2D([0], [0], color='k',linestyle='-', label='F19')
#vcg_patch=Line2D([0], [0], color='k', linestyle='--',label='VCG04')
g_patch=Line2D([0], [0], color='k', linestyle='-.',label='G23')
ax.legend(handles=[ccm_patch,f_patch,g_patch], loc = 'upper left')

#plt.xscale('log')
ax.set_xlabel(r'Inverse Wavelength ($\mu$m$^{-1}$)')
ax.tick_params(axis='both', direction='in', which='both')
#plt.title('Extinction Curve Comparison')
ax.set_ylabel('Extinction (mag)')
plt.tight_layout()
plt.show()
