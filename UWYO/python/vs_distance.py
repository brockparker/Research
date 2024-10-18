import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Plot some as upper limits
# add text to lable each star

stars = ['04','13','24','26','27', 'Zeta Oph']
dist = [223, 449, 202, 669, 499, 135]
dist_err_lower = [5, 8, 21, 19, 5, 11]
dist_err_upper = [9, 20, 6, 18, 6, 13]

av_ccm = [0.86, 1.31, 0.67, 1.15, 1.05, 0.93]
av_f = [0.88, 1.33, 0.65, 1.17, 1.01, 1.03]
av_g = [0.76, 1.26, 0.56, 1.11, 0.99, 0.77]

rv_ccm = [1.99, 2.37, 1.47, 2.48, 2.18, 2.86]
rv_f = [2.01, 2.35, 1.37, 2.50, 2.06, 3.11]
rv_g = [1.79, 2.38, 1.26, 2.48, 2.09, 2.43]

av_ccm_err = [0.07, 0.07, 0.07, 0.08, 0.06, 0.07]
av_f_err = [0.06, 0.06, 0.07, 0.07, 0.05, 0.06]
av_g_err = [0.07, 0.07, 0.07, 0.08, 0.06, 0.09]

rv_ccm_err = [0.19, 0.09, 0.08, 0.11, 0.10, 0.17]
rv_f_err = [0.22, 0.11, 0.09, 0.11, 0.12, 0.14]
rv_g_err = [0.24, 0.11, 0.10, 0.12, 0.12, 0.24]

colors = ['tab:blue', 'tab:green', 'tab:red', 'tab:orange', 'tab:cyan', 'tab:purple', 'tab:olive', 'tab:pink']

plt.rc('axes', labelsize=16)
plt.rc('figure', titlesize=16)
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('legend', fontsize=16)

fig, (ax1,ax2) = plt.subplots(2, 1, figsize=(8,6), gridspec_kw={'height_ratios': [1, 1], 'wspace':0,'hspace':0}, sharex='all', layout='tight')

for i in range(len(stars)):
	ax1.errorbar(dist[i], av_ccm[i], yerr = av_ccm_err[i], xerr = [dist_err_lower[i:i+1], dist_err_upper[i:i+1]], marker = 's', capsize = 0, linestyle = '', color = colors[i], markersize=12, alpha=0.5)
	ax1.errorbar(dist[i], av_f[i], yerr = av_f_err[i], xerr = [dist_err_lower[i:i+1], dist_err_upper[i:i+1]], marker = '^', capsize = 0, linestyle = '', color = colors[i], markersize=12, alpha=0.5)
	ax1.errorbar(dist[i], av_g[i], yerr = av_g_err[i], xerr = [dist_err_lower[i:i+1], dist_err_upper[i:i+1]], marker = 'P', capsize = 0, linestyle = '', color = colors[i], markersize=12, alpha=0.5)

	ax2.errorbar(dist[i], rv_ccm[i], yerr = rv_ccm_err[i], xerr = [dist_err_lower[i:i+1], dist_err_upper[i:i+1]], marker = 's', capsize = 0, linestyle = '', color = colors[i], markersize=12, alpha=0.5)
	ax2.errorbar(dist[i], rv_f[i], yerr = rv_f_err[i], xerr = [dist_err_lower[i:i+1], dist_err_upper[i:i+1]], marker = '^', capsize = 0, linestyle = '', color = colors[i], markersize=12, alpha=0.5)
	ax2.errorbar(dist[i], rv_g[i], yerr = rv_g_err[i], xerr = [dist_err_lower[i:i+1], dist_err_upper[i:i+1]], marker = 'P', capsize = 0, linestyle = '', color = colors[i], markersize=12, alpha=0.5)




# =============================================================================
# 	if rv_f[i]<1:
# 		el = ax2.errorbar(dist[i], rv_f[i], xerr = [dist_err_lower[i:i+1], dist_err_upper[i:i+1]], marker = '^', capsize = 0, linestyle = '', color = colors[i], markersize=12,uplims=[1], yerr=0.2, alpha=0.5)
# 		elines = el.get_children()
# 		elines[1].set_color('k')
# 		elines[3].set_color('k')
# =============================================================================



ax2.set_xlabel('Distance (pc)')
ax1.set_ylabel('$A_V$ (mag)')
ax2.set_ylabel('$R_V$')
ax1.tick_params(axis='both', direction='in', which='both')
ax2.tick_params(axis='both', direction='in', which='both', top=True)

ccm_patch = Line2D([], [], color='k',linestyle='-' , marker='s', label='CCM89', markersize=10)
f_patch = Line2D([], [], color='k',linestyle='-', marker='^', label='F19', markersize=10)
#vcg_patch=Line2D([], [], color='k', linestyle='-', marker='d', label='VCG04', markersize=10)
g_patch=Line2D([], [], color='k', linestyle='-', marker='P', label='G23', markersize=10)

ax1.text(120, 0.60, 'ZOph', c = 'k', fontsize = 14)
ax1.text(180, 0.75, '24', c = 'k', fontsize = 14)
ax1.text(225, 0.65, '04', c = 'k', fontsize = 14)
ax1.text(440, 1.12, '13', c = 'k', fontsize = 14)
ax1.text(490, 0.85, '27', c = 'k', fontsize = 14)
ax1.text(659, 0.96, '26', c = 'k', fontsize = 14)

#ax3 = ax1.twiny()
#tick_loc = list(np.array(dist) + np.array([+5, 0, 0, 0, 0, -10]))
#ha = ['right','center','left','center','center','center']
#tick_labels = stars
#ax3.set_xlim(ax1.get_xlim())
#ax3.set_xticks(tick_loc)
#ax3.set_xticklabels(tick_labels, ha='left')

ax1.legend(fancybox=True, loc = 'upper left', ncol=2, handlelength=1, handles=[ccm_patch, f_patch, g_patch])

#ax1.set_title('Extinction vs Distance')
plt.savefig("/d/users/brock/paper_figures/vs_distance.png")
plt.show()
