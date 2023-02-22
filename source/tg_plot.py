import numpy as np
import matplotlib.pyplot as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.lines as lines
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from source.tg_util import posmax

#
#
#
def getColor(i, N, colormap='jet'):
	cm = mpl.get_cmap(colormap) 
	cNorm  = colors.Normalize(vmin=0, vmax=N+1)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
	colorVal = scalarMap.to_rgba(i)
	return colorVal

#
#
#
def get_read_alignment_polygons(alignments, readlen):
	polygons = []
	p_alpha  = []
	p_color  = []
	p_text   = []
	yp       = [-0.05 * readlen, 0.05 * readlen]
	text_y   = 0.7 * yp[1]
	for i in range(len(alignments)):
		n = alignments[i]
		#
		my_refname = n[2]
		my_refspan = str(n[3])+' - '+str(n[4])
		if my_refname[:3] == 'tel':
			if my_refname[-1] == '?':	# skip these weird regions
				p_color.append('gray')
			else:
				my_refname = 'tel'
				p_color.append('red')
		else:
			p_color.append('blue')
		p_alpha.append(float(n[6]+15.)/(60.+15.))
		#
		xp               = [n[0], n[1]]
		delta_pointy     = 0.8*(n[1]-n[0])
		delta_pointy_rev = 0.2*(n[1]-n[0])
		if n[5] == 'FWD':
			polygons.append(Polygon(np.array([[xp[0],yp[0]], [xp[0],yp[1]], [xp[0]+delta_pointy,yp[1]], [xp[1],0.], [xp[0]+delta_pointy,yp[0]]]), closed=True))
		else:
			polygons.append(Polygon(np.array([[xp[0]+delta_pointy_rev,yp[0]], [xp[0],0.], [xp[0]+delta_pointy_rev,yp[1]], [xp[1],yp[1]], [xp[1],yp[0]]]), closed=True))
		
		p_text.append((n[0], text_y, my_refname+' : '+my_refspan))
		text_y -= 0.40*yp[1]
		#
	axis_val = [0, readlen+1, 2.5*yp[0], 2.5*yp[1]]
	#
	return (polygons, p_color, p_alpha, p_text, axis_val)

#
#
#
def plot_all_read_data(density_data, tl_vals, aln_dat, tel_window, f_title, fig_name, plot_for_paper=False):
	[td_p_e0, td_p_e1, td_q_e0, td_q_e1, p_vs_q_power] = density_data
	[tl_p, tl_q] = tl_vals
	rlen = len(aln_dat[0][7])
	#
	x_off   = tel_window/2
	dens_yt = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
	dens_yl = ['0%', '20%', '40%', '60%', '80%', '100%']
	pred_yt = [-1.0, -0.5, 0.0, 0.5, 1.0]
	#
	if plot_for_paper:
		subplot_nums = [311, 312, 313]
		pred_yl = ['-1', '', '0', '', '1']
		mpl.rcParams.update({'font.size': 20, 'font.weight':'bold'})
		my_figsize = (12,8)
	else:
		subplot_nums = [411, 412, 413]
		pred_yl = ['(q) -1', '', '0', '', '(p) 1']
		mpl.rcParams.update({'font.size': 16, 'font.weight':'bold'})
		my_figsize = (12,10)
	#
	fig = mpl.figure(1, figsize=my_figsize)
	ax1 = mpl.subplot(subplot_nums[0])
	mpl.plot(np.arange(len(td_p_e1))+x_off, td_p_e1, '-m', alpha=0.5)
	mpl.plot(np.arange(len(td_p_e0))+x_off, td_p_e0, '-k', alpha=0.8)
	mpl.rcParams.update({'text.usetex': True})
	mpl.legend(['$p_1(i)$', '$p_0(i)$'], loc=1)
	mpl.rcParams.update({'text.usetex': False})
	mpl.xlim([0,rlen])
	mpl.yticks(dens_yt, dens_yl)
	mpl.ylim([0,1])
	mpl.grid(linestyle='--', alpha=0.5)
	mpl.ylabel('tel density (p)')
	if plot_for_paper == False:
		mpl.title(f_title)
	#
	mpl.subplot(subplot_nums[1], sharex=ax1)
	mpl.plot(np.arange(len(td_q_e1))+x_off, td_q_e1, '-m', alpha=0.5)
	mpl.plot(np.arange(len(td_q_e0))+x_off, td_q_e0, '-k', alpha=0.8)
	mpl.rcParams.update({'text.usetex': True})
	mpl.legend(['$q_1(i)$', '$q_0(i)$'], loc=1)
	mpl.rcParams.update({'text.usetex': False})
	mpl.yticks(dens_yt, dens_yl)
	mpl.ylim([0,1])
	mpl.grid(linestyle='--', alpha=0.5)
	mpl.ylabel('tel density (q)')
	#
	mpl.subplot(subplot_nums[2], sharex=ax1)
	mpl.plot(np.arange(len(p_vs_q_power))+x_off, p_vs_q_power, '-k')
	mpl.rcParams.update({'text.usetex': True})
	mpl.legend(['$S(i)$'], loc=1)
	mpl.rcParams.update({'text.usetex': False})
	mpl.grid(linestyle='--', alpha=0.5)
	mpl.ylabel('tel score')
	mpl.yticks(pred_yt, pred_yl)
	mpl.ylim([-1.0, 1.0])
	polygons = []
	p_color  = []
	p_alpha  = []
	if tl_p > 0:
		xp = [0, tl_p]
		yp = [-1, 1]
		polygons.append(Polygon(np.array([ [xp[0],yp[0]], [xp[0],yp[1]], [xp[1],yp[1]], [xp[1],yp[0]] ]), closed=True))
		p_color.append('red')
		p_alpha.append(0.5)
	if tl_q > 0:
		xp = [rlen - tl_q, rlen]
		yp = [-1, 1]
		polygons.append(Polygon(np.array([ [xp[0],yp[0]], [xp[0],yp[1]], [xp[1],yp[1]], [xp[1],yp[0]] ]), closed=True))
		p_color.append('red')
		p_alpha.append(0.5)
	if len(polygons):
		ax = mpl.gca()
		for i in range(len(polygons)):
			ax.add_collection(PatchCollection([polygons[i]], color=p_color[i], alpha=p_alpha[i]))
	#
	if plot_for_paper == False:
		mpl.subplot(414, sharex=ax1)
		(polygons, p_color, p_alpha, p_text, axis_val) = get_read_alignment_polygons(aln_dat, rlen)
		ax = mpl.gca()
		for i in range(len(polygons)):
			ax.add_collection(PatchCollection([polygons[i]], color=p_color[i], alpha=p_alpha[i]))
		for i in range(len(p_text)):
			mpl.text(p_text[i][0], p_text[i][1], p_text[i][2], ha='left', fontsize=9)
		mpl.axis(axis_val)
		mpl.yticks([],[])
		mpl.grid(linestyle='--', alpha=0.5)
	mpl.xlabel('read position')
	#
	mpl.tight_layout()
	mpl.savefig(fig_name)
	mpl.close(fig)

#
#	kmer_dat[i] = [[kmer1_hits, kmer2_hits, ...], tlen, tel-anchor-dist, read_orientation, read_name, read_mapq]
#
def plot_kmer_hits(kmer_dat, kmer_colors, my_chr, my_pos, fig_name, clust_dat=None, draw_boundaries=None, plot_params={}):
	which_tel = my_chr[-1]
	max_tlen  = max([n[1] for n in kmer_dat])
	n_reads   = len(kmer_dat)
	#
	stock_params = {'xstep':1000,
	                'xlim':None,
	                'custom_title':None,
	                'custom_xlabel':None,
	                'number_label_rows':True,
	                'fig_width':15}
	for k in plot_params.keys():
		stock_params[k] = plot_params[k]
	#
	X_STEP = stock_params['xstep']
	xlim   = stock_params['xlim']
	#
	mpl.rcParams.update({'font.size': 12, 'font.weight':'normal', 'lines.linewidth':1.0})
	#
	if clust_dat == None:
		read_clusters      = [list(range(n_reads))]
		read_msa_offsets   = [[0]*n_reads]
		total_rows_to_plot = n_reads
		x_axis_adj         = 0
	else:
		read_clusters      = clust_dat[0]
		read_anchor_mapq   = clust_dat[1]
		read_msa_offsets   = clust_dat[2]
		x_axis_adj         = max(clust_dat[3])	# longest subtel length
		#largest_msa_off   = max([max(n) for n in read_msa_offsets])
		total_rows_to_plot = n_reads + len(read_clusters)-1
	#
	# plotting
	#
	vert_fig_size = max(3, total_rows_to_plot * 0.33)
	#
	if which_tel == 'p':
		if xlim != None:
			xtt = [xlim[1]]
			xtl = [-xlim[0]]
			while xtt[-1] > xlim[0]:
				xtt.append(xtt[-1] - X_STEP)
				xtl.append(xtl[-1] - X_STEP)
		else:
			xtt = [max_tlen]
			xtl = [0]
			while xtt[-1] > X_STEP:
				xtt.append(xtt[-1] - X_STEP)
				xtl.append(xtl[-1] - X_STEP)
	else:
		if xlim != None:
			xtt = [n for n in range(xlim[0],xlim[1]+1,X_STEP)]
			xtl = [n for n in range(xlim[0],xlim[1]+1,X_STEP)]
		else:
			xtt = [n for n in range(0,max_tlen,X_STEP)]
			xtl = [n for n in range(0,max_tlen,X_STEP)]
	#
	fig = mpl.figure(1, figsize=(stock_params['fig_width'],vert_fig_size), dpi=200)
	reads_plotted_thus_far = 0
	for clust_i in range(len(read_clusters)):
		for i in range(len(read_clusters[clust_i])):
			my_ind_all = read_clusters[clust_i][i]
			[my_kmer_hits, my_tlen, my_dbta, my_orr, my_rname, my_mapq] = kmer_dat[my_ind_all]
			msa_adj = read_msa_offsets[clust_i][i]
			plot_i  = clust_i + reads_plotted_thus_far
			if plot_i == 0:
				ax1 = mpl.subplot(total_rows_to_plot, 1, plot_i+1)
				if n_reads > 1:
					mpl.setp(ax1.get_xticklabels(), visible=False)
				if which_tel == 'p':
					ax1.yaxis.set_label_position("right")
					ax1.yaxis.tick_right()
				if stock_params['custom_title'] == None:
					if my_pos != None and my_pos != '':
						mpl.title(my_chr + ' : ' + str(my_pos))
					else:
						mpl.title(my_chr)
				else:
					mpl.title(stock_params['custom_title'])
				current_ax = ax1
			else:
				ax2 = mpl.subplot(total_rows_to_plot, 1, plot_i+1, sharex=ax1)
				if plot_i < total_rows_to_plot-1:
					mpl.setp(ax2.get_xticklabels(), visible=False)
				if which_tel == 'p':
					ax2.yaxis.set_label_position("right")
					ax2.yaxis.tick_right()
				current_ax = ax2
			#
			for ki in range(len(my_kmer_hits)):
				if len(my_kmer_hits[ki]):
					if which_tel == 'p':
						if xlim != None:
							adj = xlim[1]  - my_tlen - msa_adj + x_axis_adj + xlim[0]
						else:
							adj = max_tlen - my_tlen - msa_adj + x_axis_adj
					else:
						adj = 0 + msa_adj - x_axis_adj
					polygons = []
					p_color  = []
					p_alpha  = []
					for kmer_span in my_kmer_hits[ki]:
						xp = [kmer_span[0]+adj, kmer_span[1]+adj]
						yp = [-1, 1]
						polygons.append(Polygon(np.array([ [xp[0],yp[0]], [xp[0],yp[1]], [xp[1],yp[1]], [xp[1],yp[0]] ]), closed=True))
						#polygons.append(Polygon(np.array([ [xp[0],yp[0]], [xp[0],yp[1]], [xp[1],(yp[0]+yp[1])/2.0] ]), closed=True))
						p_color.append(kmer_colors[ki])
						p_alpha.append(0.8)
					#
					ax = mpl.gca()
					for j in range(len(polygons)):
						ax.add_collection(PatchCollection([polygons[j]], color=p_color[j], alpha=p_alpha[j], linewidth=0))
			#
			if draw_boundaries != None:
				if which_tel == 'p':
					if xlim != None:
						draw_x = xlim[1] - draw_boundaries[clust_i] + x_axis_adj + xlim[0]
					else:
						draw_x = max_tlen - draw_boundaries[clust_i] + x_axis_adj
				else:
					draw_x = draw_boundaries[clust_i] - x_axis_adj
				#mpl.plot([draw_x, draw_x], [-1,1], '-k', linewidth=3)
				current_ax.add_line(lines.Line2D([draw_x, draw_x], [1.0 ,1.5], color='black', linewidth=3, clip_on=False))

			#
			if stock_params['number_label_rows']:
				mpl.yticks([0],['('+str(my_ind_all)+') '+my_rname])
			else:
				mpl.yticks([0],[my_rname])
			mpl.ylim([-1,1])
			#
			if xlim != None:
				mpl.xlim([xlim[0], xlim[1]])
			else:
				mpl.xlim([0, max_tlen])
			mpl.xticks(xtt, xtl)
			mpl.grid(axis='x', linestyle='--', alpha=0.6)
			#
			reads_plotted_thus_far += 1
	#
	if stock_params['custom_xlabel'] == None:
		mpl.xlabel('distance from subtelomere/telomere boundary (bp)')
	else:
		mpl.xlabel(stock_params['custom_xlabel'])
	mpl.tight_layout()
	mpl.subplots_adjust(hspace=0.0)
	#if my_chr == 'chr2p':
	#	mpl.show()
	#	exit(1)
	mpl.savefig(fig_name)
	mpl.close(fig)

#
#
#
def violin_plotting(dat_p_p, dat_l_p, dat_w_p, dat_p_q, dat_l_q, dat_w_q, plot_params, plot_means=True):
	v_line_keys = ['cmeans', 'cmins', 'cmaxes', 'cbars', 'cmedians', 'cquantiles']
	#
	if len(dat_l_p) and len(dat_p_p):
		violin_parts_p = mpl.violinplot(dat_l_p, dat_p_p, points=200, widths=dat_w_p)
		for pc in violin_parts_p['bodies']:
			pc.set_facecolor(plot_params['p_color'])
			pc.set_edgecolor('black')
			pc.set_alpha(0.7)
		for k in v_line_keys:
			if k in violin_parts_p:
				violin_parts_p[k].set_color('black')
				violin_parts_p[k].set_alpha(0.3)
	#
	if len(dat_l_q) and len(dat_p_q):
		violin_parts_q = mpl.violinplot(dat_l_q, dat_p_q, points=200, widths=dat_w_q)
		for pc in violin_parts_q['bodies']:
			pc.set_facecolor(plot_params['q_color'])
			pc.set_edgecolor('black')
			pc.set_alpha(0.7)
		for k in v_line_keys:
			if k in violin_parts_q:
				violin_parts_q[k].set_color('black')
				violin_parts_q[k].set_alpha(0.3)
	#
	# plot means for each arm
	#
	if plot_means:
		for i in range(len(dat_l_p)):
			yval = np.mean(dat_l_p[i])
			xval = dat_p_p[i]
			mpl.plot([xval - 0.3, xval + 0.3], [yval, yval], '-k', linewidth=2, alpha=0.4)
		for i in range(len(dat_l_q)):
			yval = np.mean(dat_l_q[i])
			xval = dat_p_q[i]
			mpl.plot([xval - 0.3, xval + 0.3], [yval, yval], '-k', linewidth=2, alpha=0.4)

#
#
#
def tel_len_violin_plot(tel_len_dict, out_fn, plot_means=True, ground_truth_dict={}, custom_plot_params={}):
	#
	# plotting constants
	#
	mpl.rcParams.update({'font.size': 18, 'font.weight':'bold'})
	plot_params = {'p_color':'blue',
	               'q_color':'red',
	               'xlabel_rot':0,
	               'y_label':'<-- q   telomere length   p -->',
	               'p_ymax':20000,
	               'q_ymax':20000,
	               'y_step':5000,
	               'fig_size':(16,6),
	               'norm_by_readcount':True,
	               'skip_plot':[],
	               'include_unanchored':False,
	               'boxplot':False,
	               'boxfliers':False,
	               'custom_yticks':None}
	for k in custom_plot_params.keys():
		plot_params[k] = custom_plot_params[k]
	#
	xlab = ['-'] + [str(n) for n in range(1,22+1)] + ['X', 'Y']
	xtck = list(range(1,len(xlab)+1))
	ydel = plot_params['y_step']
	if plot_params['custom_yticks'] != None:
		(p_ymax, q_ymax) = (plot_params['custom_yticks'][0][-1], -plot_params['custom_yticks'][0][0])
		ytck = plot_params['custom_yticks'][0]
		ylab = plot_params['custom_yticks'][1]
	else:
		(p_ymax, q_ymax) = (plot_params['p_ymax'], plot_params['q_ymax'])
		ytck = list(range(-q_ymax, p_ymax+ydel, ydel))
		ylab = []
		for n in ytck:
			if n == 0:
				ylab.append('')
			else:
				ylab.append(str(abs(n)//1000) + 'kb')
	#
	ref_2_x = {'chr'+xlab[i]:xtck[i] for i in range(len(xlab))}
	ref_2_x['unanchored'] = xtck[0]
	ref_2_x['unanchore']  = xtck[0]
	#
	if len(tel_len_dict):
		readcount_denom = max([len(tel_len_dict[k]) for k in tel_len_dict.keys() if k != 'unanchored'])
	else:
		readcount_denom = 1
	width_max = 0.9
	width_min = 0.1
	width_box = 0.6
	#
	# read in lengths and create data structures needed for violin plot
	#
	(dat_l_p, dat_l_q) = ([], [])
	(dat_p_p, dat_p_q) = ([], [])
	(dat_w_p, dat_w_q) = ([], [])
	for k in tel_len_dict.keys():
		if len(tel_len_dict[k]) == 0 or k in plot_params['skip_plot']:
			continue
		if plot_params['norm_by_readcount']:
			my_width = min( [width_max, max([width_min, width_max*(float(len(tel_len_dict[k]))/readcount_denom)])] )
		else:
			if plot_params['boxplot']:
				my_width = width_box
			else:
				my_width = width_max
		if k[-1] == 'p' or k == 'unanchored':
			dat_p_p.append(ref_2_x[k[:-1]])
			dat_l_p.append([])
			dat_w_p.append(my_width)
		elif k[-1] == 'q':
			dat_p_q.append(ref_2_x[k[:-1]])
			dat_l_q.append([])
			dat_w_q.append(my_width)
		for n in tel_len_dict[k]:
			if k[-1] == 'p' or k == 'unanchored':
				dat_l_p[-1].append(n)
			elif k[-1] == 'q':
				dat_l_q[-1].append(-n)
	#
	# violin plot
	#
	fig = mpl.figure(1,figsize=plot_params['fig_size'])
	#
	if plot_params['boxplot']:
		box_params   = {'linewidth':2, 'facecolor':(0.9, 0.9, 0.9)}
		mean_params  = {'linewidth':2, 'linestyle':'solid', 'color':(0.5, 0.5, 0.5)}
		line_params  = {'linewidth':2}
		flier_params = {'marker':'.', 'markerfacecolor':(0.0, 0.0, 0.0), 'markersize':8, 'linestyle':'none', 'alpha':0.2}
		mpl.boxplot(dat_l_p, vert=True, positions=dat_p_p, widths=dat_w_p, patch_artist=True, showfliers=plot_params['boxfliers'], boxprops=box_params, medianprops=mean_params, whiskerprops=line_params, capprops=line_params, flierprops=flier_params)
		mpl.boxplot(dat_l_q, vert=True, positions=dat_p_q, widths=dat_w_q, patch_artist=True, showfliers=plot_params['boxfliers'], boxprops=box_params, medianprops=mean_params, whiskerprops=line_params, capprops=line_params, flierprops=flier_params)
	else:
		violin_plotting(dat_p_p, dat_l_p, dat_w_p, dat_p_q, dat_l_q, dat_w_q, plot_params, plot_means)
	#
	# plot ground-truth values (for simulated data)
	#
	for k in ground_truth_dict.keys():
		xval = ref_2_x[k[:-1]]
		if k[-1] == 'p':
			yval = ground_truth_dict[k]
		elif k[-1] == 'q':
			yval = -ground_truth_dict[k]
		else:
			print('skipping weird benchmark contig:', k, ground_truth_dict[k])
			continue
		mpl.plot([xval - 0.35, xval + 0.35], [yval, yval], '-k', linewidth=3, alpha=1.0)
	#
	mpl.plot([0,len(xlab)+1], [0,0], '-k', linewidth=3)
	if plot_params['include_unanchored']:
		mpl.xticks(xtck, xlab, rotation=plot_params['xlabel_rot'])
		mpl.xlim([0,len(xlab)+1])
	else:
		mpl.xticks(xtck[1:], xlab[1:], rotation=plot_params['xlabel_rot'])
		mpl.xlim([1,len(xlab)+1])
	mpl.yticks(ytck, ylab)
	mpl.ylim([-q_ymax, p_ymax])
	mpl.ylabel(plot_params['y_label'])
	mpl.grid(linestyle='--', alpha=0.5)
	mpl.tight_layout()
	mpl.savefig(out_fn)
	mpl.close(fig)

#
# tel_len_by_samp[i] = (samp_name, 'p'/'q', tlen_list)
#
# ground_truth_by_samp[i] = (samp_name, 'p'/'q', tlen)
#
def tel_len_violin_single_chr_multiple_samples(tel_len_by_samp, out_fn, plot_means=True, ground_truth_by_samp=[], custom_plot_params={}):
	#
	# plotting constants
	#
	mpl.rcParams.update({'font.size': 18, 'font.weight':'bold'})
	plot_params = {'p_color':'blue',
	               'q_color':'red',
	               'xlabel_rot':0,
	               'y_label':'<-- q   telomere length   p -->',
	               'p_ymax':20000,
	               'q_ymax':20000,
	               'y_step':5000,
	               'fig_size':(16,6)}
	for k in custom_plot_params.keys():
		plot_params[k] = custom_plot_params[k]
	#
	xlab = []
	for n in tel_len_by_samp:
		if n[0] not in xlab:
			xlab.append(n[0])
	xtck = list(range(1,len(xlab)+1))
	ydel = plot_params['y_step']
	(p_ymax, q_ymax) = (plot_params['p_ymax'], plot_params['q_ymax'])
	ytck = list(range(-q_ymax, p_ymax+ydel, ydel))
	ylab = []
	for n in ytck:
		if n == 0:
			ylab.append('')
		else:
			ylab.append(str(abs(n)//1000) + 'k')
	#
	samp_2_x = {xlab[i]:xtck[i] for i in range(len(xlab))}
	#
	if len(tel_len_by_samp):
		readcount_denom = max([len(n[2]) for n in tel_len_by_samp])
	else:
		readcount_denom = 1
	width_max = 1.0
	width_min = 0.1
	#
	# read in lengths and create data structures needed for violin plot
	#
	(dat_l_p, dat_l_q) = ([], [])
	(dat_p_p, dat_p_q) = ([], [])
	(dat_w_p, dat_w_q) = ([], [])
	for (samp_name, pq, tlen_list) in tel_len_by_samp:
		if len(tlen_list) == 0:
			continue
		my_width = min( [width_max, max([width_min, width_max*(float(len(tlen_list))/readcount_denom)])] )
		if pq == 'p':
			dat_p_p.append(samp_2_x[samp_name])
			dat_l_p.append([])
			dat_w_p.append(my_width)
		elif pq == 'q':
			dat_p_q.append(samp_2_x[samp_name])
			dat_l_q.append([])
			dat_w_q.append(my_width)
		for n in tlen_list:
			if pq == 'p':
				dat_l_p[-1].append(n)
			elif pq == 'q':
				dat_l_q[-1].append(-n)
	#
	# violin plot
	#
	fig = mpl.figure(1,figsize=plot_params['fig_size'])
	#
	violin_plotting(dat_p_p, dat_l_p, dat_w_p, dat_p_q, dat_l_q, dat_w_q, plot_params, plot_means)
	#
	# use the ground-truth plotting function to compare against other tl methods (e.g. denovo assembly)
	#
	for i in range(len(ground_truth_by_samp)):
		xval = samp_2_x[ground_truth_by_samp[i][0]]
		if ground_truth_by_samp[i][1] == 'p':
			yval = ground_truth_by_samp[i][2]
		elif ground_truth_by_samp[i][1] == 'q':
			yval = -ground_truth_by_samp[i][2]
		else:
			continue
		mpl.plot([xval - 0.35, xval + 0.35], [yval, yval], '-k', linewidth=2, alpha=1.0)
	#
	mpl.plot([0,len(xlab)+1], [0,0], '-k', linewidth=3)
	mpl.xticks(xtck, xlab, rotation=plot_params['xlabel_rot'])
	mpl.xlim([0,len(xlab)+1])
	mpl.yticks(ytck, ylab)
	mpl.ylim([-q_ymax, p_ymax])
	mpl.ylabel(plot_params['y_label'])
	mpl.grid(linestyle='--', alpha=0.5)
	mpl.tight_layout()
	mpl.savefig(out_fn)
	mpl.close(fig)

#
#
#
def tel_len_bar_plot(tel_len_dict, out_fn, custom_plot_params={}):
	#
	# plotting constants
	#
	mpl.rcParams.update({'font.size': 18, 'font.weight':'bold'})
	plot_params = {'p_color':'blue',
	               'q_color':'red',
	               'xlabel_rot':0,
	               'y_label':'<-- q   telomere length   p -->',
	               'p_ymax':20000,
	               'q_ymax':20000,
	               'y_step':5000,
	               'fig_size':(16,6),
	               'skip_plot':[],
	               'include_unanchored':False,
	               'ytick_suffix':'',
	               'hatch_data':None}
	for k in custom_plot_params.keys():
		plot_params[k] = custom_plot_params[k]
	#
	xlab = ['-'] + [str(n) for n in range(1,22+1)] + ['X', 'Y']
	xtck = list(range(1,len(xlab)+1))
	ydel = plot_params['y_step']
	(p_ymax, q_ymax) = (plot_params['p_ymax'], plot_params['q_ymax'])
	ytck = list(range(-q_ymax, p_ymax+ydel, ydel))
	ylab = [str(abs(n))+plot_params['ytick_suffix'] for n in ytck]
	#
	ref_2_x = {'chr'+xlab[i]:xtck[i] for i in range(len(xlab))}
	ref_2_x['unanchored'] = xtck[0]
	ref_2_x['unanchore']  = xtck[0]
	#
	# read in lengths and create data structures needed for violin plot
	#
	(dat_l_p, dat_l_q) = ([], [])
	(dat_p_p, dat_p_q) = ([], [])
	for k in tel_len_dict.keys():
		if k in plot_params['skip_plot']:
			continue
		if k[-1] == 'p' or k == 'unanchored':
			dat_p_p.append(ref_2_x[k[:-1]])
			dat_l_p.append(tel_len_dict[k])
		elif k[-1] == 'q':
			dat_p_q.append(ref_2_x[k[:-1]])
			dat_l_q.append(tel_len_dict[k])
	#
	if plot_params['hatch_data'] != None:
		(hat_l_p, hat_l_q) = ([], [])
		(hat_p_p, hat_p_q) = ([], [])
		for k in plot_params['hatch_data'].keys():
			if k in plot_params['skip_plot']:
				continue
			if k[-1] == 'p' or k == 'unanchored':
				hat_p_p.append(ref_2_x[k[:-1]])
				hat_l_p.append(plot_params['hatch_data'][k])
			elif k[-1] == 'q':
				hat_p_q.append(ref_2_x[k[:-1]])
				hat_l_q.append(plot_params['hatch_data'][k])
	#
	# violin plot
	#
	fig = mpl.figure(1,figsize=plot_params['fig_size'])
	#
	polygons = []
	p_color  = []
	p_alpha  = []
	for i in range(len(dat_p_p)):
		xp = [dat_p_p[i]-0.3, dat_p_p[i]+0.3]
		yp = [0, dat_l_p[i]]
		polygons.append(Polygon(np.array([ [xp[0],yp[0]], [xp[0],yp[1]], [xp[1],yp[1]], [xp[1],yp[0]] ]), closed=True))
		p_color.append((0.6, 0.6, 0.6))
		p_alpha.append(1.0)
	for i in range(len(dat_p_q)):
		xp = [dat_p_q[i]-0.3, dat_p_q[i]+0.3]
		yp = [0, -dat_l_q[i]]
		polygons.append(Polygon(np.array([ [xp[0],yp[0]], [xp[0],yp[1]], [xp[1],yp[1]], [xp[1],yp[0]] ]), closed=True))
		p_color.append((0.6, 0.6, 0.6))
		p_alpha.append(1.0)
	#
	polygons2 = []
	p_color2  = []
	if plot_params['hatch_data'] != None:
		for i in range(len(hat_p_p)):
			xp = [hat_p_p[i]-0.3, hat_p_p[i]+0.3]
			yp = [0, hat_l_p[i]]
			polygons2.append(Polygon(np.array([ [xp[0],yp[0]], [xp[0],yp[1]], [xp[1],yp[1]], [xp[1],yp[0]] ]), closed=True))
			p_color2.append((0.45, 0.45, 0.45))
		for i in range(len(hat_p_p)):
			xp = [hat_p_q[i]-0.3, hat_p_q[i]+0.3]
			yp = [0, -hat_l_q[i]]
			polygons2.append(Polygon(np.array([ [xp[0],yp[0]], [xp[0],yp[1]], [xp[1],yp[1]], [xp[1],yp[0]] ]), closed=True))
			p_color2.append((0.45, 0.45, 0.45))
	#
	ax = mpl.gca()
	for j in range(len(polygons)):
		ax.add_collection(PatchCollection([polygons[j]], color=p_color[j], alpha=p_alpha[j], linewidth=0))
	for j in range(len(polygons2)):
		ax.add_collection(PatchCollection([polygons2[j]], facecolor=p_color2[j], hatch='//', linewidth=0))
	#
	mpl.plot([0,len(xlab)+1], [0,0], '-k', linewidth=3)
	if plot_params['include_unanchored']:
		mpl.xticks(xtck, xlab, rotation=plot_params['xlabel_rot'])
		mpl.xlim([0,len(xlab)+1])
	else:
		mpl.xticks(xtck[1:], xlab[1:], rotation=plot_params['xlabel_rot'])
		mpl.xlim([1,len(xlab)+1])
	mpl.yticks(ytck, ylab)
	mpl.ylim([-q_ymax, p_ymax])
	mpl.ylabel(plot_params['y_label'])
	mpl.grid(linestyle='--', alpha=0.5)
	mpl.tight_layout()
	mpl.savefig(out_fn)
	mpl.close(fig)
