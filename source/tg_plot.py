import numpy as np
import matplotlib.pyplot as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.lines as lines
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from source.tg_align import UNKNOWN
from source.tg_util  import exists_and_is_nonzero

MAX_PLOT_SIZE = 65000 # the actual max is 65535 but it's risky
DEFAULT_DPI   = 100
KMER_HITS_DPI = 200


def getColor(i, N, colormap='jet'):
    cm = mpl.get_cmap(colormap)
    cNorm  = colors.Normalize(vmin=0, vmax=N+1)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    colorVal = scalarMap.to_rgba(i)
    return colorVal


def get_read_alignment_polygons(alignments, readlen):
    polygons = []
    p_alpha  = []
    p_color  = []
    p_text   = []
    yp       = [-1.0, 1.0]
    text_y   = 0.7 * yp[1]
    for i in range(len(alignments)):
        n = alignments[i]
        #
        my_refname = n[2]
        my_refspan = str(n[3])+' - '+str(n[4])
        if my_refname[:3] == 'tel':
            if my_refname[-1] == '?':   # skip these weird regions
                p_color.append('gray')
            else:
                my_refname = 'tel'
                p_color.append('red')
        else:
            p_color.append('blue')
        p_alpha.append(float(n[6]+15.)/(60.+15.))
        #
        xp               = [n[0], n[1]]
        delta_pointy     = 0.75*(n[1]-n[0])
        delta_pointy_rev = 0.25*(n[1]-n[0])
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


def plot_tel_signal(density_data, f_title, fig_name, tl_vals=None, interstitial_tels=[], plot_for_paper=False):
    #
    [td_p_e0, td_p_e1, td_q_e0, td_q_e1, p_vs_q_power, rlen, tel_window] = density_data
    #
    x_off   = tel_window / 2
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
    #mpl.rcParams.update({'text.usetex': True})
    #mpl.legend(['$p_1(i)$', '$p_0(i)$'], loc=1)
    #mpl.rcParams.update({'text.usetex': False})
    mpl.legend(['p1(i)', 'p0(i)'], loc=1)
    mpl.xlim([0,rlen])
    mpl.yticks(dens_yt, dens_yl)
    mpl.ylim([0,1])
    mpl.grid(linestyle='--', alpha=0.5)
    mpl.ylabel('tel density (p)')
    if plot_for_paper is False:
        mpl.title(f_title)
    #
    mpl.subplot(subplot_nums[1], sharex=ax1)
    mpl.plot(np.arange(len(td_q_e1))+x_off, td_q_e1, '-m', alpha=0.5)
    mpl.plot(np.arange(len(td_q_e0))+x_off, td_q_e0, '-k', alpha=0.8)
    #mpl.rcParams.update({'text.usetex': True})
    #mpl.legend(['$q_1(i)$', '$q_0(i)$'], loc=1)
    #mpl.rcParams.update({'text.usetex': False})
    mpl.legend(['q1(i)', 'q0(i)'], loc=1)
    mpl.yticks(dens_yt, dens_yl)
    mpl.ylim([0,1])
    mpl.grid(linestyle='--', alpha=0.5)
    mpl.ylabel('tel density (q)')
    #
    mpl.subplot(subplot_nums[2], sharex=ax1)
    mpl.plot(np.arange(len(p_vs_q_power))+x_off, p_vs_q_power, '-k')
    #mpl.rcParams.update({'text.usetex': True})
    #mpl.legend(['$S(i)$'], loc=1)
    #mpl.rcParams.update({'text.usetex': False})
    mpl.legend(['S(i)'], loc=1)
    mpl.grid(linestyle='--', alpha=0.5)
    mpl.ylabel('tel score')
    mpl.yticks(pred_yt, pred_yl)
    mpl.ylim([-1.0, 1.0])
    #
    polygons = []
    p_color  = []
    p_alpha  = []
    yp = [-1, 1]
    if tl_vals is not None:
        [tl_p, tl_q] = tl_vals
        if tl_p > 0:
            xp = [0, tl_p]
            polygons.append(Polygon(np.array([[xp[0],yp[0]], [xp[0],yp[1]], [xp[1],yp[1]], [xp[1],yp[0]]]), closed=True))
            p_color.append('red')
            p_alpha.append(0.5)
        if tl_q > 0:
            xp = [rlen - tl_q, rlen]
            polygons.append(Polygon(np.array([[xp[0],yp[0]], [xp[0],yp[1]], [xp[1],yp[1]], [xp[1],yp[0]]]), closed=True))
            p_color.append('red')
            p_alpha.append(0.5)
    for xp in interstitial_tels:
        if xp is not None:
            polygons.append(Polygon(np.array([[xp[0],yp[0]], [xp[0],yp[1]], [xp[1],yp[1]], [xp[1],yp[0]]]), closed=True))
            p_color.append('purple')
            p_alpha.append(0.5)
    if len(polygons):
        ax = mpl.gca()
        for i in range(len(polygons)):
            ax.add_collection(PatchCollection([polygons[i]], color=p_color[i], alpha=p_alpha[i]))
    #
    mpl.xlabel('read position')
    mpl.tight_layout()
    mpl.savefig(fig_name)
    mpl.close(fig)


def plot_kmer_hits(kmer_dat, kmer_colors, my_chr, my_pos, fig_name, clust_dat=None, draw_boundaries=None, plot_params={}):
    #
    #   kmer_dat[i] = [[kmer1_hits, kmer2_hits, ...], tlen, tel-anchor-dist, read_orientation, read_name, read_mapq, read_fasta]
    #
    if my_chr:
        which_tel = my_chr[-1]
    else:
        which_tel = 'q'
    max_tlen  = max([n[1] for n in kmer_dat])
    n_reads   = len(kmer_dat)
    #
    stock_params = {'xstep':1000,
                    'xlim':None,
                    'custom_title':None,
                    'custom_xlabel':None,
                    'number_label_rows':True,
                    'fig_width':15,
                    'fig_height':None,
                    'font.size':12,
                    'font.weight':'normal',
                    'draw_grid':True,
                    'mpl_fignum':1,
                    'dpi':KMER_HITS_DPI}
    for k in plot_params.keys():
        stock_params[k] = plot_params[k]
    #
    X_STEP = stock_params['xstep']
    xlim   = stock_params['xlim']
    #
    mpl.rcParams.update({'font.size':stock_params['font.size'],
                         'font.weight':stock_params['font.weight'],
                         'lines.linewidth':1.0})
    #
    if clust_dat is None:
        read_clusters      = [list(range(n_reads))]
        read_msa_offsets   = [[0]*n_reads]
        total_rows_to_plot = n_reads
        x_axis_adj         = 0
    else:
        read_clusters      = clust_dat[0]
        read_msa_offsets   = clust_dat[2]
        x_axis_adj         = max(clust_dat[3])  # longest subtel length
        total_rows_to_plot = n_reads + len(read_clusters)-1
    #
    # plotting
    #
    vert_fig_size = max(3, total_rows_to_plot * 0.35)
    vert_fig_size = min(vert_fig_size, MAX_PLOT_SIZE / stock_params['dpi'])
    if stock_params['fig_height'] is not None:
        vert_fig_size = stock_params['fig_height']
    #
    if which_tel == 'p':
        if xlim is not None:
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
        if xlim is not None:
            xtt = [n for n in range(xlim[0],xlim[1]+1,X_STEP)]
            xtl = [n for n in range(xlim[0],xlim[1]+1,X_STEP)]
        else:
            xtt = [n for n in range(0,max_tlen,X_STEP)]
            xtl = [n for n in range(0,max_tlen,X_STEP)]
    #
    fig = mpl.figure(stock_params['mpl_fignum'], figsize=(stock_params['fig_width'],vert_fig_size), dpi=stock_params['dpi'], layout="constrained")
    fig.get_layout_engine().set(h_pad=0.0, hspace=0.0)
    #
    reads_plotted_thus_far = 0
    for clust_i in range(len(read_clusters)):
        for i in range(len(read_clusters[clust_i])):
            my_ind_all = read_clusters[clust_i][i]
            [my_kmer_hits, my_tlen, my_dbta, my_orr, my_rname, my_mapq, my_fastadat] = kmer_dat[my_ind_all]
            msa_adj = read_msa_offsets[clust_i][i]
            plot_i  = clust_i + reads_plotted_thus_far
            if plot_i == 0:
                ax1 = mpl.subplot(total_rows_to_plot, 1, plot_i+1)
                if n_reads > 1:
                    for tick in ax1.xaxis.get_major_ticks():
                        tick.tick1line.set_visible(False)
                        tick.tick2line.set_visible(False)
                        tick.label1.set_visible(False)
                        tick.label2.set_visible(False)
                if which_tel == 'p':
                    ax1.yaxis.set_label_position("right")
                    ax1.yaxis.tick_right()
                if stock_params['custom_title'] is None:
                    if my_pos is not None and my_pos != '':
                        mpl.title(my_chr + ' : ' + str(my_pos))
                    else:
                        mpl.title(my_chr)
                else:
                    mpl.title(stock_params['custom_title'])
                current_ax = ax1
            else:
                ax2 = mpl.subplot(total_rows_to_plot, 1, plot_i+1, sharex=ax1)
                if plot_i < total_rows_to_plot-1:
                    for tick in ax2.xaxis.get_major_ticks():
                        tick.tick1line.set_visible(False)
                        tick.tick2line.set_visible(False)
                        tick.label1.set_visible(False)
                        tick.label2.set_visible(False)
                if which_tel == 'p':
                    ax2.yaxis.set_label_position("right")
                    ax2.yaxis.tick_right()
                current_ax = ax2
            #
            my_adj_params = [my_tlen, max_tlen, msa_adj, x_axis_adj]
            plot_kmer_hits_on_current_plot(my_kmer_hits, kmer_colors, my_adj_params, which_tel=which_tel, xlim=xlim)
            #
            if draw_boundaries is not None:
                if which_tel == 'p':
                    if xlim is not None:
                        draw_x = xlim[1] - draw_boundaries[clust_i] + x_axis_adj + xlim[0]
                    else:
                        draw_x = max_tlen - draw_boundaries[clust_i] + x_axis_adj
                else:
                    draw_x = draw_boundaries[clust_i] - x_axis_adj
                current_ax.add_line(lines.Line2D([draw_x, draw_x], [1.0,1.5], color='black', linewidth=3, clip_on=False))

            #
            if stock_params['number_label_rows']:
                mpl.yticks([0],['('+str(my_ind_all)+') '+my_rname])
            else:
                mpl.yticks([0],[my_rname])
            mpl.ylim([-1,1])
            #
            if xlim is not None:
                mpl.xlim([xlim[0], xlim[1]])
            else:
                mpl.xlim([0, max_tlen])
            mpl.xticks(xtt, xtl)
            if stock_params['draw_grid']:
                mpl.grid(axis='x', linestyle='--', alpha=0.6)
            #
            reads_plotted_thus_far += 1
    #
    if stock_params['custom_xlabel'] is None:
        mpl.xlabel('distance from subtelomere boundary (bp)')
    else:
        mpl.xlabel(stock_params['custom_xlabel'])
    #
    mpl.savefig(fig_name)
    mpl.close(fig)


def plot_kmer_hits_on_current_plot(kmer_hits, kmer_colors, adj_params, which_tel='q', xlim=None):
    [my_tlen, max_tlen, msa_adj, x_axis_adj] = adj_params
    for ki in range(len(kmer_hits)):
        if len(kmer_hits[ki]):
            if which_tel == 'p':
                if xlim is not None:
                    adj = xlim[1]  - my_tlen - msa_adj + x_axis_adj + xlim[0]
                else:
                    adj = max_tlen - my_tlen - msa_adj + x_axis_adj
            else:
                adj = 0 + msa_adj - x_axis_adj
            polygons = []
            p_color  = []
            p_alpha  = []
            for kmer_span in kmer_hits[ki]:
                xp = [kmer_span[0]+adj, kmer_span[1]+adj]
                yp = [-1, 1]
                polygons.append(Polygon(np.array([[xp[0],yp[0]], [xp[0],yp[1]], [xp[1],yp[1]], [xp[1],yp[0]]]), closed=True))
                #polygons.append(Polygon(np.array([[xp[0],yp[0]], [xp[0],yp[1]], [xp[1],(yp[0]+yp[1])/2.0]]), closed=True))
                p_color.append(kmer_colors[ki])
                p_alpha.append(0.8)
            #
            ax = mpl.gca()
            for j in range(len(polygons)):
                ax.add_collection(PatchCollection([polygons[j]], color=p_color[j], alpha=p_alpha[j], linewidth=0))


def plot_fusion(readname, anchordat1, anchordat2, kmer_dat_tuple, kmer_metadata, fig_name):
    (kmer_dat_fwd, kmer_dat_rev) = kmer_dat_tuple
    [my_kmer_hits_fwd, my_tlen, my_dbta, my_orr, my_rname, my_mapq, my_fastadat] = kmer_dat_fwd
    [my_kmer_hits_rev, my_tlen, my_dbta, my_orr, my_rname, my_mapq, my_fastadat] = kmer_dat_rev
    [KMER_LIST, KMER_COLORS, KMER_LETTER, KMER_FLAGS] = kmer_metadata
    my_adj_params = [my_tlen, my_tlen, 0, 0]
    #
    mpl.rcParams.update({'font.size':12,
                         'font.weight':'normal',
                         'lines.linewidth':1.0})
    #
    polygons = []
    p_color  = []
    p_alpha  = []
    p_text   = []
    [rspan_start, rspan_end, refname, refspan_start, refspan_end, orientation, mapq] = anchordat1[6]
    my_chr = refname.split('_')[-1]
    xp = [rspan_start, rspan_end]
    yp = [-1, 1]
    if refspan_end > refspan_start:
        polygons.append(Polygon(np.array([[xp[0],yp[0]], [xp[0],yp[1]], [xp[1],(yp[0]+yp[1])/2.0]]), closed=True))
    else:
        polygons.append(Polygon(np.array([[xp[1],yp[0]], [xp[1],yp[1]], [xp[0],(yp[0]+yp[1])/2.0]]), closed=True))
    p_color.append('purple')
    p_alpha.append(float(mapq + 15) / (60 + 15))
    p_text.append((xp[0], -0.05, f'{my_chr}:{refspan_start}-{refspan_end}'))
    #
    [rspan_start, rspan_end, refname, refspan_start, refspan_end, orientation, mapq] = anchordat2[6]
    my_chr = refname.split('_')[-1]
    offset_rightread = my_tlen - anchordat2[2][0]
    xp = [rspan_start + offset_rightread, rspan_end + offset_rightread]
    yp = [-1, 1]
    if refspan_end > refspan_start:
        polygons.append(Polygon(np.array([[xp[0],yp[0]], [xp[0],yp[1]], [xp[1],(yp[0]+yp[1])/2.0]]), closed=True))
    else:
        polygons.append(Polygon(np.array([[xp[1],yp[0]], [xp[1],yp[1]], [xp[0],(yp[0]+yp[1])/2.0]]), closed=True))
    p_color.append('purple')
    p_alpha.append(float(mapq + 15) / (60 + 15))
    p_text.append((xp[0], -0.05, f'{my_chr}:{refspan_start}-{refspan_end}'))
    #
    fig = mpl.figure(1, figsize=(16,4), dpi=200)
    ax1 = mpl.subplot(311)
    for i in range(len(polygons)):
        ax1.add_collection(PatchCollection([polygons[i]], color=p_color[i], alpha=p_alpha[i], linewidth=0))
    for i in range(len(p_text)):
        mpl.text(p_text[i][0], p_text[i][1], p_text[i][2], ha='left', fontsize=9)
    mpl.title(readname)
    mpl.xticks([],[])
    mpl.yticks([],[])
    mpl.xlim([0, my_tlen])
    mpl.ylim([-1,1])
    #
    mpl.subplot(312)
    plot_kmer_hits_on_current_plot(my_kmer_hits_fwd, KMER_COLORS, my_adj_params)
    mpl.xticks([],[])
    mpl.yticks([],[])
    mpl.xlim([0, my_tlen])
    mpl.ylim([-1,1])
    #
    mpl.subplot(313)
    plot_kmer_hits_on_current_plot(my_kmer_hits_rev, KMER_COLORS, my_adj_params)
    mpl.xlim([0, my_tlen])
    mpl.yticks([],[])
    mpl.ylim([-1,1])
    mpl.xlabel('read coordinate')
    #
    mpl.tight_layout()
    mpl.savefig(fig_name)
    mpl.close(fig)


def violin_plotting(dat_p_p, dat_l_p, dat_w_p, dat_p_q, dat_l_q, dat_w_q, plot_params, plot_means=True):
    v_line_keys = ['cmeans', 'cmins', 'cmaxes', 'cbars', 'cmedians', 'cquantiles']
    #
    if len(dat_l_p) and len(dat_p_p):
        violin_parts_p = mpl.violinplot(dat_l_p, dat_p_p, points=200, widths=dat_w_p)
        for pc in violin_parts_p['bodies']:
            pc.set_facecolor(plot_params['p_color'])
            pc.set_edgecolor('black')
            pc.set_alpha(0.85)
        for k in v_line_keys:
            if k in violin_parts_p:
                violin_parts_p[k].set_color('black')
                violin_parts_p[k].set_alpha(0.4)
    #
    if len(dat_l_q) and len(dat_p_q):
        violin_parts_q = mpl.violinplot(dat_l_q, dat_p_q, points=200, widths=dat_w_q)
        for pc in violin_parts_q['bodies']:
            pc.set_facecolor(plot_params['q_color'])
            pc.set_edgecolor('black')
            pc.set_alpha(0.85)
        for k in v_line_keys:
            if k in violin_parts_q:
                violin_parts_q[k].set_color('black')
                violin_parts_q[k].set_alpha(0.4)
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


def tel_len_violin_plot(tel_len_dict_list, out_fn, plot_means=True, custom_plot_params={}):
    #
    mpl.rcParams.update({'font.size': 18, 'font.weight':'normal'})
    plot_params = {'p_colors':[(70/255, 130/255, 180/255)], # steelblue
                   'q_colors':[(186/255, 85/255, 211/255)], # mediumorchid
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
                   'custom_yticks':None,
                   'custom_title':'',
                   'spacer_between_alleles':2,
                   'legend':[]}
    for k in custom_plot_params.keys():
        if k in plot_params:
            plot_params[k] = custom_plot_params[k]
        else:
            print('Error: tel_len_violin_plot() was given an invalid plot_param:', k)
            exit(1)
    #
    alleles_per_arm = len(tel_len_dict_list)
    spacer_between_alleles = plot_params['spacer_between_alleles']
    samp_names = []
    if plot_params['legend']:
        samp_names = plot_params['legend'][:min(alleles_per_arm, len(plot_params['legend']))]
    #
    xlab_temp = [str(n) for n in range(1,22+1)] + ['X', 'Y']
    xlab = ['-']
    for n in xlab_temp:
        xlab.append(n)
        if alleles_per_arm > 1:
            for i in range(1,alleles_per_arm+spacer_between_alleles):
                xlab.append('')
    xtck = list(range(1,len(xlab)+1))
    ydel = plot_params['y_step']
    if plot_params['custom_yticks'] is not None:
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
    ref_2_x['unanchore'] = xtck[0]
    if plot_params['include_unanchored'] is False:
        plot_params['skip_plot'].append('unanchored')
    #
    width_max = 0.85
    width_min = 0.10
    width_box = 0.70
    #
    # read in lengths and create data structures needed for violin plot, then do the plotting
    #
    fig = mpl.figure(1, figsize=plot_params['fig_size'], dpi=200) # default dpi = 100
    # plot some line off screen for legend entries
    for tldi in range(len(tel_len_dict_list)):
        mpl.plot([-2,-2], [0,1], color=plot_params['p_colors'][min(tldi, len(plot_params['p_colors'])-1)])
    #
    for tldi,tel_len_dict in enumerate(tel_len_dict_list):
        dat_l_p, dat_l_q = [], []
        dat_p_p, dat_p_q = [], []
        dat_w_p, dat_w_q = [], []
        #
        if len(tel_len_dict):
            tel_readcounts = [len(tel_len_dict[k]) for k in tel_len_dict.keys() if k != 'unanchored']
            if len(tel_readcounts):
                readcount_denom = max(tel_readcounts)
            else:
                readcount_denom = 1
        else:
            readcount_denom = 1
        #
        for k in tel_len_dict.keys():
            if len(tel_len_dict[k]) == 0 or k in plot_params['skip_plot']:
                continue
            if plot_params['norm_by_readcount']:
                my_width = min([width_max, max([width_min, width_max*(float(len(tel_len_dict[k]))/readcount_denom)])])
            else:
                if plot_params['boxplot']:
                    my_width = width_box
                else:
                    my_width = width_max
            if k[-1] == 'p' or k == 'unanchored':
                dat_p_p.append(ref_2_x[k[:-1]] + tldi)
                dat_l_p.append([])
                dat_w_p.append(my_width)
            elif k[-1] == 'q':
                dat_p_q.append(ref_2_x[k[:-1]] + tldi)
                dat_l_q.append([])
                dat_w_q.append(my_width)
            for n in tel_len_dict[k]:
                if k[-1] == 'p' or k == 'unanchored':
                    dat_l_p[-1].append(n)
                elif k[-1] == 'q':
                    dat_l_q[-1].append(-n)
        #
        plot_params['p_color'] = plot_params['p_colors'][min(tldi, len(plot_params['p_colors'])-1)]
        plot_params['q_color'] = plot_params['q_colors'][min(tldi, len(plot_params['q_colors'])-1)]
        #
        if plot_params['boxplot']:
            mean_params  = {'linewidth':1, 'linestyle':'dotted', 'color':(0.1, 0.1, 0.1)}
            line_params  = {'linewidth':1}
            flier_params = {'marker':'.', 'markerfacecolor':(0.0, 0.0, 0.0), 'markersize':8, 'linestyle':'none', 'alpha':0.2}
            box_params = {'linewidth':1, 'facecolor':plot_params['p_color']}
            mpl.boxplot(dat_l_p, vert=True, positions=dat_p_p, widths=dat_w_p, patch_artist=True, showfliers=plot_params['boxfliers'], boxprops=box_params, medianprops=mean_params, whiskerprops=line_params, capprops=line_params, flierprops=flier_params)
            box_params = {'linewidth':1, 'facecolor':plot_params['q_color']}
            mpl.boxplot(dat_l_q, vert=True, positions=dat_p_q, widths=dat_w_q, patch_artist=True, showfliers=plot_params['boxfliers'], boxprops=box_params, medianprops=mean_params, whiskerprops=line_params, capprops=line_params, flierprops=flier_params)
        else:
            violin_plotting(dat_p_p, dat_l_p, dat_w_p, dat_p_q, dat_l_q, dat_w_q, plot_params, plot_means)
    #
    xt_packed = [(xtck[i]+0.5*(alleles_per_arm-1), xlab[i]) for i in range(len(xtck)) if len(xlab[i])]
    num_blank = len([n for n in xlab if len(n) == 0])
    if alleles_per_arm > 1:
        num_blank -= spacer_between_alleles  # don't need right-padding for chrY because it's the last
    xtck = [n[0] for n in xt_packed]
    xlab = [n[1] for n in xt_packed]
    mpl.plot([0,len(xlab)+num_blank+1], [0,0], '-k', linewidth=3)
    if 'chrYp' in plot_params['skip_plot'] and 'chrYq' in plot_params['skip_plot']:
        xtck = xtck[:-1]
        xlab = xlab[:-1]
    if plot_params['include_unanchored']:
        mpl.xticks(xtck, xlab, rotation=plot_params['xlabel_rot'], fontweight='bold')
        mpl.xlim([0,len(xlab)+num_blank+1])
    else:
        mpl.xticks(xtck[1:], xlab[1:], rotation=plot_params['xlabel_rot'], fontweight='bold')
        mpl.xlim([1,len(xlab)+num_blank+1])
    mpl.yticks(ytck, ylab)
    mpl.ylim([-q_ymax, p_ymax])
    mpl.ylabel(plot_params['y_label'], fontweight='bold')
    if plot_params['custom_title']:
        mpl.title(plot_params['custom_title'])
    if samp_names:
        mpl.legend(samp_names, prop={'size':12})
    mpl.grid(linestyle='--', alpha=0.5)
    mpl.tight_layout()
    mpl.savefig(out_fn)
    mpl.close(fig)


def tel_len_violin_single_chr_multiple_samples(tel_len_by_samp, out_fn, plot_means=True, ground_truth_by_samp=[], custom_plot_params={}):
    #
    # tel_len_by_samp[i] = (samp_name, 'p'/'q', tlen_list)
    #
    # ground_truth_by_samp[i] = (samp_name, 'p'/'q', tlen)
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
        my_width = min([width_max, max([width_min, width_max*(float(len(tlen_list))/readcount_denom)])])
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


def tel_len_bar_plot(tel_len_dict, out_fn, custom_plot_params={}):
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
    if plot_params['hatch_data'] is not None:
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
        polygons.append(Polygon(np.array([[xp[0],yp[0]], [xp[0],yp[1]], [xp[1],yp[1]], [xp[1],yp[0]]]), closed=True))
        p_color.append((0.6, 0.6, 0.6))
        p_alpha.append(1.0)
    for i in range(len(dat_p_q)):
        xp = [dat_p_q[i]-0.3, dat_p_q[i]+0.3]
        yp = [0, -dat_l_q[i]]
        polygons.append(Polygon(np.array([[xp[0],yp[0]], [xp[0],yp[1]], [xp[1],yp[1]], [xp[1],yp[0]]]), closed=True))
        p_color.append((0.6, 0.6, 0.6))
        p_alpha.append(1.0)
    #
    polygons2 = []
    p_color2  = []
    if plot_params['hatch_data'] is not None:
        for i in range(len(hat_p_p)):
            xp = [hat_p_p[i]-0.3, hat_p_p[i]+0.3]
            yp = [0, hat_l_p[i]]
            polygons2.append(Polygon(np.array([[xp[0],yp[0]], [xp[0],yp[1]], [xp[1],yp[1]], [xp[1],yp[0]]]), closed=True))
            p_color2.append((0.45, 0.45, 0.45))
        for i in range(len(hat_p_q)):
            xp = [hat_p_q[i]-0.3, hat_p_q[i]+0.3]
            yp = [0, -hat_l_q[i]]
            polygons2.append(Polygon(np.array([[xp[0],yp[0]], [xp[0],yp[1]], [xp[1],yp[1]], [xp[1],yp[0]]]), closed=True))
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


def make_tvr_plots(kmer_hit_dat, read_clust_dat, my_chr, my_pos, telcompplot_fn, telcompcons_fn, mtp_params, dpi=None):
    #
    [KMER_METADATA, KMER_COLORS, MIN_READS_PER_PHASE, PLOT_FILT_CVECS, DUMMY_TEL_MAPQ, DO_NOT_OVERWRITE] = mtp_params
    #
    # adjust kmer_hit_dat based on the filters and etc that were applied during clustering
    #
    my_consensus_vecs = read_clust_dat[4]
    my_color_vectors  = read_clust_dat[5]
    my_end_err_lens   = read_clust_dat[6]
    my_tvr_tel_bounds = read_clust_dat[7]
    redrawn_consensus = convert_colorvec_to_kmerhits(my_consensus_vecs, KMER_METADATA)
    # do we want to plot denoised tvrs of individual reads?
    if PLOT_FILT_CVECS:
        redrawn_kmerhits = convert_colorvec_to_kmerhits(my_color_vectors, KMER_METADATA)
        for rdki in range(len(redrawn_kmerhits)):
            kmer_hit_dat[rdki][0]  = redrawn_kmerhits[rdki] # replace kmer hit tuples for plotting
            kmer_hit_dat[rdki][1] -= my_end_err_lens[rdki]  # subtract size of artifacts at end of reads
    #
    consensus_kmer_hit_dat = []
    consensus_clust_dat    = [[],[],[],[0]] # fake data so that plot_kmer_hits doesn't complain
    consensus_tvr_tel_pos  = []
    for rdki in range(len(redrawn_consensus)):
        cons_readcount = len(read_clust_dat[0][rdki])
        cons_readname  = 'consensus-' + str(rdki) + ' [' + str(cons_readcount)
        cons_tvrlen    = my_tvr_tel_bounds[rdki]
        if cons_readcount == 1:
            cons_readname += ' read]'
        else:
            cons_readname += ' reads]'
        if cons_readcount >= MIN_READS_PER_PHASE:
            consensus_kmer_hit_dat.append([redrawn_consensus[rdki], len(my_consensus_vecs[rdki]), 0, 'FWD', cons_readname, DUMMY_TEL_MAPQ, None])
            consensus_clust_dat[0].append([rdki])
            consensus_clust_dat[1].append([DUMMY_TEL_MAPQ])
            consensus_clust_dat[2].append([0])
            consensus_tvr_tel_pos.append(cons_tvrlen)
    #
    # TVR plotting (clustered reads + consensus for each allele)
    #
    custom_plot_params = {'xlim':[-1000,15000]}
    if dpi is not None:
        custom_plot_params['dpi'] = dpi
    if telcompplot_fn is not None:
        if DO_NOT_OVERWRITE is False or exists_and_is_nonzero(telcompplot_fn) is False:
            plot_kmer_hits(kmer_hit_dat, KMER_COLORS, my_chr, my_pos, telcompplot_fn,
                           clust_dat=read_clust_dat,
                           plot_params=custom_plot_params)
    if telcompcons_fn is not None and len(consensus_clust_dat[0]):
        if DO_NOT_OVERWRITE is False or exists_and_is_nonzero(telcompcons_fn) is False:
            plot_kmer_hits(consensus_kmer_hit_dat, KMER_COLORS, my_chr, my_pos, telcompcons_fn,
                           clust_dat=consensus_clust_dat,
                           draw_boundaries=consensus_tvr_tel_pos,
                           plot_params=custom_plot_params)


def convert_colorvec_to_kmerhits(colorvecs, repeats_metadata):
    #
    [kmer_list, kmer_colors, kmer_letters, kmer_flags] = repeats_metadata
    #
    amino_2_kmer_ind = {}
    for i in range(len(kmer_colors)):
        amino_2_kmer_ind[kmer_letters[i]] = i
    #
    out_kmerhits = []
    for i in range(len(colorvecs)):
        out_kmerhits.append([[] for n in range(len(kmer_colors))])
        if len(colorvecs[i]) == 0:
            continue
        current_block = colorvecs[i][0]
        current_start = 0
        colorvecs[i] += UNKNOWN
        for j in range(1,len(colorvecs[i])):
            if colorvecs[i][j] != current_block:
                if current_block != UNKNOWN:
                    my_ind = amino_2_kmer_ind[current_block]
                    out_kmerhits[-1][my_ind].append((current_start, j))
                current_block = colorvecs[i][j]
                current_start = j
    return out_kmerhits


def plot_some_tvrs(tvrs, labels, kmer_metadata, plot_fn, custom_plot_params={}):
    #
    # helper function for quickly plotting a set of tvrs (assumes 'q' orientation)
    #
    [KMER_LIST, KMER_COLORS, KMER_LETTER, KMER_FLAGS] = kmer_metadata
    redrawn_tvrs = convert_colorvec_to_kmerhits(tvrs, kmer_metadata)
    clust_khd = []
    for i in range(len(redrawn_tvrs)):
        clust_khd.append([redrawn_tvrs[i], len(tvrs[i]), 0, 'FWD', labels[i], 60, None])
    clust_dat = []
    clust_dat.append([list(range(len(clust_khd)))])
    clust_dat.append([[60]*len(clust_dat[0][0])])
    clust_dat.append([[0]*len(clust_dat[0][0])])
    clust_dat.append([0])
    plot_kmer_hits(clust_khd, KMER_COLORS, '', 0, plot_fn, clust_dat=clust_dat, plot_params=custom_plot_params)


def readlen_plot(readlens_all, readlens_tel, readlens_final, plot_fn, xlim=(1000,1000000), nbins=60):
    mpl.rcParams.update({'font.size': 16, 'font.weight':'bold'})
    fig = mpl.figure(1,figsize=(10,8))
    #
    mpl.subplot(211)
    logbins = np.geomspace(xlim[0], xlim[1], nbins)
    mpl.hist(readlens_all, bins=logbins, color='gray')
    mpl.xscale('log')
    mpl.xlim(xlim[0], xlim[1])
    mpl.grid(which='both', linestyle='--', alpha=0.5)
    mpl.ylabel('read count')
    mpl.legend([f'{len(readlens_all)} total reads'])
    #
    mpl.subplot(212)
    logbins = np.geomspace(xlim[0], xlim[1], nbins)
    mpl.hist(readlens_tel, bins=logbins, color='blue')
    mpl.hist(readlens_final, bins=logbins, color='green')
    mpl.xscale('log')
    mpl.xlim(xlim[0], xlim[1])
    mpl.grid(which='both', linestyle='--', alpha=0.5)
    mpl.xlabel('length (bp)')
    mpl.ylabel('read count')
    mpl.legend([f'{len(readlens_tel)} initial tel reads', f'{len(readlens_final)} final tel reads'])
    #
    mpl.tight_layout()
    mpl.savefig(plot_fn)
    mpl.close(fig)
