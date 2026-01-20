from matplotlib.colors import ListedColormap
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# import this script by the following lines
# import sys
# sys.path.insert(1, "/mnt/x/common/Lab_python_functions/")
# from plot_gsea import plot_GSEA

def plot_GSEA(df_orig, qval_thresh=0.25, fig_size=(5, 5), yticks_fontsize = 8, flip_colors = False,
               xticks_fontsize = 8, xlabel_fontsize = 10, dpi=None, verbose=True, return_figure = False,
               bubble_scale=2, x_padding=0.1, row_spacing=1.0, legend_round=2, 
               pos_color = sns.color_palette("hls", 5)[3], neg_color = sns.color_palette("hls", 5)[0]):
    if verbose:
        print("Q-value thresh: " + str(qval_thresh))
    
    df = df_orig.copy()
    

    # sort df by NES
    df = df.sort_values(by='NES', ascending=True)

    df['Name'] = df['Term'].str.replace('_', ' ')
    df['Gene %'] = df['Gene %'].str.rstrip('%').astype('float')

    # Filter indices based on q-value threshold
    ind = df['FDR q-val'] < qval_thresh

    if sum(ind) == 0:
        print("No gene sets pass the q-value threshold.")
        if return_figure:
            return [df, None]
        return df

    df = df.loc[ind]
    df['FDR q-val'] = np.array(df['FDR q-val'], dtype=np.float64)

    # add the minimum non zero fdr q value to the zero values to avoid log(0)
    min_non_zero_fdr_q_val = df['FDR q-val'][df['FDR q-val'] > 0].min()
    df['FDR q-val'] = df['FDR q-val'].replace(0, min_non_zero_fdr_q_val)
    df['nlog10_qval'] = -np.log10(df['FDR q-val'])
    # vmin, vmax = 0, max(df['nlog10_qval'])

    # Calculate bubble sizes for legend
    bubsizes = [df['Gene %'].min(), df['Gene %'].median(), df['Gene %'].max()]
    
    # Determine the color map
    my_cmap = ListedColormap([neg_color, pos_color])
    if flip_colors:
        my_cmap = ListedColormap([pos_color, neg_color])
    if len(np.unique(np.sign(df.NES))) == 1:
        if np.unique(np.sign(df.NES))[0] == 1:
            my_cmap = ListedColormap([pos_color]) # blue
        else:
            my_cmap = ListedColormap([neg_color]) # red

    # Create the plot
    fig, ax = plt.subplots(figsize=fig_size, dpi=dpi)
    y_positions = np.arange(len(df.NES)) * row_spacing
    scatter = ax.scatter(
        df.NES, 
        y_positions, 
        s=df['Gene %'] * bubble_scale, 
        c=np.sign(df.NES),
        cmap=my_cmap, 
        alpha=1,
        zorder=3,
    )

    # Add bubble size legend
    decimals = 2 if legend_round is None else max(int(legend_round), 0)
    for size in bubsizes:
        ax.scatter([], [], s=size * bubble_scale, color='k', alpha=1,
                   label=f'{float(size):.{decimals}f}')
    ax.legend(title='% Genes in set',fontsize=6, title_fontsize=6,bbox_to_anchor=(1.04, 0.5),
                loc="center left", frameon=False)
    
    # Set axis labels and title
    ax.set_yticks(y_positions)
    ax.set_yticklabels(df.Name, fontsize=yticks_fontsize)
    int_start = int(np.ceil(df.NES.min()))
    int_end = int(np.floor(df.NES.max()))
    if int_start <= int_end:
        if int_start % 2 != 0:
            int_start += 1
        ticks = np.arange(int_start, int_end + 1, 2)
        if len(ticks) > 0:
            ax.set_xticks(ticks)
    ax.tick_params(axis='x', labelsize=xticks_fontsize)
    ax.set_xlabel('Normalized Enrichment Score (NES)', fontsize=xlabel_fontsize)
    nes_min, nes_max = df.NES.min(), df.NES.max()
    nes_range = nes_max - nes_min
    pad = nes_range * x_padding if nes_range > 0 else 0.5
    ax.set_xlim(nes_min - pad, nes_max + pad)
    ax.grid(True, linestyle='-', alpha=0.5, zorder=-1)

    plt.tight_layout()
    plt.show()

    if return_figure:
        return [df, fig]
    return df
