import os
import matplotlib.pyplot as plt
import numpy as np
import nibabel as nb


def show_slices(slices,ax = None, cmap = 'gray'):
    """ Function to display row of image slices """
    # slice_0 = epi_img_data[26, :, :]
    # slice_1 = epi_img_data[:, 30, :]
    # slice_2 = epi_img_data[:, :, 16]
    # show_slices([slice_0, slice_1, slice_2])
    if ax is None:
        fig, ax = plt.subplots(1, len(slices),figsize=(len(slices)*10,10))
    for i, slice in enumerate(slices):
       ax[i].imshow(slice.T, cmap=cmap, origin="lower")
    return fig,ax



def plot_dvars(filename,outpath):

    if os.path.exists(os.path.join(outpath,filename+'_summary-dvars.png')) is True:
        print('removing old summary plot')
        os.remove(os.path.join(outpath,filename+'_summary-dvars.png'))


    dvars = {}
    for key in  os.listdir(outpath):
        if key == 'Pre_reg':
            dvars['Orig'] = np.load(os.path.join(outpath,key,'pre_dvars.npy'))
            fd=np.load(os.path.join(outpath,key,'framewise_displacement.npy'))
        else:
            dvars[key] = np.load(os.path.join(outpath,key,'post_dvars.npy'))

    fig, ax = plt.subplots(len(dvars)+1, 1,figsize=(10,2*(len(dvars)+1) ))
    ax[0].plot(fd)
    ax[0].set_ylabel('FD')
    ax[0].text(0,sorted(fd)[-1]*0.9,
               "No.of Bad vol: " + str(np.count_nonzero(fd > 0.2)))
    ax[0].hlines(0.2, 0, fd.shape[0], colors='r', linestyles='dashed')
    i = 1
    for key in dvars:
        ax[i].plot(dvars[key])
        ax[i].set_ylabel(key)
        ax[i].text(0, sorted(dvars[key])[-1] * 0.9,
                   "No.of Bad vol: " + str(np.count_nonzero(dvars[key] > 5)))
        ax[i].hlines(5, 0, dvars[key].shape[0], colors='r', linestyles='dashed')
        i += 1

    fig.tight_layout()
    fig.savefig(os.path.join(outpath,filename+'_summary-dvars.png'))
    print('saving dvars plot to '+os.path.join(outpath,filename+'_summary-dvars.png'))
    plt.close(fig)
    return fig,ax


def plot_ts_matrix(in_ts_matrix,in_label, ax = None):
    '''
    plot timeseires matrix in a single plot, and label outliner
    :param in_ts_matrix:np.array, shape = (in_labels,timepoints)
    :param ax:  optional,axes. if not exist, then create a fig
    :return: matplotlib.figure.Figure
    '''
    max_ind = np.where(in_ts_matrix == np.max(in_ts_matrix))[0]
    min_ind = np.where(in_ts_matrix == np.min(in_ts_matrix))[0]
    ts_mean = np.mean(in_ts_matrix,axis=1)
    mean_max_ind = np.where(ts_mean == np.max(ts_mean))[0]
    mean_min_ind = np.where(ts_mean == np.min(ts_mean))[0]
    ts_sd = np.std(in_ts_matrix,axis=1)
    sd_max_ind = np.where(ts_sd == np.max(ts_sd))[0]
    sd_min_ind = np.where(ts_sd == np.min(ts_sd))[0]

    if ax is None:
        fig,ax = plt.subplots(1,1,figsize = (20,4))

    for ind,ts in enumerate(in_ts_matrix):
        label = str(in_label[ind])
        if ind in max_ind:
            ax.plot(ts,label = "Max signal: "+label)
        elif ind in min_ind:
            ax.plot(ts,label = 'Min signal: '+label)
        elif ind in mean_max_ind:
            ax.plot(ts,label = 'Max Mean: '+label)
        elif ind in mean_min_ind:
            ax.plot(ts, label='Min Mean: ' + label)
        elif ind in sd_max_ind:
            ax.plot(ts, label='Max std: ' + label)
        elif ind in sd_min_ind:
            ax.plot(ts, label='Min std: ' + label)
        else:
            ax.plot(ts)
    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
              ncol=6, mode="expand", borderaxespad=0.)
    return ax


def FS_LUT2Cm(in_file):
    '''
    convert FreeSurfer style lookup table to a DataFrame
    :param in_file:path to look up table
    :return:
    '''
    import pandas as pd
    from matplotlib import cm
    filename = os.path.splitext(in_file)[0]
    if os.path.splitext(in_file)[-1] =='.txt':
        print('receiving FS original LookUp txt, converting into csv')
        color_list = []
        with open(in_file) as f:
            for color_line in f.readlines():
                if color_line == '\n' or color_line == '':
                    pass
                elif '#' in color_line:
                    pass
                else:
                    color_list.append(color_line.replace('\n', ''))

        color_table = []
        for i in color_list:
            line = []
            for item in i.split(' '):
                if item != '':
                    line.append(item)
            color_table.append(line)

        color_df = pd.DataFrame(color_table)
        color_df.columns = ['label', 'region', 'R', 'G', 'B', 'A']
        color_df[['label', 'R', 'G', 'B', 'A']] = color_df[['label', 'R', 'G', 'B', 'A']].astype(int)
        color_df[[ 'R', 'G', 'B', 'A']] = color_df[['R', 'G', 'B', 'A']]/ 256
        color_df['A'] = 1
        color_df.loc[color_df[['R', 'G', 'B']].sum(axis = 1) ==0, 'A'] = 0
        print('saving csv for refer '+filename+'.csv')
        color_df.to_csv(filename+'.csv',index=False)
        in_file = filename+'.csv'
    else:
        print('unrecognized file type')
    return filename+'.csv'

def show_ROI(ROI,bg,view = 'in',n=3, position = 'sag',colormap = 'FreeSurfer'):

    '''
    plot regional ROI on top of background
    :param ROI: roi img
    :param bg: background img
    :param view: 'all' or 'zoom-in'
    :param fig: plt.figure
    :param axe: plt.subplots list
    :param position: sagittal, coronal or axial plane
    :param colormap:
    :return: axe list
    '''
    from matplotlib import cm
    import pandas as pd
    if type(ROI) == nb.nifti1.Nifti1Image:
        ROI_img = ROI
    elif type(ROI) == str:
        ROI_img = nb.load(ROI)

    ROI_img = nb.as_closest_canonical(ROI_img)
    ROI_matrix = ROI_img.get_fdata().astype(np.int_)
    labels = np.unique(ROI_matrix).tolist()
    viridis = cm.get_cmap('viridis', len(labels)+1)
    colors = viridis(np.linspace(0, len(labels), len(labels) + 1))

    from matplotlib.colors import ListedColormap
    if colormap == 'FreeSurfer':
        LUT = pd.read_csv(os.path.join(os.path.split(os.path.realpath(__file__))[0],'FreeSurferColorLUT.csv'))
    else:
        print('load input colormap: '+os.path.abspath(colormap))
        LUT = pd.read_csv(os.path.abspath(colormap))

    # convert ROI_matrix
    for i, label in enumerate(labels):
        ROI_matrix[np.where(ROI_matrix == label)] = i
        colors[i] = LUT.loc[LUT['label'] == label, ['R', 'G', 'B', 'A']].to_numpy()
        labels[i] = LUT.loc[LUT['label'] == label,['region']].iat[0,0]
    colormap = ListedColormap(colors)

    if type(bg) == nb.nifti1.Nifti1Image:
        bg_img = bg
    elif type(bg) == str:
        bg_img = nb.load(bg)
    bg_img = nb.as_closest_canonical(bg_img)
    bg_matrix = bg_img.get_fdata()

    ind = np.nonzero(ROI_matrix)
    if 'in' in view:
        #print('generate zoom in view')
        x_l = 2 * np.min(ind[0]) - np.max(ind[0])
        x_u = 2 * np.max(ind[0]) - np.min(ind[0])
        y_l = 2 * np.min(ind[1]) - np.max(ind[1])
        y_u = 2 * np.max(ind[1]) - np.min(ind[1])
        z_l = 2 * np.min(ind[2]) - np.max(ind[2])
        z_u = 2 * np.max(ind[2]) - np.min(ind[2])
        ROI_matrix_FOV = ROI_matrix[x_l:x_u, y_l:y_u, z_l:z_u]
        bg_matrix_FOV = bg_matrix[x_l:x_u, y_l:y_u, z_l:z_u]
        ind = np.nonzero(ROI_matrix_FOV)
    else:
        ROI_matrix_FOV = ROI_matrix
        bg_matrix_FOV = bg_matrix


    if 'sag' in position:
        a = 0
        ind[a].sort
        slice_number = np.around(np.linspace(ind[a][0], ind[a][-1], n + 2))[1:-1].astype(np.int)
        slices_ROI = [ROI_matrix_FOV[i, :, :] for i in slice_number]
        slices_bg = [bg_matrix_FOV[i, :, :] for i in slice_number]
    elif 'cor' in position:
        a = 1
        ind[a].sort
        slice_number = np.rint(np.linspace(ind[a][0], ind[a][-1], n + 2))[1:-1].astype(np.int)
        slices_ROI = [ROI_matrix_FOV[:, i, :] for i in slice_number]
        slices_bg = [bg_matrix_FOV[:, i, :] for i in slice_number]
    elif 'ax' in position:
        a = 2
        ind[a].sort
        slice_number = np.around(np.linspace(ind[a][0], ind[a][-1], n + 2))[1:-1].astype(np.int)
        slices_ROI = [ROI_matrix_FOV[:, :, i] for i in slice_number]
        slices_bg = [bg_matrix_FOV[:, :, i] for i in slice_number]


    # Todo: add fig/axis input
    fig, axis = show_slices(slices_bg)
    for ax,slice in zip(axis,slices_ROI):
        ax.imshow(slice.T, cmap=colormap, origin="lower")
    '''
    mappable = plt.cm.ScalarMappable(cmap = colormap)
    cb = fig.colorbar(mappable)
    cb.set_ticks(np.arange(0,len(labels))+0.5)
    cb.set_ticklabels(labels)
    '''
    return fig,axis






