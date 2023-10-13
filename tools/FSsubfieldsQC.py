import os
import re
import pandas as pd
from .utils import show_ROI
import matplotlib.pyplot as plt
from joblib import Parallel, delayed


def find_latest_version(filelist):
    versions = [re.search(r'v\d+',filename).group() for filename in filelist]
    versions = list(set(versions))
    versions.sort()
    return versions[-1]
'''
def get_Thalamic_filelist(root,subject):
    mri_dir = os.path.join(root, subject, 'mri')
    assert os.path.exists(mri_dir), 'mri_dir does not exist:{}'.format(mri_dir)

    filelist = [os.path.join(mri_dir,i) for i in os.listdir(mri_dir) if 'Thalamic' in i]
    version = find_latest_version(filelist)
    filelist = [i for i in filelist if version in i]

    output_dict = {'subject':subject,
                   'region':'Thalamic',
                   'version':version,
                   'mri_dir':mri_dir,
                   'volume':[os.path.join(mri_dir,i) for i in filelist if i.endswith('txt')],
                   'img':[os.path.join(mri_dir,i) for i in filelist if i.endswith('FSvoxelSpace.mgz')]}
    return output_dict
'''
def find_file_list(in_list,pattern):
    out_list = []
    for file in in_list:
        match = re.search(pattern,file)
        if match:
            out_list.append(match.group())
    return out_list

def get_HA_filelist(root,subject,logfile):
    '''
    check if all paths exist
    find highest version of output, and generate a dict with paths to use in later functions
    sink subject id that didn't have all paths
    '''
    mri_dir = os.path.join(root, subject, 'mri')
    if os.path.exists(mri_dir) is False:
        print('mri_dir does not exist:{}, ignore if it is not a subject'.format(mri_dir))
        with open(logfile,'a') as f:
            f.write('recon-all -s {}'.format(subject))

    else:
        filelist = find_file_list(os.listdir(mri_dir),'[lr]h.\w+Volumes-T\w+.v\d+.txt')
        if len(filelist) == 0:
            print(r"subject {} didn't run segment_HA, please check and re-run.".format(subject))
            with open(logfile,'a') as f:
                f.write(r'segmentHA_T1.sh -s {} '.format(subject))
                f.write('\n')
        else:
            version = find_latest_version(filelist)

            volume_files_pattern = '[lr]h.\w+Volumes-T\w+.{}.txt'.format(version)
            volume_files = find_file_list(os.listdir(mri_dir), volume_files_pattern)

            img_files_pattern = '[lr]h.hippoAmygLabels-T1.{}.[CH]\w+.FSvoxelSpace.mgz'.format(version)
            img_files = find_file_list(os.listdir(mri_dir),img_files_pattern)

            if len(volume_files) != 4 or len(img_files) != 4:
                print(r"number of output file for subject {} is wrong , please check and re-run.".format(subject))
                with open(logfile,'a') as f:
                    f.write(r'segmentHA_T1.sh -s {} '.format(subject))
                    f.write('\n')
            else:
                output_dict = {'subject':subject,
                               'region':'HA',
                               'version':version,
                               'mri_dir':mri_dir,
                               'root_dir':root,
                               'volume':[os.path.join(mri_dir,i) for i in volume_files],
                               'img':[os.path.join(mri_dir,i) for i in img_files]}
                return output_dict

def get_eTIV(root,subject):
    '''

    :param root: str, /path/to/$SUBJECTS_DIR
    :param subject: subject ID
    :return: float, eTIV or Estimated Total Intracranial Volume in mm3
    '''
    aseg_path = os.path.join(root,subject,'stats/aseg.stats')
    with open (aseg_path,'r') as f:
        aseg_content = f.readlines()
    # Measure EstimatedTotalIntraCranialVol, eTIV, Estimated Total Intracranial Volume, 1573720.795071, mm^3
    eTIV_line = [i for i in aseg_content if '# Measure EstimatedTotalIntraCranialVol' in i][0]
    eTIV = eTIV_line.split(',')[-2]
    eTIV = float(eTIV)
    return eTIV

def extract_HA_SF_volume(in_dict):
    '''
    extract volume outputs from 'hippoSfVolumes-T1.v21.txt' and 'amygNucVolumes-T1.v21.txt'
    add head and body volumes together, and return volume outcomes of single session
    :param in_dict: output from get_HA_filelist or get_Thalamic_filelist
    :return: pd.DataFrame, volume outcomes from a single subject/session
    '''
    #get eTIV
    eTIV = get_eTIV(root=in_dict['root_dir'], subject=in_dict['subject'])

    session_data = pd.DataFrame({'Subject': in_dict['subject'],'QA':0,'eTIV':eTIV}, index=['volume'])

    filelists = in_dict['volume']
    for filepath in filelists:
        side = re.search('[lr]h',filepath).group()
        data = pd.read_csv(filepath, delimiter=' ', header=None, names=['region', 'volume'])
        data.set_index('region', inplace=True)
        data = data.T
        columns = ['_'.join([side, i]) for i in data.columns]
        data.columns = columns
        session_data = session_data.merge(data, how='outer', left_index=True, right_index=True)

    for i in session_data.columns:
        # add head+body volumes together
        if 'head' in i and 'Whole' not in i:
            body = i.replace('head', 'body')
            session_data[i.replace('-head', '')] = session_data[i] + session_data[body]
    return session_data

def plot_subfield(in_dict,outpath):
    '''
    generate fig for qa
    :param in_dict: output from get_HA_filelist or get_Thalamic_filelist
    '''

    FileList = in_dict['img']
    fig_name = in_dict['subject']
    bg = os.path.join(in_dict['mri_dir'], 'nu.mgz', )

    if os.path.exists(bg) is False:
        print(bg,'not exist')
        return None
    else:
        print('{}: generating QA plots ...'.format(in_dict['subject']))
        for file in FileList:
            filename = os.path.split(file)[-1]
            png = os.path.join(outpath, fig_name+'_'+filename + '.png')
            ROI = file

            fig, ax = show_ROI(ROI, bg)
            fig.suptitle('_'.join([fig_name, file]), fontsize=30)
            fig.savefig(png)
            plt.close(fig)
        print('{}: QA plot generated!'.format(in_dict['subject']))


def post_SF_segment(root, num_cores, outpath = None, QA = True):
    # set output directory
    if outpath is None:
        outpath = os.path.join(root,'Segment_HA_QA')
    print('saving outputs to',outpath)

    if os.path.exists(outpath) is True:
        from shutil import rmtree
        print('outpath exist, deleting')
        rmtree(outpath)

    # get subjects list in root
    subjects = [i for i in os.listdir(root) if os.path.isdir(os.path.join(root, i)) is True]

    # check if all files exist, write log file if subject needs re-run
    os.mkdir(outpath)
    logfile = os.path.join(outpath,'QA_log.txt')
    with open(logfile,'w') as f:
        f.write(r'#! /bin/bash')
        f.write('\n')
        f.write(r'export SUBJECTS_DIR={}'.format(root))
        f.write('\n')
        f.write(r"# comment out lines that don't need")
        f.write('\n')

    subjects_dict =Parallel(n_jobs=num_cores)(delayed(get_HA_filelist)(root,i,logfile) for i in subjects)
    subjects_dict = list(filter(None,subjects_dict))
    SF_volume = pd.DataFrame()
    for subject_dict in subjects_dict:
        # for each subject, collect session data
        session_data = extract_HA_SF_volume(subject_dict)
        SF_volume = SF_volume.append(session_data)

    if QA is True:
        print('generating QA plots ...')
        #pool = multiprocessing.Pool(num_cores)
        #[pool.apply_async(plot_subfield,args=(subject_dict,outpath)) for subject_dict in subjects_dict]

        Parallel(n_jobs=num_cores)(delayed(plot_subfield)(subject_dict,outpath) for subject_dict in subjects_dict)

    volume_path = os.path.join(outpath,'SF_volume.csv')
    SF_volume.to_csv(volume_path,index=False)

    print('post segment_HA done!')
    print('Volume table saved to {}'.format(volume_path))
    print('please check QA plots, and fill in QA column in Volume table, then do the next step:')
    print('')
    return SF_volume



