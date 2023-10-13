import os
import pandas as pd
from os.path import join as osj

def merge_table(path_clinical,path_SF_volume,out_dir):
        '''

        1. take input from clinical table + volume table,
        merge them with 'Subject', and save
        2. prepare for the next statistic analysis

        :param path_clinical: str, path to clinical data table
        :param path_SF_volume: str, path to volume data table
        :param out_dir: str, output directory
        :return:
        '''

        #path_clinical = r'D:\hmrrc\Lmri\SampleData\Segment_HA_QA\SF_clinical.csv'
        #path_SF_volume = r'D:\hmrrc\Lmri\SampleData\Segment_HA_QA\SF_volume.csv'
        #out_path = r'D:\hmrrc\Lmri\SampleData\statisitcs\merged_data.csv'

        clinical_df = pd.read_csv(path_clinical)
        SF_volume = pd.read_csv(path_SF_volume)
        merged_df = pd.merge(on='Subject',left=clinical_df,right=SF_volume,how='outer')
        merged_df.drop(merged_df.filter(regex="Unname"),axis=1, inplace=True)
        out_filename = 'merged_data.csv'
        out_path = osj(out_dir,out_filename)
        merged_df.to_csv(out_path,index=False)

        print('data merged! please check.')
        print('merged data saved to {}'.format(out_path))
        generate_template(out_dir,out_path)
        print('then run the next step: statistics')

def generate_template(out_dir,out_path):
        '''
        generate template for statistic design
        :param out_dir: str, path to output directory
        :return:
        '''
        dict = {"DesignName":"ExampleDesign1",
                "group":"diagnosis",
                "filter":"Age>=18",
                "covariates":"ICV__age__sex",
                "x_list":"YBOCS1__YBOCS2",
                "inter_terms":"diagnosis*age__diagnosis*sex",
                "out_dir":out_dir,
                "data_path":out_path,
                "pic_size_box":"3__7",
                "legend.position.box":"none",
                "BoxColorSetting":"none",
                "pic_size_cor":"7__3",
                "legend.position.cor":"top"}

        out_template = osj(out_dir,'design.csv')
        with open(out_template, 'w') as f:
                [f.write('{0},{1}\n'.format(key, value)) for key, value in dict.items()]
        print("design template saved to {}, please fill in design and run: ".format(out_template))

def get_design(design_csv):
        '''
        get design from design.csv file
        :param design_csv: str, path to design table
        :return: a list with args
        '''
        design_df = pd.read_csv(design_csv,header=None,index_col=0)
        out_args_list = []
        for col in design_df.columns:
                args = design_df[col].to_list()
                out_args_list.append(args)

        return out_args_list


def run_R_script(Rcmd, args):
        '''

        :param Rcmd: str, cmd to run r script. e.g. 'C:/PROGRA~1/R/R-41~1.1/bin/Rscript'
        :param args: a list of args to run r script
        :return:
        '''
        import subprocess
        out_dir = osj(args[6],args[0])
        os.mkdir(out_dir)
        args[6] = out_dir
        rscript_path =os.path.join(os.path.split(os.path.realpath(__file__))[0], 'main.R')
        print('start r session: {}'.format(args[0]))

        args[2] = '"{}"'.format(args[2])
        cmd = ' '.join([Rcmd, "--vanilla",rscript_path]+args[1:])
        print(cmd)
        #subprocess.call([Rcmd, "--vanilla",rscript_path]+args[1:])
        os.system(cmd)

def Wrap_Statistic(design_csv,Rcmd,num_cores):
        '''
        run statistic design
        :param Rcmd: string, cmd used to call Rscript.
        :param design_csv: str, path to design.csv
        :return:
        '''
        from joblib import Parallel, delayed

        args_list = get_design(design_csv)

        Parallel(n_jobs=num_cores)(delayed(run_R_script)(Rcmd,args) for args in args_list)




