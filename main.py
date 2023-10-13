from configs import make_cfg_args
from tools import create_task


def main(cfg, args):
	func = create_task(args.task_name)
	if args.task_name == 'QA':
		root = args.data_root
		outdir = args.output
		num_cores = args.n_cpus
		plot = args.no_plot
		sub = func(root=root,num_cores=num_cores,outpath=outdir,QA = plot)

	elif args.task_name == 'merge':
		path_clinical = args.input_table
		path_SF_volume = args.SF_volume
		outdir = args.output
		sub=func(path_clinical,path_SF_volume,outdir)

	elif args.task_name == 'statistic':
		design_csv = args.design
		Rcmd = args.Rcmd
		num_cores = args.n_cpus
		sub=func(design_csv,Rcmd,num_cores)



if __name__ == '__main__':
	cfg, args = make_cfg_args()
	main(cfg, args)