from configs import make_cfg_args
from tools import create_task


def main(cfg, args):
	func = create_task(args.task_name)
	if args.task_name == 'QA':
		root = args.data_root
		outdir = args.output
		sub = func(root=root,outpath=outdir)

	elif args.task_name == 'merge':
		path_clinical = args.input_table
		path_SF_volume = args.SF_volume
		outdir = args.output
		sub=func(path_clinical,path_SF_volume,outdir)

	elif args.task_name == 'statistic':
		design_csv = args.design
		Rcmd = args.Rcmd
		sub=func(design_csv,Rcmd)



if __name__ == '__main__':
	cfg, args = make_cfg_args()
	main(cfg, args)