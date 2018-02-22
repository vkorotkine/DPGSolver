#!/usr/bin/env python3

import sys,os
import subprocess,shlex
import glob

class Cmd_Line_Params:
	""" Container for command line related parameters. """

	def __init__ (self):
		self.exec_cmd = "" ###< The name of the executable command (without arguments).
		### The root of the error output directory.
		self.output_err_dir_root = "@CMAKE_BINARY_DIR@/output/errors"

	def set_exec_cmd (self,n_procs):
		""" Set the executable command based on the number of processors to be used. """
		self.exec_cmd = "@MPIEXEC@ -n "+str(n_procs)+" @CMAKE_BINARY_DIR@/bin/test_integration_convergence"

def execute_commands (cmd_line_params,err_output_dir,test_case,petsc_opts,ctrl_spec,ctrl_matrix):
	""" Execute the commands corresponding to the input parameters. """
	for ctrl_s in ctrl_spec:
		# Remove any existing error files from previous runs of the same test case.
		files = glob.glob(err_output_dir+"/*.txt")
		for f in files:
			subprocess.run(shlex.split("rm "+f))

		for ctrl_m in ctrl_matrix:
			# Run jobs corresponding to the current settings.
			results_name = "_".join([test_case,ctrl_s])
			ctrl_name = ctrl_m+'/'+"_".join(["RESULTS",results_name])

			cmd = " ".join([cmd_line_params.exec_cmd,ctrl_name,petsc_opts,results_name])

			subprocess.call(shlex.split(cmd),cwd=str(os.getcwd()))

def run_euler_joukowski (cmd_line_params,ar_val):
	""" Run the Euler Joukowski cases. """
	err_output_dir = "/".join([cmd_line_params.output_err_dir_root,"euler","steady","joukowski"])
	test_case  = "euler_joukowski_full_dg_2D"
	petsc_opts = "petsc_options_gmres_r120"

	if (ar_val == 1):
		ctrl_spec   = ["ar1_iso", "ar1_super", ]
	elif (ar_val == 4):
		ctrl_spec   = ["ar4_iso", "ar4_super", ]
	else:
		sys.exit("Unsupported ar_val: \""+str(ar_val)+"\".")
	ctrl_matrix = ["p0-0", "p1-1", "p2-3", ]

	execute_commands(cmd_line_params,err_output_dir,test_case,petsc_opts,ctrl_spec,ctrl_matrix)

def run_navier_stokes_taylor_couette (cmd_line_params,ar_val):
	""" Run the Navier-Stokes Taylor-Couette cases. """
	err_output_dir = "/".join([cmd_line_params.output_err_dir_root,"navier_stokes","steady","taylor_couette"])


	test_case  = "navier_stokes_taylor_couette_diabatic_dg_tri"
	petsc_opts = "petsc_options_gmres_default"

	if (ar_val == 1):
		ctrl_spec   = ["ar1_iso", "ar1_super", ]
	elif (ar_val == 8):
		ctrl_spec   = ["ar8_iso", "ar8_super", ]
	else:
		sys.exit("Unsupported ar_val: \""+str(ar_val)+"\".")
	ctrl_matrix = ["p0-1", "p2-3", ]

	execute_commands(cmd_line_params,err_output_dir,test_case,petsc_opts,ctrl_spec,ctrl_matrix)


if __name__ == "__main__":
	"""
	This script can be used to execute all jobs needed to generate data for the isoparametric vs. superparametric
	investigation for the Euler and Navier-Stokes equations.

	One command line argument must be provided, specifying which jobs should be run:
	- all:                              All jobs listed below;
	- euler_joukowski_ar1:              Aspect ratio ~ 1 only.
	- euler_joukowski_ar4:              Aspect ratio ~ 4 only.
	- navier_stokes_taylor_couette_ar1: Aspect ratio ~ 1 only.
	- navier_stokes_taylor_couette_ar8: Aspect ratio ~ 8 only.
	"""

	argc = len(sys.argv)
	assert (argc == 2),"Invalid number of input arguments: "+str(argc)+" (2 required)."

	clp = Cmd_Line_Params()
	clp.set_exec_cmd(1)

	jobs_name = sys.argv[1]

	if (jobs_name == "all" or "euler_joukowski_ar1" in jobs_name):
		run_euler_joukowski(clp,1)
	if (jobs_name == "all" or "euler_joukowski_ar4" in jobs_name):
		run_euler_joukowski(clp,4)
	if (jobs_name == "all" or "navier_stokes_taylor_couette_ar1" in jobs_name):
		run_navier_stokes_taylor_couette(clp,1)
	if (jobs_name == "all" or "navier_stokes_taylor_couette_ar8" in jobs_name):
		run_navier_stokes_taylor_couette(clp,8)
