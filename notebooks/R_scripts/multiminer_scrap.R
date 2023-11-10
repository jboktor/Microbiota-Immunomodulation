
#  module avail alphafold/2.2.0


slurm_run_alphafold <- function(jobname,
                                slurm_out,
                                input_fasta,
                                output_dir,
                                mode = "multimer",
                                walltime = "1-00:00",
                                use_gpu = "True",
                                mem = "128G",
                                gpus = 1,
                                cpus_per_task = 8) {
  shell_do(
    glue(
      "sbatch",
      " --job-name={jobname}",
      " --output={slurm_out}/{jobname}.out",
      " --error={slurm_out}/{jobname}.err",
      " --time={walltime}",
      " --gres=gpu:{gpus}",
      " --mem={mem}",
      " --cpus-per-task={cpus_per_task}",
      " {src_dir}/shell-scripts/alphafold.sub",
      " -g {use_gpu}",
      " -i {input_fasta}",
      " -o {output_dir}",
      " -m {mode}", 
      " -n 5" # number of models to generate
    )
  )
}


cluster_reports <- glue(
  "{wkdir}/.cluster_runs/",
  "{Sys.Date()}_AlphaFold-Multimer_TGFB-trimmedseqs_{rand_string()}"
)
message("\n\n CREATING:  ", cluster_reports, "\n")
shell_do(glue("mkdir -p {cluster_reports}"))

complex_fastas <-
  list.files(
    glue("{wkdir}/data/interim/fastas/processed/complexes"),
    # glue("{wkdir}/data/interim/alphafold-multimer/TGFB_fasta"),
    full.names = TRUE
  )

for (complex_path in complex_fastas) {
  jname <- fs::path_ext_remove(basename(complex_path))
  slurm_run_alphafold(
    jobname = jname,
    slurm_out = cluster_reports,
    input_fasta = complex_path,
    mem = "128G",
    gpus = 1,
    cpus_per_task = 8,
    output_dir = glue("{wkdir}/data/interim/alphafold-multimer/2023-02-28_TGFB_complexes/{jname}")
  )
}




# # Setting up SignalP 6.0 python package
# SIGNALP_DIR = "/home/jboktor/miniconda3/envs/signalp6/lib/python3.8/site-packages/signalp"
# slurm_shell_do(
#   glue(
#     "cp -r {homedir}/software/signalp6_fast/signalp-6-package/models/*",
#     " {SIGNALP_DIR}/model_weights/"
#   ),
#   memory = "5G",
#   ncpus = 4
# )


#______________________________________________________________________________
#  PIPER RUNS

# ./piper --help






# Piper methods

#' This algorithum is the underlying docking engine underneath ClusPro. This t
#' tools 

# $wkdir/piper -vv \ # maximim verbosity
#-c 1.0 \ # gridcell size (default) 
# -k4 \ Use first <k> eigenvalues from prm file (default: all)
# --msur_k=1.0 \  Use <k> for msur solvent and prm radii scaling (default: 1.0)
# --maskr=1.0 \  Mask all atoms within <k> angstroms of mask atoms
# -T FFTW_EXHAUSTIVE \ Use <plan type> as planning type rigor for fftw.
# -R 70000 \   Use first <nrots> rotation matrices from rotation file (default: all)
# -t 1 \  Use <t> for poisson boltzmann extrema threshold (default: 40.0)
# -p $wkdir/prms/atoms.prm \   atom parameter file
# -f $wkdir/prms/coeffs.prm   coefficient parameter file
# -r $wkdir/prms/rots.prm \   rotations parameter file
# $wkdir/example/1qqu_pnon.pdb \  receptor pdb ()
# $wkdir/example/1ba7_pnon.pdb


# 

