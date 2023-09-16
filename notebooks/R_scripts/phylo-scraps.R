


# # Phylogenomic scrapyard

# Formating our fasta file to number entries, and replace all spaces and commas with underscores

# # # shell command to number fasta entries
# cd /central/groups/MazmanianLab/joeB/Microbiota-Immunomodulation/Microbiota-Immunomodulation/data/interim/refseq_proteins/genus_sampling

# awk '/^>/ {print ">" ++count "_" substr($0, 2)}; !/^>/ {print}' RefSeq-Genus-sample_2023-05-18.faa > RefSeq-Genus-sample_2023-05-18_numbered.faa

# # command to remove space in fasta headers and replace with underscore
# ~/bbmap/reformat.sh in=RefSeq-Genus-sample_2023-05-18_numbered.faa out=RefSeq-Genus-sample_2023-05-18_numbered-headerfixed.faa addunderscore fixjunk 

# # command to replace commas with underscore
# sed '/^>/ s/\,/_/g' RefSeq-Genus-sample_2023-05-18_numbered-headerfixed.faa > RefSeq-Genus-sample_2023-05-18_numbered-headerfixed-comma.faa




# ______________________________________________________________________________
# BULID REFSEQ DATABASE
# ______________________________________________________________________________

# mkdir -p ../data/interim/refseq_genomes/genus_sampling
# python-scripts/refseq_build.py --output ../data/interim/refseq_genomes/genus_sampling \
#     --genbank \
#     --sample 1 \
#     --rank genus \
#     --cats archaea,bacteria,fungi,invertebrate,plant,protozoa,vertebrate_mammalian,vertebrate_other \
#     --manual \
#     --above \
#     --reference \
#     --represent


# blastp -query data/interim/fastas/processed/monomers/immport/ENSG00000278567/processed_entries.fasta \
# -db /central/groups/MazmanianLab/joeB/Microbiota-Immunomodulation/Microbiota-Immunomodulation/data/interim/refseq_proteins/genus_sampling/BLASTDB/DB \
# -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
# -evalue 10 \
# -num_threads 16 \
# -max_target_seqs 100000 \
# -out realtest.out

# # res <- shell_do(
# #   glue(
# #     "mamba run -n pdmbsR seqkit stats",
# #     " /central/scratch/jbok/cleaned_fasta/GCA_000004795.1_Nlon_1.0_cleaned-protein.fasta"
# #   ),
# #   stdout_path = glue("{wkdir}/seqkit_stats.txt")
# # )


# conda activate /home/jboktor/miniconda3/envs/pdmbsR
# python-scripts/refseq_build.py --output ../data/interim/refseq_genomes/genus_sampling \
#     --genbank \
#     --sample 1 \
#     --rank genus \
#     --cats archaea,bacteria,fungi,invertebrate,plant,protozoa,vertebrate_mammalian,vertebrate_other \
#     --above \
#     --reference \
#     --represent


# ___________________________________________________________
# #### TRANSEQ EUKARYOTES

# ```{r, eval = FALSE}
# batch_input_df %>% glimpse
# batch_input_df_transeq  <- batch_input_df %>%
#   filter(is.na(refseq_prot_path)) %>%
#   filter(kingdom %nin% c("k__Bacteria", "k__Archaea")) %>%
#   filter(prodigal_unprocessed_lgr) %>% 
#   filter(!transeq_path_exists | transeq_filesize == 0 )
# batch_input_df_transeq %>% glimpse

# future::plan(
#   future.batchtools::batchtools_slurm,
#   template = glue("{wkdir}/batchtools_templates/batchtools.slurm.tmpl"),
#   resources = list(
#     name = "transeq-euk",
#     memory = "1000G",
#     ncpus = 1,
#     walltime = 36000
#   )
# )
# # Chunk files (N per job) and download
# n_jobs <- ceiling(nrow(batch_input_df_transeq) / 1)
# transeq_runs <- listenv()
# tic()
# for (job in 1:n_jobs) {
#   input_list <- chunk_func(
#     batch_input_df_transeq$refseq_fna_path, n_jobs)[[job]]
#   output_list <- chunk_func(
#     batch_input_df_transeq$transeq_path, n_jobs)[[job]]
#   transeq_runs[[job]] %<-% shell_do_transeq(
#     input_genome_list = input_list,
#     output_fasta_list = output_list
#   )
# }
# toc()

# ```












# #### PRODIGAL BACTERIA/ARCHAEA
# ```{r, eval = FALSE}
# shell_do(glue("mkdir -p /central/scratch/jbok/prodigal_proteins"))
# batch_input_df %>% glimpse
# batch_input_df_ba  <- batch_input_df %>%
#   filter(is.na(refseq_prot_path)) %>% 
#   filter(kingdom %in% c("k__Bacteria", "k__Archaea")) %>%
#   filter(prodigal_unprocessed_lgr)
# batch_input_df_ba %>% glimpse
# future::plan(
#   future.batchtools::batchtools_slurm,
#   template = glue("{wkdir}/batchtools_templates/batchtools.slurm.tmpl"),
#   resources = list(
#     name = glue("{get_time()}_prodigal-refseq_genomes"),
#     memory = "100G",
#     ncpus = 1,
#     walltime = 21600
#   )
# )
# # Chunk files (500 per job) and downlo testad
# n_jobs <- ceiling(length(batch_input_df_ba$refseq_fna_path) / 20)
# prodigal_runs <- listenv()
# tic()
# for (job in 1:n_jobs) {
#   input_list <- chunk_func(
#     batch_input_df_ba$refseq_fna_path, n_jobs
#   )[[job]]
#   output_list <- chunk_func(
#     batch_input_df_ba$prodigal_path, n_jobs
#   )[[job]]
#   prodigal_runs[[job]]  %<-% shell_do_prodigal_list(
#     input_genome_list = input_list,
#     output_fasta_list = output_list
#   )
# }
# toc()

# ```

# ___________________________________________________________
# PRODIGAL EUKARYOTES

# ```{r, eval = FALSE}
# batch_input_df %>% glimpse
# batch_input_df_euk  <- batch_input_df %>%
#   filter(is.na(refseq_prot_path)) %>% 
#   filter(kingdom %nin% c("k__Bacteria", "k__Archaea")) %>%
#   filter(prodigal_unprocessed_lgr)
# batch_input_df_euk %>% glimpse

# # Initiate future.batchtools backend for parallel processing
# future::plan(
#   future.batchtools::batchtools_slurm,
#   template = glue("{wkdir}/batchtools_templates/batchtools.slurm.tmpl"),
#   resources = list(
#     name = glue("{get_time()}_prodigal-refseq_genomes-euk"),
#     memory = "100G",
#     ncpus = 1,
#     walltime = 21600
#   )
# )
# # Chunk files (500 per job) and downlo testad
# n_jobs <- ceiling(length(batch_input_df_euk$refseq_fna_path) / 20)
# prodigal_runs <- listenv()
# tic()
# for (job in 1:n_jobs) {
#   input_list <- chunk_func(
#     batch_input_df_euk$refseq_fna_path, n_jobs
#   )[[job]]
#   output_list <- chunk_func(
#     batch_input_df_euk$prodigal_path, n_jobs
#   )[[job]]
#   prodigal_runs[[job]]  %<-% shell_do_prodigal_list(
#     input_genome_list = input_list,
#     output_fasta_list = output_list,
#     addtn_flags = " -c -g 1" # for eukaryotes only
#   )
# }
# toc()


# ```


# Files that fail prodigal processing proceed to transeq.
# ```{r, eval = FALSE}
# # # remove empty protein files
# batch_input_df %>%
#   filter(prodigal_filesize == 0) %>% #View
#   pull(prodigal_path) %>% 
#   purrr::map( ~shell_do(glue("rm {.}")))
# ```


# Checking that numbers match up

# ```{r, eval = FALSE}
# # should be 30,957 genomes / rows
# num_ncbi <- batch_input_df %>%
#   filter(!is.na(refseq_prot_path)) %>%
#   pull(refseq_prot_path) %>%
#   unique() %>%
#   length()
# num_prodigal <- batch_input_df %>%
#   drop_na(prodigal_filesize) %>%
#   filter(prodigal_filesize > 0) %>%
#   pull(prodigal_path) %>%
#   unique() %>%
#   length()
# num_transeq <- batch_input_df %>%
#   filter(transeq_path_exists) %>%
#   pull(transeq_path) %>%
#   unique() %>%
#   length()

# num_ncbi + num_prodigal + num_transeq

# all_genomes <- batch_input_df$genome %>% unique 
# a <- batch_input_df %>%
#   filter(!is.na(refseq_prot_path)) %>%
#   pull(genome) %>%
#   unique()
# b <- batch_input_df %>%
#   drop_na(prodigal_filesize) %>%
#   filter(prodigal_filesize > 0) %>%
#   pull(genome) %>%
#   unique()
# c <- batch_input_df %>%
#   filter(transeq_path_exists) %>%
#   pull(genome) %>%
#   unique()

# # NOTE: These four genomes fail to load even into transeq (too large) and are tossed
# # batch_input_df %>%
# #   filter(genome %in% setdiff(all_genomes, c(a, b, c)))

# ```