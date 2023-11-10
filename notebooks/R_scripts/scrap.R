

# ______________________________________________________________________________
# Alignment Algorithms ----
# ______________________________________________________________________________

# load in data
keyname_hits <- readRDS(
  glue(
    "{wkdir}/data/interim/tmp/",
    "2023-03-09_gget_immport.rds"
  )
)

# adding column for raw fasta path
raw_fasta_paths <-
  list.files(glue(
    "{wkdir}/data/interim/fastas/raw/monomers/immport"
  ), full.names = TRUE)

names(raw_fasta_paths) <-
  fs::path_ext_remove(basename(raw_fasta_paths))

raw_fasta_df <- list_items_to_df(
  list(raw_fasta_paths), "ensembl_id", "fasta_path_raw"
)

workflow_meta <- keyname_hits %>%
  ungroup() %>%
  left_join(raw_fasta_df, by = "ensembl_id") 
  
workflow_meta %>%
  qc_check_files(col = "fasta_path_raw")


#______________________________________________________________________________

processed_fasta_paths <-
  list.files(
    glue(
      "{wkdir}/data/interim/fastas/processed/monomers/immport"
    ),
    full.names = TRUE, recursive = TRUE, pattern = "fasta"
  )

names(processed_fasta_paths) <- processed_fasta_paths %>%
  strex::str_before_last(., pattern = "/") %>%
  strex::str_after_last(., pattern = "/")

processed_fastas_df <- list_items_to_df(
  list(processed_fasta_paths), "ensembl_id", "fasta_path_processed"
)
workflow_meta %<>%
  left_join(processed_fastas_df, by = "ensembl_id")
workflow_meta %>% glimpse

workflow_meta %>%
  qc_check_files(col = "fasta_path_processed")

saveRDS(workflow_meta, file = glue(
  "{wkdir}/data/interim/tmp/",
  "{Sys.Date()}_workflow-meta-1_FastaQC.rds"
))



slurm_shell_do(
  cmd =
    glue(
      "hhblits -cpu 8",
      " -i {wkdir}/data/interim/fastas/processed/monomers/immport/ENSG00000176697/processed_entries.fasta",
      " -d /central/groups/MazmanianLab/joeB/Downloads/RefDBs/",
        "HHsuite3db/UniRef30_2022_02/UniRef30_2022_02",
      " -o ENSG00000176697.hhr",
      " -oa3m ENSG00000176697.a3m",
      " -n 2", # number of iterations
      " -diff 1000", # maximum number of sequences to keep in MSA
      " -id 90" # maximum similarity to input sequence
    ),
  memory = "10G",
  ncpus = 4,
  walltime = 3600
)







# ______________________________________________________________________________
# Visualizing MSA and fasta files for quality
outdir <- glue("{wkdir}/data/interim/hhblits_results/immport")


# remove files that failed runs and re-run
# list.files(
#   glue("{wkdir}/data/interim/hmmersearch_results/hhblits_pHMMs/immport"),
#   full.names = TRUE #, pattern = "_tblout.txt"
# ) %>%
#   purrr::set_names() %>%
#     purrr::map(~ file.size(.)) %>%
#     keep(~ . < 1) %>%
#     names() %>%
#     purrr::map(~ file.remove(.))





hmmer_results_fsizes <- list.files(
  glue("{wkdir}/data/interim/hmmersearch_results/hhblits_pHMMs/immport"),
  full.names = TRUE, pattern = "_tblout.txt"
) %>%
  purrr::set_names() %>%
    purrr::map(~ file.size(.)) %>%
    purrr::map(~ data.frame("file_size" = .)) %>%
  dplyr::bind_rows(.id = "hmmsearch_path")

hmmer_results_fsizes %>%
  ggplot(aes(x = file_size)) +
  scale_x_log10() +
  geom_density()













# p_statint_all_hits <- hmmer_results_df %>% 
#   ggplot(aes(
#     x = -log10(sequence_evalue + impt),
#     y = fct_reorder(gene_name, -sequence_evalue)
#     )) +
#   ggdist::stat_interval(
#     .width = c(.25, 0.5, .75, 1),
#     show.legend = T,
#     height = 0.0001
#   ) +
#   scale_color_brewer(palette = "PuBuGn") +
#   theme_light() +
#   facet_grid(transdomain_blastp_signal ~ catalog, 
#     space =  "free", scales = "free_y") +
#   labs(x = "-log10(E-value)", y = NULL) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.grid = element_blank()
#     )


# module load singularity/3.8.0
# singularity shell \
#  -B /central:/central \
#  /central/groups/MazmanianLab/joeB/docking/masif_seed
# cd /central/groups/MazmanianLab/joeB/git/masif_seed/masif/data/masif_site
# ./data_prepare_one.sh 2PJY_BC
# ./predict_site.sh 2PJY_BC
# ./color_site.sh 2PJY_BC

# ------- with AFmultimer model
# ./data_prepare_one.sh \
# --file /central/groups/MazmanianLab/joeB/Microbiota-Immunomodulation/Microbiota-Immunomodulation/data/interim/alphafold-multimer/2023-02-28_TGFB_complexes/TGFB3_TGFBR1_TGFBR2/TGFB3_TGFBR1_TGFBR2/ranked_0.pdb AF2mTGFBR1TGFBR2_BC
# ./predict_site.sh AF2mTGFBR1TGFBR2_BC && ./color_site.sh AF2mTGFBR1TGFBR2_BC



## AN example masif_seed search
# cd /central/groups/MazmanianLab/joeB/git/masif_seed/masif_seed_search/data/masif_targets
# ./run_target_protocol.sh 2PJY_BC
# cd targets/2PJY_BC/
# ./run.sh 2PJY_BC


# 5G file
# data/interim/jackhmmer_results/immport/ENSG00000162747/ENSG00000162747_msa.fasta

# input fasta
# data/interim/fastas/processed/monomers/immport/ENSG00000162747/processed_entries.fasta



fasta_paths_df <- list(fastas_processed_keyname_hits) %>% 
  bind_rows(.id = "id") %>%
  tidyr::pivot_longer(!id, names_to = "ensembl_id", values_to = "fasta_path") %>% 
  dplyr::select(-id) %>% 
  mutate(fasta_path = glue("{fasta_path}/processed_entries.fasta"))

saveRDS(fasta_paths_df, glue(
  "{wkdir}/data/interim/tmp/",
  "{Sys.Date()}_fasta-paths_processed_keyname_hits-data-frame.rds"
))

protein_catalog_path <- "/central/groups/MazmanianLab/joeB/Downloads/protein_catalogs"
chunked_catalog_paths <- list.files(
  glue("{protein_catalog_path}/clustered_catalogs/merged/chunked_fasta"),
  full.names = TRUE) %>% 
  keep(grepl("chunk", .))

# create a datatable mapping each ensembl_id to each chunked catalog
binned_eids <- data.frame(
  "ensembl_id" = rep(
    keyname_hits$ensembl_id,
    each = length(chunked_catalog_paths)
  ),
  "refdb_path" = rep(
    chunked_catalog_paths,
    times = length(keyname_hits$ensembl_id)
  )
)

alignment_df_list <- bind_rows(
  keyname_hits %>%
    dplyr::select(ensembl_id) %>%
    mutate(method = "needle"),
  keyname_hits %>%
    dplyr::select(ensembl_id) %>%
    mutate(method = "water")
    ) %>%
    ungroup() %>%
    full_join(binned_eids, by = "ensembl_id") %>% 
    full_join(fasta_paths_df, by = "ensembl_id") %>%
    select(-c(gene_name, sequences_aa))

alignment_df_list %>% glimpse()

saveRDS(
  alignment_df_list,
  glue("{wkdir}/data/interim/tmp/emboss-alignment-execution-df.rds")
)

alignment_df_list %>% glimpse()
alignment_df_list$ensembl_id %>% unique() %>% length()
alignment_df_list$fasta_path %>% unique() %>% length()
alignment_df_list$method %>% unique() %>% length()
alignment_df_list$refdb_path %>% unique() %>% length()

alignment_df_list$fasta_path[1]


#______________________________________________________________________________
# EXECUTE ALIGNMENT JOBS ----


alignment_df_list <- readRDS(
  glue("{wkdir}/data/interim/tmp/emboss-alignment-execution-df.rds")
)
# PRIORITIZING RUNS
keyname_hits <- readRDS(
  glue("{wkdir}/data/interim/tmp/",
  "2023-03-09_gget_immport.rds")
  )
# View(keyname_hits)

gene_df <- read.delim(
  glue("{wkdir}/data/input/Immport_gene_lists/ImmportGeneList.txt"),
      stringsAsFactors = FALSE, header = TRUE
      )
go_gene_df <- read.delim(
  glue("{wkdir}/data/input/Immport_gene_lists/ImmportGeneListGOAnnotation.txt"),
      stringsAsFactors = FALSE, header = TRUE
      )
immune_goi <- gene_df %>%   
  filter(
    !grepl("HLA", Symbol),
    !grepl("immunoglobulin", Name),
    !grepl("IGK|IGL", Symbol),
    !grepl("T cell receptor", Name)
    ) %>% 
  distinct()




# immune_goi %>% filter(Symbol == "IL10")
# priority_genes %>% filter(Symbol == "IL10")

priority_genes <- immune_goi %>% 
  filter(
    Category %in% c(
      "Cytokines" #, 
      # "Chemokines" #, 
      # "Interferons", 
      # "Interleukins", 
      # "TGFb_Family_Member" #, 
      # "TNF_Family_Members"
      )) %>% 
    # filter(grepl("IL|CXC", Symbol )) %>%
    # grepl("IL10|TGFB|IL1RN|IL17|IL35|TNF", Symbol )) %>%
    dplyr::select(-c(Category, Chromosome)) %>% 
    distinct() %>% 
    left_join(keyname_hits, by = c("Symbol" = "gene_name")) %>% 
    drop_na(ensembl_id)

priority_genes$ensembl_id %>% unique()
priority_genes$Symbol %>% unique()
# priority_genes %>% 
#   filter(Symbol == "IL1RN")

# Removing runs that have already been completed
emboss_output_dir <- "/central/scratch/jbok/emboss_alignments"
emboss_output_files <- list.files(emboss_output_dir, full.names = TRUE)

# make a temporary column with the file output name
alignment_df_list_tmp <- 
  alignment_df_list %>% 
    mutate(
      output_file = glue(
        "{emboss_output_dir}/",
        "{ensembl_id}_",
        "{fs::path_ext_remove(basename(refdb_path))}.{method}"
        ),
      output_rds = glue("{output_file}.rds")
    ) %>% 
    # SELECT PRIORITIZED RUNS
    filter(ensembl_id %in% priority_genes$ensembl_id) %>%
      # REMOVE RUNS THAT HAVE ALREADY BEEN COMPLETED
      filter(output_file %nin% emboss_output_files) %>%
      filter(output_rds %nin% emboss_output_files)


alignment_df_list_tmp %>% glimpse()
alignment_df_list_tmp$ensembl_id %>% unique() 
alignment_df_list_tmp %>%
  left_join(priority_genes, by = "ensembl_id") %>% # View()
  pull(Symbol) %>%
  table()

# arrange so that genes that are the closest to being completed are first
gene_rank_priority <- alignment_df_list_tmp %>%
  group_by(ensembl_id) %>%
  summarize(n = n())


# Finalized batchtools execution dataframe
batchtools_params <- alignment_df_list_tmp %>%
  select(ensembl_id, fasta_path, refdb_path, method) %>%
  left_join(gene_rank_priority, by = "ensembl_id") %>%
  arrange(n) %>%
  select(-n)
batchtools_params %>% glimpse()

# configure registry ----
cluster_run <- glue("{Sys.Date()}_EMBOSS-Alignment_ID-{rand_string()}/")
message("\n\nRUNNING:  ", cluster_run, "\n")
breg <- makeRegistry(
  file.dir = glue(
    "{wkdir}/.cluster_runs/",
    cluster_run
  ),
  seed = 42
)
breg$cluster.functions <- batchtools::makeClusterFunctionsSlurm(
  template = glue("{wkdir}/batchtools_templates/batchtools.slurm.tmpl"),
  scheduler.latency = 0.1,
  fs.latency = 65
)
# Submit Jobs ----
jobs <- batchMap(
  fun = align_fasta_sequences,
  args = batchtools_params,
  reg = breg
)
setJobNames(jobs,
  paste("EMBOSS",
  batchtools_params$ensembl_id,
  batchtools_params$method,
  sep = "_"),
  reg = breg
)
getJobNames(jobs, reg = breg)
submitJobs(jobs,
  resources = list(
    walltime = 172800,
    memory = 75000,
    ncpus = 1,
    max.concurrent.jobs = 9999
  )
)

# waitForJobs(sleep = 3)


#______________________________________________________________________________
# Aggregate Alignment Files ----

keyname_hits <- readRDS(
  glue("{wkdir}/data/interim/tmp/",
  "2023-03-09_gget_immport.rds")
  )
# Collect list of files to aggregate
chunked_ouput_files <-list.files(glue("/central/scratch/jbok/emboss_alignments"), pattern = "rds$")

chunked_ouput_df <- 
  data.frame(
    file_path = chunked_ouput_files,
    ensembl_id = chunked_ouput_files %>% map_chr(~str_split(., "_")[[1]][1])
  ) %>% 
  mutate(
    method = case_when(
    grepl("needle", file_path) ~ "needle",
    grepl("water", file_path) ~ "water",
    TRUE ~ "ERROR"
  )) 


View(chunked_ouput_stats_df)

# View(chunked_ouput_df)
# priority_genes
#' Create a datatable consisting of all files to aggregate 
#' needs to have columns: id, output_path_results, output_dir_db, input_files
#' This dataframe will be loaded available output files that match a criteria (have all 50 rds files completed)
#' and then filtered to remove existing processed files
#' There will need to be an intermidate data frame that contains the desired output files 
#' and a logical filter to remove rows that have already been processed

# list of output files that have been processed
agg_processed_files <- list.files(
  glue("{wkdir}/data/processed/emboss_alignments/immport"), 
  pattern = "rds$"
  )

# Incomplete runs
chunked_ouput_stats_df %>% 
  select(ensembl_id, n_chunks, method) %>% 
  filter(n_chunks != 50) %>% 
  print(n = Inf)

# intermediate_df <- 
#   # chunked_ouput_df %>% 
#   # nest(input_files = c("file_path")) %>%
#   # left_join(chunked_ouput_stats_df, by = "ensembl_id") %>% 
#   chunked_ouput_df %>% 
#     filter(n_chunks == 50)  %>% 
#     select(ensembl_id, method, file_path) %>%
#     dplyr::rename(input_files = file_path) %>%
#     # tidyr::nest(.by = c("ensembl_id", "method")) %>% 
#     mutate(id =
#       glue(
#         "{ensembl_id}_2023-02-13_catalog_",
#         "MMSeq2-95_rep_seq.{method}"
#       ),
#       expected_output_path = glue(
#         "{results_dir}/{id}.rds"
#       ),
#       output_dir_db = sql_dir,
#       output_path_results = results_dir
#       ) %>% 
#       filter(expected_output_path %nin% agg_processed_files) 



future::plan(
  list(
    tweak(
      future.batchtools::batchtools_slurm,
      template = glue(
        "{wkdir}/batchtools_templates/",
        "batchtools.slurm.tmpl"
      ),
      resources = list(
        name = glue("aggregate-rds_{rand_string()}"),
        memory = "8G",
        ncpus = 4,
        walltime = 108000
      )
    ),
    multisession
  )
)

chunked_ouput_stats_df <-
  chunked_ouput_df %>%
  group_by(ensembl_id, method) %>%
  dplyr::summarise(n_chunks = n()) %>%
  left_join(priority_genes, by = "ensembl_id") %>%
  filter(n_chunks == 50) %>%
  mutate(
    expected_output_file = glue(
      "{wkdir}/data/processed/emboss_alignments/immport/",
      "{ensembl_id}_2023-02-13_catalog_MMSeq2-95_rep_seq.{method}.rds"
    ),
    output_exists = map_lgl(expected_output_file, ~ (file.exists(.)))
  )


eids <- chunked_ouput_stats_df$ensembl_id
file_dir <- glue("/central/scratch/jbok/emboss_alignments")
file_paths <- list.files(file_dir, pattern = "rds$", full.names = TRUE)
sql_dir <- "/central/scratch/jbok/emboss_alignments-SQL-DB"
results_dir <- glue({"{wkdir}/data/processed/emboss_alignments/immport"})
shell_do(glue("mkdir -p {sql_dir}"))

for (eid in eids) {
  message("Processing: ", eid)
  for (method in c("needle", "water")) {
    name_out <-
      glue(
        "{eid}_2023-02-13_catalog_",
        "MMSeq2-95_rep_seq.{method}"
      )
    # check if output results exist
    output_check <- file.exists(glue("{results_dir}/{name_out}.rds"))
    if (output_check) {
      message("Skipping: ", name_out, " file exists ...")
      next
    }
    message("Aggregating: ", method, " ", eid)
    files <- file_paths %>%
      keep(grepl(
        glue(
          "{eid}_2023-02-13_catalog_",
          "MMSeq2-95_rep_seq_.*{method}.rds"
        ), .
      ))
    if (length(files) < 50) {
      message(
        "ERROR: ", length(files),
        " files to aggregate for: ", name_out, "\n"
      )
      next
    }
    agg_files_needle %<-% {
      aggregate_alignment_files(
        input_files = files,
        output_dir_db = sql_dir,
        output_path_results = results_dir,
        id = name_out
      )
    }
  }
}











# List of processed microbial alignment result files
microbial_alignment_files <- list.files(
  glue("{wkdir}/data/processed/emboss_alignments/immport"),
  pattern = "rds$", full.names = TRUE
) %>% as.list()
names(microbial_alignment_files) <- 
  microbial_alignment_files %>% map_chr(~  str_split(basename(.), "_")[[1]][1])

# List of processed self alignment result files
self_alignment_files <- list.files(
  glue("/central/scratch/jbok/emboss_alignments_immport_self"),
  pattern = "rds$", full.names = TRUE
) %>% as.list()
names(self_alignment_files) <- 
  self_alignment_files %>% map_chr(~ str_split(basename(.), "_")[[1]][1])

# Join the lists
processed_alignment_results <-
  full_join(
    microbial_alignment_files %>%
      grep("needle", ., value = TRUE) %>%
      list_items_to_df("ensembl_id", "needle_path_microbial"),
    microbial_alignment_files %>%
      grep("water", ., value = TRUE) %>%
      list_items_to_df("ensembl_id", "water_path_microbial")
  ) %>%
  full_join(
    self_alignment_files %>%
      grep("needle", ., value = TRUE) %>%
      list_items_to_df("ensembl_id", "needle_path_self")
  ) %>%
    full_join(
      self_alignment_files %>%
        grep("water", ., value = TRUE) %>%
        list_items_to_df("ensembl_id", "water_path_self")
    )

workflow_meta <- readRDS(
  glue(
    "{wkdir}/data/interim/tmp/",
    "2023-03-22_workflow-meta-3_HmmersearchQC.rds"
  )
)

alignment_paths_df <- workflow_meta %>%
  left_join(processed_alignment_results) %>%
  select(ensembl_id, gene_name, contains(c("needle", "water"))) %>%
  drop_na(contains("path"))

eids <- alignment_paths_df %>%
  filter(grepl("IL27", gene_name)) %>%
  pull(ensembl_id)

alignment_paths_df$gene_name %>% unique()
top_immport_hits
# intersect(top_immport_hits, alignment_paths_df$gene_name %>% unique())

for (eid in eids) {
  goi <- alignment_paths_df %>% 
  filter(ensembl_id == eid) %>%
  pull(gene_name)
  message("Processing: ", goi)
  output_fig <- glue(
  "{wkdir}/figures/sequence-alignment/",
  "{Sys.Date()}_{goi}_{eid}_patch-alignment_data.png"
  )
  if (file.exists(output_fig)) {
    message("Skipping: ", goi, " file exists ...")
    next
  }

# Loading in data
needle_data <- readRDS(
  filter(alignment_paths_df, gene_name == goi)$needle_path_microbial
)
water_data <- readRDS(
  filter(alignment_paths_df, gene_name == goi)$water_path_microbial
)
needle_data_self <- readRDS(
  filter(alignment_paths_df, gene_name == goi)$needle_path_self
) %>%
  mutate(method = "Needleman-Wunsch (Global)") %>% 
  slice_sample(n = 10000)
water_data_self <- readRDS(
  filter(alignment_paths_df, gene_name == goi)$water_path_self
) %>%
  mutate(method = "Smith-Waterman (Local)") %>% 
  slice_sample(n = 10000)

  # Wrangling data
  needle_df <- needle_data$subsample %>%
    mutate(method = "Needleman-Wunsch (Global)")
water_df <- water_data$subsample %>%
  mutate(method = "Smith-Waterman (Local)")

needle_all_df <- bind_rows(
  needle_data$subsample,
  needle_data$top_percent_identity,
  needle_data$top_percent_similarity
) %>%
  distinct() %>%
  mutate(method = "Needleman-Wunsch (Global)")

water_all_df <- bind_rows(
  water_data$subsample,
  water_data$top_percent_identity,
  water_data$top_percent_similarity
) %>%
  distinct() %>%
  mutate(method = "Smith-Waterman (Local)")

  # Plotting 
p_hist_needle <- ggplot() +
  geom_histogram(
  data = needle_data_self, aes(x = Score),
  bins = 150, fill = "grey", alpha = 0.7) +
  ggpk_dist(data = needle_df, aes(x = Score)) +
  scale_x_log10()

p_hist_water <- ggplot() +
  geom_histogram(
  data = water_data_self, aes(x = Score),
  bins = 150, fill = "grey", alpha = 0.7) +
  ggpk_dist(data = water_df, aes(x = Score)) +
  scale_x_log10()

p_lvs_needle <- 
  ggplot(needle_all_df, aes(percent_similarity, Length)) +
    ggpk_length_vs_sim()

p_lvs_water <- 
  ggplot(water_all_df, aes(percent_similarity, Length)) +
    ggpk_length_vs_sim()

target_cols <-
  c("Length", "Score", "percent_identity", "percent_similarity", "percent_gaps")

alignment_data <- bind_rows(water_all_df, needle_all_df)
alignment_wide <- alignment_data %>%
  distinct() %>%
  group_by(protein_1, protein_2, method) %>%
  dplyr::slice_max(order_by = Score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(protein_1, protein_2, method, target_cols) %>%
  pivot_wider(names_from = method, values_from = target_cols)


p_needle_vs_water <- alignment_wide %>%
  filter(`Length_Smith-Waterman (Local)` >= 9) %>%
  ggplot(aes(
    x = `percent_similarity_Smith-Waterman (Local)`, 
    y = `percent_similarity_Needleman-Wunsch (Global)`)) +
  geom_pointdensity(adjust = 1, size = 0.25) +
  scale_color_viridis(option = "G") +
  geom_ysidedensity(aes(x = stat(density))) +
  geom_xsidedensity(aes(y = stat(density))) +
  scale_ysidex_continuous(guide = "none") +
  scale_xsidey_continuous(guide = "none") +
  labs(
    x = "Percent Similarity (Local)",
    y = "Percent Similarity (Global)",
    title = "") +
  theme_light() +
    theme(
      ggside.panel.scale.x = .17,
      ggside.panel.scale.y = .085,
      legend.position = "bottom",
      legend.key.width = unit(1, 'cm')
    )

patch <-
  (p_hist_needle + p_hist_water) /
    (p_lvs_needle + p_lvs_water) /
    p_needle_vs_water

patch_edit <- patch +
  plot_layout(heights = c(1.5, 1.5, 2)) +
  plot_annotation(
    tag_levels = "A",
    title = glue("{eid} - {goi} Microbial Catalog Alignment Profile")
  )
ggsave(
  output_fig,
  patch_edit,
  width = 10, height = 15, dpi = 600
)
  }



# _______________________________________________________________________________
# HMMER Analysis ----
# _______________________________________________________________________________







# _____________________________________________________________________________
# HMMER Search ----
# Hmmsearch run on MSA profile against custom database


# Visualizing HMMER Search Results ----



hmmer_results_paths <- list.files(
  "/central/scratch/jbok/mim_temp/hmmersearch_results/immport",
  recursive = TRUE, full.names = TRUE, pattern = "_tblout.txt"
  ) %>%
  keep(file.size(.) > 0)

# load in hmmersearch results
hmmer_results <- map_df(hmmer_results_paths, read_tblout)

hmmer_results_df <- hmmer_results %>%
  select(-c(domain_accession, query_accession)) %>%
  mutate(catalog = map_chr(
    strex::str_before_nth(domain_name, "_", 2),
    ~ str_remove(., "CATID_")
  ))

saveRDS(
  hmmer_results_df,
  glue(
    "{wkdir}/data/processed/hmmersearch_results/",
    "{Sys.Date()}_hmmersearch_results.rds"
  )
)


workflow_meta %<>%
  mutate(
    hmmersearch_output_expected = glue(
      "/central/scratch/jbok/mim_temp/",
      "hmmersearch_results/immport/{ensembl_id}/{ensembl_id}_tblout.txt"
    )) %>%
    mutate(hmmersearch_output_created = 
      map_lgl(hmmersearch_output_expected, ~ file.exists(.))) %>%
    mutate(hmmersearch_output_exists = 
      map_lgl(hmmersearch_output_expected, ~ as.logical(file.size(.) > 0)))

saveRDS(
  workflow_meta,
  glue(
    "{wkdir}/data/interim/tmp/",
    "{Sys.Date()}_workflow-meta-3_HmmersearchQC.rds"
  )
) 

#______________________________________________________________________________

View(top_immport_hits)

# What are the top hits, across all catalogs?
top_immport_hits <- analysis_df %>%
  filter(catalog != "Wormbase") %>%
  group_by(query_name) %>%
  slice_min(order_by = sequence_evalue, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  slice_min(order_by = sequence_evalue, n = 30) %>%
  pull(gene_name)

top_5_immport_hits_per_catalog <- analysis_df %>%
  group_by(query_name, catalog) %>%
  slice_min(order_by = sequence_evalue, n = 5, with_ties = FALSE) %>%
  # ungroup() %>%
  group_by(catalog) %>%
  slice_min(order_by = sequence_evalue, n = 10, with_ties = FALSE) %>%
  pull(gene_name) %>% unique()
top_5_immport_hits_per_catalog

top_immport_hits_top_sig <- analysis_df %>%
  filter(gene_name %in% top_5_immport_hits_per_catalog) %>%
  group_by(gene_name, catalog) %>%
  slice_min(order_by = sequence_evalue, n = 10)

protein_catalog_path <-
  glue(
    "/central/groups/MazmanianLab/joeB/Downloads/",
    "protein_catalogs/clustered_catalogs/merged/",
    "2023-02-13_catalog_MMSeq2-95_rep_seq.fasta"
  )
# Extracting top hits from protein catalog
for (gene_hit in top_immport_hits) {
  domain_ids <- top_immport_hits_top_sig %>%
    filter(gene_name == gene_hit) %>%
    pull(domain_name) %>%
    unique()

  write.table(
    as.data.frame(domain_ids),
    file = "temp_headerlist.txt",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  output_fasta <- glue(
    "{wkdir}/data/interim/alignment_hits/",
    "2023-03-23_{gene_hit}_top-10-per-subcatalog.fasta"
  )
  if (file.exists(output_fasta)) {
    next
  }
  shell_do(
    glue(
      "conda run -n pdmbsR seqtk subseq {protein_catalog_path}",
      " temp_headerlist.txt > {output_fasta}"
    )
  )
}




hits <- analysis_df %>%
  filter(catalog == "UHGP") %>%
  slice_min(order_by = sequence_evalue, n = 10) %>%
  pull(domain_name) %>% unique()






# hits <-
analysis_df %>%
  group_by()
  slice_min(order_by = sequence_evalue, n = 10) %>%
  pull(domain_name) %>% unique()

write.table(
  as.data.frame(hits),
  file = "fasta_header_top-UHGP.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


workflow_meta$hmmersearch_output_exists


# library(ggdist)

# analysis_df %>% glimpse
# analysis_df %>%
#       ggplot(aes(
#         x = -log10(sequence_evalue), 
#         y = -log10(best_domain_evalue),
#         group = description
#         )) +
#       scale_color_manual(values = catalogs_cols) +
#       geom_point(aes(color = catalog)) +
#       theme_bw()


# library(tidybayes)

p_statint <- analysis_df %>%
  ggplot(aes(
    x = -log10(sequence_evalue + impt),
    y = fct_reorder(gene_name, -sequence_evalue)
    )) +
  ggdist::stat_interval(
    .width = c(.25, 0.5, .75, 1),
    show.legend = T,
    height = 0.0001
  ) +
  scale_color_brewer(palette = "PuBuGn") +
    theme_light() +
    facet_wrap(~catalog, nrow = 1) +
    labs(x = "-log10(E-value)", y = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
    )
ggsave("hmmsearch_statintervals.png",
  p_statint, width = 12, height = 80, dpi = 600,
  limitsize = FALSE
)


p_statint_tophits <- analysis_df %>%
  filter(gene_name %in% top_immport_hits) %>% 
  ggplot(aes(
    x = -log10(sequence_evalue + impt),
    y = fct_reorder(gene_name, -sequence_evalue)
    )) +
  ggdist::stat_interval(
    .width = c(.25, 0.5, .75, 1),
    show.legend = T,
    height = 0.0001
  ) +
  scale_color_brewer(palette = "PuBuGn") +
    theme_light() +
    facet_wrap(~catalog, nrow = 1) +
    labs(x = "-log10(E-value)", y = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
    )

ggsave(
  glue(
    "{wkdir}/figures/hmmersearch/",
    "{Sys.Date()}_hmmsearch_statintervals-top30-immport-hits.png"
    ),
  p_statint_tophits,
  width = 12, height = 7, dpi = 600
)


p_statint_tophits_percat <- analysis_df %>%
  filter(gene_name %in% top_5_immport_hits_per_catalog) %>% 
  ggplot(aes(
    x = -log10(sequence_evalue + impt),
    y = fct_reorder(gene_name, -sequence_evalue)
    )) +
  ggdist::stat_interval(
    .width = c(.25, 0.5, .75, 1),
    show.legend = T,
    height = 0.0001
  ) +
  scale_color_brewer(palette = "PuBuGn") +
    theme_light() +
    facet_wrap(~catalog, nrow = 1) +
    labs(x = "-log10(E-value)", y = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
    )

ggsave(
  glue(
    "{wkdir}/figures/hmmersearch/",
    "{Sys.Date()}_hmmsearch_statintervals-top5-immport-hits-per-catalog.png"
    ),
  p_statint_tophits_percat,
  width = 12, height = 7, dpi = 600
)




analysis_df %>% 
  # filter(
  #   Category_Chemokines == 1 |
  #     Category_Chemokine_Receptors == 1 |
  #     Category_Cytokines == 1 #|
  #     Category_Cytokine_Receptors == 1
  #     ) %>% #
  filter(grepl("IL|TGFB|TLR", gene_name)) %>% 
  ggplot(aes(
    x = -log10(sequence_evalue + impt),
    y = fct_reorder(gene_name, -sequence_evalue)
    )) +
  ggdist::stat_interval(
    .width = c(.25, 0.5, .75, 1),
    show.legend = TRUE,
    height = 0.0001
  ) +
  scale_color_brewer(palette = "PuBuGn") +
    theme_light() +
    facet_wrap(~catalog, nrow = 1) +
  theme(
    # axis.text.y = element_blank(),
    panel.grid = element_blank()
    )

ggsave("hmmsearch_statintervalsasdfasdf.png",
  p_statint, width = 9, height = 80, dpi = 600,
  limitsize = FALSE
)


p_pointint <- analysis_df %>%
  ggplot(aes(
    x = log10(sequence_score + 1),
    y = fct_reorder(gene_name, sequence_score)
    )) +
    stat_pointinterval(stroke = 0.1)
ggsave("hmmsearch_pointintervals.png",
  p_pointint, width = 9, height = 40, dpi = 600,
  limitsize = FALSE
)




# _____________________________________________________________________________

# Pulling a specific sequence from catalog


seq_df <- readRDS(
  glue(
    "{wkdir}/data/interim/tmp/",
    "2023-04-26_signalp-trimmed-sequence-df.rds"
  )
)

hit_headers <- analysis_df %>%
  filter(catalog == "UHGP") %>%
  filter(gene_name == "C5") %>% 
  slice_min(order_by = sequence_evalue, n = 10) %>%
  pull(domain_name)

# extract top hits into a fasta file
fasta_dir <- glue("{wkdir}/data/interim/alignment_hits")
catalog_hits_fasta <- glue("{fasta_dir}/{Sys.Date()}_C5-UHGP_10hits.fasta")
extract_seq_from_catalog(
  fasta_header_list = hit_headers,
  output_fasta_path = catalog_hits_fasta,
  catalog_path = protein_catalog_path
)

# append original sequence into file
trimmed_seq <- seq_df %>% filter(grepl("^C5_", id))
seqinr::write.fasta(
  sequences = trimmed_seq$sequences_aa_signalp_trimmed,
  names = trimmed_seq$id,
  file.out = catalog_hits_fasta,
  open = "a",
  nbchar = 1e100,
  as.string = FALSE
)

library(ggmsa)
tic()
msa_plot <- ggmsa(
  msa = glue("{wkdir}/data/interim/alignment_hits/2023-04-26_C5-UHGP_10hits.fasta"),
  # char_width = 0.25, 
  seq_name = TRUE) +
  geom_seqlogo() +
  geom_msaBar()
toc()

ggsave(
  glue("{wkdir}/figures/MSA/2023-04-26_C5-UHGP_10hits.png"),
  msa_plot, 
  width = 80, height = 6, dpi = 600, limitsize = FALSE
)

seq_list <- protr::readFASTA(catalog_hits_fasta)

parseq <- protr::parSeqSim(
  seq_list,
  cores = 2,
  batches = 1,
  verbose = FALSE,
  type = "local",
  submat = "BLOSUM62",
  gap.opening = 10,
  gap.extension = 4
)
rownames(parseq) <- names(seq_list)
colnames(parseq) <- names(seq_list)


long_df <- parseq %>%
  as.data.frame() %>%
  rownames_to_column("seq1") %>%
  pivot_longer(
   !seq1,
    names_to = "seq2",
    values_to = "similarity"
  )

  long_df %>%
    ggplot(aes(x = seq1, y = seq2, fill = similarity)) +
    geom_tile() +
    geom_text(aes(label = round(similarity, 2)), size = , color = "white") +
    scale_fill_viridis_c(option = "F") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )


# _____________________________________________________________________________

# Run Self-Alignment

# Select MSAs that are complete and convert them into a fasta file


# Apply alignment algorithm to the jackhmmer MSA fasta using SignalP processed fastas




# need to make a data-table with the following columns:
# 1. ensemble_id
# 2. fasta_path
# 3. refdb_path (fasta from msa)
# 4. method (needle or water)

workflow_meta <- readRDS(glue(
  "{wkdir}/data/interim/tmp/",
  "2023-03-21_workflow-meta-2_JackhmmerQC.rds"
))

workflow_meta %>% glimpse()

prep_batch_input_df <- workflow_meta %>%
  filter(jackhmmer_msa_exists) %>%
  dplyr::rename(
    fasta_path = fasta_path_processed,
    refdb_path = jackhmmer_msa_path
  ) %>%
  mutate(
    output_dir =
      "/central/scratch/jbok/emboss_alignments_immport_self"
      ) %>%
  select(ensembl_id, fasta_path, refdb_path, output_dir)
prep_batch_input_df %>% glimpse()

prep_batch_input_df_filter <-
  bind_rows(
    prep_batch_input_df %>% mutate(method = "needle"),
    prep_batch_input_df %>% mutate(method = "water")
  ) %>%
  mutate(
    output_rds = glue(
      "{output_dir}/{ensembl_id}",
      "_{fs::path_ext_remove(basename(refdb_path))}.{method}.rds"
    ),
    output_rds_exists = map_lgl(output_rds, ~ file.exists(.))
  ) %>% 
  filter(!output_rds_exists) %>%
    select(-c(output_rds, output_rds_exists))

prep_batch_input_df_filter %>% glimpse()

batchtools_params <- 
  prep_batch_input_df_filter %>%
  arrange(ensembl_id)

batchtools_params %>% glimpse()
# batchtools_params %>% View()


# configure registry ----
cluster_run <- glue("{Sys.Date()}_EMBOSS-SELF_ID-{rand_string()}/")
message("\n\nRUNNING:  ", cluster_run, "\n")
breg <- makeRegistry(
  file.dir = glue(
    "{wkdir}/.cluster_runs/",
    cluster_run
  ),
  seed = 42
)
breg$cluster.functions <- batchtools::makeClusterFunctionsSlurm(
  template = glue("{wkdir}/batchtools_templates/batchtools.slurm.tmpl"),
  scheduler.latency = 0.1,
  fs.latency = 65
)
# Submit Jobs ----
jobs <- batchMap(
  fun = align_fasta_sequences,
  args = batchtools_params,
  reg = breg
)
setJobNames(jobs,
  paste("EMBOSS",
  batchtools_params$ensembl_id,
  batchtools_params$method,
  sep = "_"),
  reg = breg
)
getJobNames(jobs, reg = breg)
submitJobs(jobs,
  resources = list(
    walltime = 172800,
    memory = 50000,
    ncpus = 1,
    max.concurrent.jobs = 9999
  )
)






aur <- glue(
  "MRKFSRYAFTSMATVTLLSSLTPAALASDTNHKPATSDINFEITQKSDAVKALKELPKSE",
  "NVKNHYQDYSVTDVKTDKKGFTHYTLQPSVDGVHAPDKEVKVHADKSGKVVLINGDTDAK",
  "KVKPTNKVTLSKDEAADKAFNAVKIDKNKAKNLQDDVIKENKVEIDGDSNKYIYNIELIT",
  "VTPEISHWKVKIDADTGAVVEKTNLVKEAAATGTGKGVLGDTKDININSIDGGFSLEDLT",
  "HQGKLSAYNFNDQTGQATLITNEDENFVKDDQRAGVDANYYAKQTYDYYKNTFGRESYDN",
  "HGSPIVSLTHVNHYGGQDNRNNAAWIGDKMIYGDGDGRTFTNLSGANDVVAHELTHGVTQ",
  "ETANLEYKDQSGALNESFSDVFGYFVDDEDFLMGEDVYTPGKEGDALRSMSNPEQFGQPS",
  "HMKDYVYTEKDNGGVHTNSGIPNKAAYNVIQAIGKSKSEQIYYRALTEYLTSNSNFKDCK",
  "DALYQAAKDLYDEQTAEQVYEAWNEVGVE"
)


tst <- blastr(aur)
tst_nr <- blastr(aur, db = "nr")




# workflow_meta <- readRDS(
#   glue(
#     "{wkdir}/data/interim/tmp/",
#     "2023-03-22_workflow-meta-3_HmmersearchQC.rds"
#   )
# )

# target_df <- workflow_meta %>%
#   filter(gene_name == "HSP90AA1")


# tst <- blastr(
#   sequences = target_df$sequences_aa,
#   ensembl_id = target_df$ensembl_id,
#   output_limit = 10,
#   evalue = 1
# )

# tst %>% glimpse
# tst %>% View



# library(tidyverse)
# library(ggtree)

# tree <- read.tree(
#   glue("{wkdir}/data/input/banfield-tree-of-life-2016/banfield_tree_MLE-concat-rRNA-newick.tree.txt")
# )

# ggtree(tree) + 
#   theme_tree2()

# ggtree(tree) + geom_text(aes(label=node), hjust=-.3)
# ggtree(tree, layout="circular")

# # Finally, add tip labels and adjust axis
# ggtree(tree) + 
#   theme_tree2() + 
#   geom_tiplab(align=TRUE, linesize=.5)





# # configure registry ----
# breg <- makeRegistry(
#   file.dir = glue(
#     "{wkdir}/.cluster_runs/",
#     "{Sys.Date()}_BLAST_{sample(10000:99999, 1)}/"
#   ),
#   seed = 42
# )
# breg$cluster.functions <- batchtools::makeClusterFunctionsSlurm(
#   template = glue("{wkdir}/batchtools_templates/batchtools.slurm.tmpl"),
#   scheduler.latency = 2,
#   fs.latency = 65
# )
# # Submit Jobs ----
# jobs <- batchMap(
#   fun = blastr,
#   args = dplyr::select(keyname_hits, sequences_aa, ensembl_id),
#   reg = breg
# )
# setJobNames(jobs, paste0("BLAST_", keyname_hits$ensembl_id), reg = breg)
# submitJobs(jobs,
#   resources = list(walltime = 10800,
#     memory = 1024,
#     ncpus = 8,
#     max.concurrent.jobs = 9999)
)

# getJobNames()$job.id
# waitForJobs(sleep = 5)
# tst <- loadResult(1)
# tst %>% glimpse

# proteins_of_interest <- c("interleukin", "interferons", "TGFB", "TNFA")
# keyname_search <- gget$search(proteins_of_interest, "homo_sapiens")
# keyname_hits <- keyname_search %>% filter(
#   biotype == "protein_coding",
#   !grepl("receptor", ensembl_description)
# )
# # collect amino acid sequences for each protein using the ensembl_id
# gget_seq <- function(ensembl_id, amino_acid = TRUE) {
#   seq_results <- gget$seq(ensembl_id, translate = amino_acid)
#   return(seq_results[2])
# }
# keyname_hits %<>% mutate(sequences_aa =  map_chr(ensembl_id, gget_seq))
# keyname_hits %>% glimpse
# saveRDS(keyname_hits, glue("{wkdir}/data/interim/tmp/",
# "{Sys.Date()}_gget_proteins-of-interest.rds")
# )






# Downloading MGnify Database
# http://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2023_02/

# wget_download_slurm(
#   jobname = "mgy_clusters_fasta",
#   slurm_out = "/central/scratch/jbok/slurmdump/",
#   walltime = "4-00:00:00",
#   download_link = "http://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2023_02/mgy_clusters.fa.gz",
#   output_dir = "/central/groups/MazmanianLab/joeB/Downloads/RefDBs/MGnify_2023_02/"
# )

# for (n in 1:30){
#   jname <- glue("mgy_proteins_{n}")
#   dlink <- glue(
#     "http://ftp.ebi.ac.uk/pub/databases/metagenomics/",
#     "peptide_database/2023_02/mgy_proteins_{n}.fa.gz"
#   )
#   wget_download_slurm(
#     jobname = jname,
#     slurm_out = "/central/scratch/jbok/slurmdump/",
#     walltime = "20:00:00",
#     download_link = dlink,
#     output_dir = "/central/groups/MazmanianLab/joeB/Downloads/RefDBs/MGnify_2023_02/"
#   )
# }



# formatting_params <- "query,target,fident,alnlen,mismatch,qstart,qend,tstart,tend,evalue,bits,prob,lddt,alntmscore,taxname"




# max_prob_df <- esm_results_df %>%
#   group_by(id) %>%
#   slice_max(order_by = prob, n = 1) %>%
#   slice_max(order_by = evalue, n = 1)

# max_prob_df %>% glimpse
# # View(max_prob_df)

# max_prob_df %>%
#   filter(prob > 0.5) %>%
#   ggplot(aes(prob, fct_reorder(gene_target_id, prob))) +
#   geom_point(aes(color = alntmscore)) +
#   theme_bw() +
#     scale_color_viridis_c(option = "H") +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1)
#     )

  # esm_results_df %>% 
  # ggplot(aes(alntmscore, lddt)) +
  # geom_point(aes(color = prob), size = 0.3, alpha = 0.7) +
  # theme_bw() +
  # scale_color_viridis_c(option = "H") +
  # geom_smooth(method = "loess", se = FALSE) +
  # coord_fixed() +
  # # labs(x = "Probability", y = "E-value", color = "Alignment score") +
  # theme(
  #   legend.position = "top",
  #   axis.text.x = element_text(angle = 45, hjust = 1)
  # )



# output_name <- "IL10_formatted"
# input_pdb <- glue("{wkdir}/il10.pdb")
# for (DB in input_db_list){
#   db_name <- basename(DB)
#   pdb_name <- fs::path_ext_remove(basename(input_pdb))
#   output_name <- glue("{pdb_name}_{db_name}")
#   print(db_name)
#   foldseek_cmd <- glue(
#     "conda run -n foldseek",
#     " foldseek easy-search",
#     " {input_pdb} {DB} {output_name}",
#     " tmp --threads {threads}",
#     " --format-output {formatting_params}" #,
#     # " --taxon-list 2,4751,2157,10239" # bacteria, fungi, archaea and viruses
#   )
#   slurm_shell_do(
#     cmd = foldseek_cmd, memory = "100G", walltime = 14400, ncpus = threads
#   )
# }


colnames(results_df) <- unlist(str_split(formatting_params, ","))
View(results_df)


library(ggmsa)

ggmsa(
  msa = glue("{wkdir}/data/interim/alignment_hits/2023-03-23_NFKBIE_top-10-per-subcatalog.fasta"),
  char_width = 0.5, seq_name = TRUE
  ) +
  geom_seqlogo() +
  geom_msaBar()

# # # convert hmmer stockholm to fasta
# # esl-reformat fasta msacopy.fasta > aln.fasta

# p_msa <-
#   ggmsa(
#   "/central/scratch/jbok/mim_temp/jackhmmer_results/immport/ENSG00000034510/aln.a2m",
#   # start = 408, end = 500,
#   # char_width = 0.5, 
#   seq_name = TRUE) +
#   scale_y_discrete(labels = msa_seq_labels) +
#   geom_seqlogo() +
#   geom_msaBar()

# ggsave(
#   glue("{analysis_dir}/figures/test_msa.png"),
#   plot = p_msa, width = 20, height = 10, dpi = 600
#   )


#   library(inferrnal)
#   msa <- read_stockholm_msa(
#     "/central/scratch/jbok/mim_temp/jackhmmer_results/immport/ENSG00000034510/ENSG00000034510_msa.fasta",
#     type = "AA"
#   )

# msa %>% str
















# open the file connection
con <- file(
    glue(
    "/central/scratch/jbok/mim_temp/clustering/",
    "2023-04-17_MMSeq2-90_MGnify90_custom90_cluster.tsv"),
    open = "r"
    )

# initialize an empty list to store the selected representative sequences
selected <- list()
# initialize a list to store the first column values
values <- list()
blacklist <- list()
nrow <- 0
# loop over each line in the file
while (length(line <- readLines(con, n = 1)) > 0) {
  message(nrow)
  nrow <- nrow + 1
  # split the line into columns
  cols <- strsplit(line, "\t")[[1]]
  # check if the value has already appeared in the first column
  if (cols[[1]] %in% blacklist) {
    next
  }
  else if (cols[[1]] %in% values) {
    # remove the value from the list of values
    values <- values[values != cols[[1]]]
    blacklist <- c(blacklist, cols[[1]])
  } else {
    # add the value to the list of values
    values <- c(values, cols[[1]])
  }
}
# close the file connection
close(con)
# extract the selected representative sequences
message("Length Unique: ", length(values))
message("Length Blacklist: ", length(blacklist))
selected <- values
saveRDS(selected, "MGnify90_custom90_singleton-repseqs.rds")



# open the file connection
con <- file(
    glue(
    "/central/scratch/jbok/mim_temp/clustering/",
    "2023-04-17_MMSeq2-90_UniRef90_custom90_cluster.tsv"),
    open = "r"
    )

# initialize an empty list to store the selected representative sequences
selected <- list()
# initialize a list to store the first column values
values <- list()
blacklist <- list()
nrow <- 0
# loop over each line in the file
while (length(line <- readLines(con, n = 1)) > 0) {
  message(nrow)
  nrow <- nrow + 1
  # split the line into columns
  cols <- strsplit(line, "\t")[[1]]
  # check if the value has already appeared in the first column
  if (cols[[1]] %in% blacklist) {
    next
  }
  else if (cols[[1]] %in% values) {
    # remove the value from the list of values
    values <- values[values != cols[[1]]]
    blacklist <- c(blacklist, cols[[1]])
  } else {
    # add the value to the list of values
    values <- c(values, cols[[1]])
  }
}
# close the file connection
close(con)
# extract the selected representative sequences
message("Length Unique: ", length(values))
message("Length Blacklist: ", length(blacklist))
selected <- values
saveRDS(selected, "UniRef90_custom90_singleton-repseqs.rds")





























blastp_res_plot <- readRDS(
  glue(
    "{wkdir}/data/interim/phylogenetic_analysis/",
    "2023-06-04_blastp_res_trimmed.rds"
  )
)
blastp_res_plot %>% glimpse()


# set.seed(42)
# rand_eids <- blastp_res_plot$qaccver %>%
#   unique() %>%
#   sample(20)

# blastp_res_plot_wk <- blastp_res_plot %>%
#   filter(qaccver %in% rand_eids)
# blastp_res_plot_wk %>% glimpse

# blastp_res_wide <- blastp_res_plot %>%
#   select(qaccver, genome, evalue) %>%
#   mutate(evalue = -log10(evalue + 1)) %>% 
#   pivot_wider(names_from = qaccver, values_from = evalue, values_fill = 0) %>%
#   column_to_rownames(var = "genome")
# blastp_res_wide %>% glimpse
# eid_order <-
#     blastp_res_wide %>%
#     t() %>%
#     dist(method = "euclidean") %>%
#     seriate() %>%
#     get_order()
# eid_order
# genome_order <-
#     blastp_res_wide %>%
#     dist(method = "euclidean") %>%
#     seriate() %>%
#     get_order()
# genome_order
# saveRDS(
#   colnames(blastp_res_wide)[eid_order],
#   glue(
#     "{wkdir}/data/interim/phylogenetic_analysis/",
#     "{Sys.Date()}_seriated-eid-order.rds"
#   )
# )
# saveRDS(
#   rownames(blastp_res_wide)[genome_order],
#   glue(
#     "{wkdir}/data/interim/phylogenetic_analysis/",
#     "{Sys.Date()}_seriated-genome-order.rds"
#   )
# )


blastp_res_plot %>%
  ggplot(aes(x = root, y = organism_name)) +
  geom_tile(aes(fill = kingdom)) +
  scale_fill_brewer(palette = "Set3") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank()
  )

 colorbar <-
   merged_meta %>%
   mutate(root = ".") %>% 
    ggplot(aes(x = root, y = genome, fill = kingdom)) +
    geom_tile() +
    labs(fill = NULL) +
    scale_y_discrete(position = "right") +
    theme_minimal() +
    # scale_fill_manual(values = comm_humz_cols) +
    theme(
      axis.text.x = element_blank(),
      axis.title = element_blank(),
      axis.text.y = element_text(size = 1)
    )

final_plot <- heatmap_plot +
  colorbar + plot_layout(widths = c(40, 1), guides = "collect") & theme(legend.position = "top")











# ![](../figures/phylogenomics/kingdom_piechart.png){#fig-kingdom-piechart width=50%}


# loop through data.frames and select top hit per organism, bind to rest of list

# dt <- blastp_res
# library(stringi)
# # Assuming your data table is named 'dt' and the column containing the cell values is named 'cell_value'

# # Define a function to extract the content within the last set of brackets
# extract_last_bracket_content <- function(saccver) {
#   pattern <- "\\[([^\\[\\]]+)\\](?!.*\\[)"
#   match <- stri_match_last_regex(saccver, pattern)
#   if (!is.na(match) && !is.na(match[,2])) {
#     gsub("_", " ", match[,2])
#   } else {
#     NA
#   }
# }
# possibly_extract_last_bracket_content <-
#   possibly(extract_last_bracket_content, otherwise = NA_character_)

# dt <- blastp_res %>%
#   ungroup() %>%
#   slice_sample(n = 10000) %>%
#   mutate(new_column = map_chr(saccver, possibly_extract_last_bracket_content)) %>%
#   pull(new_column) %>%
#   unique()

# setdiff(dt, tax_metadata$organism_name)

# # extract headers from fasta file
# cd /central/groups/MazmanianLab/joeB/Microbiota-Immunomodulation/Microbiota-Immunomodulation/data/interim/refseq_proteins/genus_sampling
# awk 'sub(/^>/, "")' RefSeq-Genus-sample_2023-05-18_numbered.faa > RefSeq-Genus-sample_2023-05-18_numbered-header.txt

# # chunk header into 50 files
# mkdir -p header_chunks
# split -l 3000000 RefSeq-Genus-sample_2023-05-18_numbered-header.txt header_chunks/chunk_




# # read in header txt file as a data.frame
# library(data.table)
# genus_dir <- glue("{wkdir}/data/interim/refseq_proteins/genus_sampling")
# db_headers <- read.delim(
#   glue("{genus_dir}/RefSeq-Genus-sample_2023-05-18_numbered-header.txt"),
#   header = FALSE, stringsAsFactors = FALSE
# ) %>%
#   as.data.table()

# t %>% glimpse










# # rezip genomes that failed to rezip for some reason
# future::plan(multisession, workers = 12)
# batch_input_df$refseq_fna_path %>%
#   keep(!grepl("\\.gz$", .) & grepl("\\.fna", .)) %>%
#   furrr::future_map( ~ shell_do_zip(.))

# future::plan(multisession, workers = 24)
# list.files(
#   glue("{genus_dir}/download/genomic_fna"),
#   full.names = TRUE,
#   pattern = ".fna"
#   ) %>%
#   keep(!grepl("\\.gz$", .) & grepl("\\.fna", .)) %>%
#   furrr::future_map( ~ shell_do_zip(.))


# get list of fna file that already have compressed files
genomes_dir <- 
  glue("{scratch_dir}/refseq_proteins/genus_sampling/download/genomic_fna")
fna_paths <- list.files(
  genomes_dir,
  full.names = TRUE,
  pattern = ".fna"
  ) %>%
  keep(!grepl("\\.gz$", .) & grepl("\\.fna", .)) %>%
    paste0(., ".gz")

gzip_paths <- list.files(
    genomes_dir,
    full.names = TRUE) %>%
    keep(grepl("\\.gz$", .)) #%>% gsub("\\.gz$", "", .)

# # remove decompressed files that also have a compressed version
# purrr::map(
#   intersect(gzip_paths, fna_paths),
#   ~ shell_do(glue("rm -rf {.}"))
# )

# # collecting protein path availability info
# future::plan(multisession, workers = 15)
# batch_input_df_euk <-
#   batch_input_df %>%
#   filter(kingdom %nin% c("k__Bacteria", "k__Archaea")) %>%
#   mutate(
#     local_faa_path = purrr::map_chr(
#       url,
#       ~ glue(
#         "{wkdir}/data/interim/refseq_proteins/",
#         "genus_sampling/download/protein_faa/",
#         "{strex::str_after_last(., '/')}_protein.faa.gz"
#       )
#     ),
#     local_faa_path_exists = purrr::map_lgl(
#       local_faa_path, file.exists
#     ),
#     web_faa_path = purrr::map_chr(
#       url,
#       ~ glue("{.}/{strex::str_after_last(., '/')}_protein.faa.gz")
#     ),
#     web_faa_path_exists = furrr::future_map_lgl(
#       web_faa_path, RCurl::url.exists
#     )
#   )


# saveRDS(
#   batch_input_df_euk,
#   glue("{genus_dir}/Eukaryota-Metazoa-Fungi-Viridiplantae_fastainfo.rds")
# )

# genus_dir <-
#   "/central/groups/MazmanianLab/joeB/Microbiota-Immunomodulation/Microbiota-Immunomodulation/data/interim/refseq_genomes/genus_sampling"
# batch_input_df_euk <- readRDS(
#   glue("{genus_dir}/Eukaryota-Metazoa-Fungi-Viridiplantae_fastainfo.rds")
# )
# batch_input_df_euk %>% glimpse


# # Piechart of kingdom distribution
# kingdom_summary_stats <-
#   batch_input_df %>%
#   select(kingdom) %>%
#   group_by(kingdom) %>%
#   summarize(value = n()) %>%
#   mutate(kingdom = gsub("k__", "", kingdom))

# # Compute the position of labels
# kingdom_summary_stats %<>%
#   arrange(desc(value)) %>%
#   mutate(prop = value / sum(kingdom_summary_stats$value) *100) %>%
#   mutate(ypos = cumsum(prop)- 0.5*prop )

# p_pie <- ggplot(kingdom_summary_stats, aes(x="", y=prop, fill=kingdom)) +
#   geom_bar(stat="identity", width=1, color="white") +
#   coord_polar("y", start=0) +
#   theme_void() +
#   scale_fill_brewer(palette="Set1")

# ggsave(
#   filename = glue("{wkdir}/figures/phylogenomics/kingdom_piechart.png"),
#   plot = p_pie,
#   width = 5,
#   height = 5
# )




# # cleanup
# fna_genomes <- list.files(
#     glue("/central/scratch/jbok/genomic_fna"),
#     pattern = ".fna", full.names = TRUE
# )
# fasta_prod <- list.files(
#   glue("/central/scratch/jbok/prodigal_proteins"),
#   pattern = ".fasta", full.names = TRUE
# )

# fna_genomes %>% length
# fasta_prod %>% length

# remove decompressed files that also have a compressed version
# purrr::map(
#   intersect(fasta_prod, fna_prod),
#   ~ shell_do(glue("rm -rf {.}"))
# )
# purrr::map(
#   setdiff(fna_prod, fasta_prod),
#   ~ shell_do(glue("mv {.} {gsub('_genomic.fna', '_protein.fasta', .)}"))
# )





#  list.files(
#   glue("/central/scratch/jbok/prodigal_proteins"),
#   pattern = ".fasta", full.names = TRUE
# )  %>% length


# _genomic.fna
# _protein.fasta






























# esm_results_df <-
#   foldkseek_results_paths %>%
#   keep(grepl("_foldseek_ESMAtlas30.m8", .)) %>%
#   purrr::set_names() %>%
#   furrr::future_map(
#     ~ possibly_read_delim(., sep = "\t", header = FALSE)) %>%
#   bind_rows(.id = "id")
# colnames(esm_results_df) <-
#   unlist(strsplit(glue("id,", output_params_base), ","))
# esm_results_df %<>%
#   as.data.table() %>%
#   mutate(gene_target_id =
#     strex::str_before_nth(id, "/", -2) %>%
#     strex::str_after_last("/"))
# saveRDS(
#   esm_results_df,
#   glue(
#     "{wkdir}/data/interim/foldseek_results/",
#     "{Sys.Date()}_ESMAtlas30_foldseek_results_df.rds"
#   )
# )













# p <- 
# foldseek_meta_res %>%
#   filter(gene_target_id %in% top_bits_overall) %>%
#   # group_by(gene_target_id, foldseek_DB) %>%
#   # slice_min(order_by = evalue, n = 20, with_ties = FALSE) %>%
#   ggplot(aes(x = bits, y = fct_reorder(gene_target_id, bits))) +
#   facet_grid(
#     transdomain_blastp_signal~foldseek_DB, 
#     space = "free_y", scales = "free_y"
#   ) +
#   geom_point(aes(color = lddt), position = position_jitter(height = 0.2)) +
#   scale_color_viridis_c(option = "F", direction = -1) +
#   scale_x_log10() +
#   theme_bw()

# ggsave(
#   glue("{wkdir}/figures/foldseek/{Sys.Date()}_top-gene-targets.png"),
#   p,
#   width = 12,
#   height = 12
# )



foldseek_meta_res %>% glimpse 



# foldseek_meta_res %>% glimpse
# foldseek_meta_res %>%
#   filter(gene_target_id == "C5_ENSG00000106804") %>%
#   # filter(foldseek_DB == "PDB") %>%
#   arrange(desc(bits)) %>%
#   View

# foldseek_meta_res %>%
#   filter(gene_target_id == "FSHB_ENSG00000131808") %>%
#   filter(foldseek_DB == "Alphafold UniProt50") %>%
#   arrange(desc(bits)) %>%
#   View

# foldseek_meta_res %>%
#   filter(gene_target_id == "C5_ENSG00000106804") %>%
#   # filter(foldseek_DB == "PDB") %>%
#   arrange(desc(bits)) %>%
#   View

# foldseek_meta_res %>%
#   filter(gene_target_id == "IL10_ENSG00000136634") %>%
#   arrange(desc(bits)) %>%
#   View



# p_esm_dist <-
foldseek_meta_res %>%
  filter(foldseek_DB == "ESMAtlas30") %>%
  mutate(signal6p_Prediction = factor(
    signal6p_Prediction,
    levels = c("OTHER", "SP")
  )) %>%
    mutate(local_global_delta = lddt - alntmscore) %>%
    mutate(poi = if_else(
      local_global_delta > 0.3 & 
      lddt > 0.7 & 
      alntmscore < 0.8 & alntmscore > 0.2, TRUE, FALSE
      )
    ) %>% 
  # filter(lddt > 0.7 & alntmscore < 0.7) %>%
    slice_sample(n = 200000) %>%
    ggplot(aes(x = alntmscore, y = lddt)) +
    geom_point(alpha = 1, aes(color = signal6p_Prediction)) +
    coord_fixed() +
    theme_bw() +
    scale_color_d3() +
    theme(legend.position = "bottom") +
    facet_zoom(xy = poi == TRUE)

# foldseek_meta_res %>%
#   filter(lddt > 0.7 & alntmscore < 0.7 & alntmscore > 0.2) %>%
#   filter(foldseek_DB == "ESMAtlas30") %>%
#   mutate(local_global_delta = lddt - alntmscore) %>%
#   ggplot(aes(x = local_global_delta, y = lddt, color = gene_category)) +
#   geom_point()

# foldseek_meta_res %>%
#   filter(lddt > 0.7 & alntmscore < 0.7 & alntmscore > 0.2) %>%
#   filter(foldseek_DB == "ESMAtlas30") %>%
#   mutate(local_global_delta = lddt - alntmscore) %>%
#   ggplot(
#     aes(
#       x = local_global_delta,
#       y = fct_reorder(gene_category, local_global_delta)
#     )) +
#   geom_point(
#     aes(color = signal6p_Prediction),
#     position = position_jitterdodge()
#   )

#' it looks like antimicrobials, cytokines/chemokines + their receptors, and antigen procressing and presentation
#' have some of the largest deltas between local and global scores, as well as the a large density of points


# foldseek_meta_res %>%
#   filter(lddt > 0.7 & alntmscore < 0.7 & alntmscore > 0.2) %>%
#   filter(foldseek_DB == "ESMAtlas30") %>%
#   filter(grepl("Cytokine|Chemokine", gene_category)) %>%
#   filter(!grepl("Receptor", gene_category)) %>%
#   mutate(local_global_delta = lddt - alntmscore) %>%
#   pull(gene_target_id) %>% unique

# p_esm_dist <- foldseek_meta_res %>%
#   filter(foldseek_DB == "ESMAtlas30") %>%
#   filter(lddt > 0.7 & alntmscore < 0.7) %>%
#   # slice_sample(n = 100000) %>%
#   ggplot(aes(x = alntmscore, y = lddt, color = gene_category)) +
#   geom_point() +
#   coord_fixed() +
#   theme_bw() +
#   theme(legend.position = "bottom")

# ggplotly(p_esm_dist)

# ggsave(
#   glue("{wkdir}/figures/foldseek/{Sys.Date()}_ESMAtlas30_alntmscore-lddt.png"),
#   p_esm_dist,
#   width = 30, height = 20
# )




top_hits_list <- foldseek_meta_res_summ %>%
  left_join(meta) %>%
  filter(grepl("Cytokine", gene_category)) %>%
  filter(grepl("IL", gene_target_id)) %>%
  filter(!grepl("Receptor", gene_category)) %>%
  arrange(desc(m_bits))

# make list of genes of interest
cytokine_gois <- c(
  "TGFA_ENSG00000163235",
  top_hits_list$gene_target_id[1:10]
)

foldseek_meta_wide <- foldseek_meta_res %>%
  filter(prob > 0.5) %>%
  group_by(gene_target_id, taxname) %>%
  summarize(eid_hits = n()) %>%
  mutate(eid_hits = 1) %>%
  pivot_wider(
    names_from = gene_target_id,
    values_from = eid_hits, values_fill = 0
  ) %>%
  column_to_rownames("taxname")





# eids <- top_bits_overall %>%
#   strex::str_after_first(., "_")






# # collect amino acid sequences for each protein using the ensembl_id
# gget_seq <- function(ensembl_id, amino_acid = TRUE) {
#   seq_results <- gget$seq(ensembl_id, translate = amino_acid)
#   return(seq_results[2])
# }
# keyname_hits %<>% mutate(sequences_aa =  map_chr(ensembl_id, gget_seq))
# keyname_hits %>% glimpse
# saveRDS(keyname_hits, glue("{wkdir}/data/interim/tmp/",
# "{Sys.Date()}_gget_proteins-of-interest.rds")
# )




# taxcat_df <- foldseek_res_df %>%
#   left_join(meta, by = "gene_target_id") %>%
#   group_by(gene_category, taxname) %>%
#   summarize(taxa_counts = n()) %>%
#   mutate(gene_category = fct_reorder(gene_category, taxa_counts))

# taxcat_wide_df <- taxcat_df %>%
#   pivot_wider(
#     names_from = gene_category,
#     values_from = taxa_counts, values_fill = 0
#   ) %>%
#   column_to_rownames("taxname")

# # Seriating Heatmap
# cat_order <-
#   taxcat_wide_df %>%
#   t() %>%
#   dist(method = "euclidean") %>%
#   seriate(method = "OLO") %>%
#   get_order() %>%
#   sort()
# cat_order
# tax_order <-
#   taxcat_wide_df %>%
#   dist(method = "euclidean") %>%
#   seriate(method = "OLO") %>%
#   get_order() %>%
#   sort()
# tax_order

# p_heatmap <- taxcat_df %>%
#   mutate(taxname = factor(taxname,
#     levels = names(tax_order), ordered = TRUE)) %>%
#   mutate(gene_category = factor(gene_category,
#     levels = names(cat_order), ordered = TRUE)
#     ) %>%
#   ggplot(aes(y = gene_category, x = taxname, fill = log2(taxa_counts + 1))) +
#     geom_tile() +
#     theme_bw() +
#     scale_fill_viridis_c(option = "G", direction = -1) +
#     labs(y = "Immune Gene Category", x = "Taxa",
#     fill = expression(log[2] ~ "[Taxa Hits + 1]")) +
#     theme(
#       legend.position = "top",
#       axis.text.x = element_blank(),
#       axis.text.y = element_text(size = 8),
#     )
# p_heatmap
#
# ggsave(
#   glue("{wkdir}/figures/foldseek/{Sys.Date()}_taxa-gene-category_heatmap.png"),
#   p_heatmap,
#   width = 18,
#   height = 7
# )




# /central/groups/MazmanianLab/joeB/Microbiota-Immunomodulation/Microbiota-Immunomodulation/data/interim/fastas/processed/monomers/immport/ENSG00000170345/prediction_results.txt


# foldseek_meta_res[sample(1:nrow(foldseek_meta_res), 10000), ] %>%
#   ggplot(aes(x=lddt, y=evalue)) +
#   scale_y_log10() +
#   geom_point()

# foldseek_meta_res[sample(1:nrow(foldseek_meta_res), 10000), ] %>%
#   ggplot(aes(x=lddt, y=prob)) +
#   geom_point()

# foldseek_meta_res %>%
#   filter(prob >= 0.95) %>%
#   filter(evalue < 0.05) %>%
#   slice_sample(n = 10000) %>%
#   ggplot(aes(x=lddt, y=alntmscore,
#     size = log10(bits),
#     fill = -log10(evalue)
#   )) +
#   scale_fill_viridis_c(option = "F", direction = -1) +
#   geom_point(alpha = 0.6, shape = 21)

# foldseek_bin_df <- foldseek_res_df %>%
#   group_by(gene_target_id, taxname) %>%
#   summarize(eid_hits = n()) %>%
#   mutate(eid_hits = 1) %>%
#   dplyr::left_join(meta, by = "gene_target_id") %>%
#   left_join(foldseek_meta_res) %>%
#   mutate(gene_target_id = fct_reorder(gene_target_id, -eid_hits)) %>%
#   mutate(taxname = fct_reorder(taxname, -eid_hits))

# top_heatmap_plot <- foldseek_bin_df %>%
#   filter(prob > 0.5) %>% 
#   filter(gene_target_id %in% cytokine_gois) %>%
#   ggplot(aes(x = gene_target_id, y = taxname, fill = bits)) +
#   # facet_grid(cols = vars(transdomain_blastp_signal),
#   #   scales = "free_x", space = "free") +
#   geom_tile() +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4))
# top_heatmap_plot

# ggsave(
#   glue("{wkdir}/figures/foldseek/{Sys.Date()}_top-gene-targets_heatmap.png"),
#   width = 15,
#   height = 62,
#   limitsize = FALSE
# )




# # Seriating Heatmap
# tax_order <-
#   foldseek_meta_wide %>%
#   t() %>%
#   dist(method = "binary") %>%
#   seriate(method = "OLO") %>%
#   get_order() %>%
#   sort()
# tax_order

# eid_order <-
#   foldseek_meta_wide %>%
#   dist(method = "binary") %>%
#   seriate(method = "OLO") %>%
#   get_order() %>%
#   sort()
# eid_order


# foldseek_bin_df <- foldseek_res_df %>%
#   filter(prob > 0.5) %>%
#   group_by(gene_target_id, taxname) %>%
#   summarize(eid_hits = n()) %>%
#   mutate(eid_hits = 1) %>%
#   dplyr::left_join(meta, by = "gene_target_id") %>%
#   mutate(gene_target_id = fct_reorder(gene_target_id, -eid_hits)) %>%
#   mutate(taxname = fct_reorder(taxname, -eid_hits))

# foldseek_bin_df %>%
#   ggplot(aes(x = gene_target_id, y = taxname, fill = eid_hits)) +
#   # facet_grid(cols = vars(transdomain_blastp_signal),
#   #   scales = "free_x", space = "free") +
#   geom_tile() +
#   theme_bw()





html_paths <- list.files(
  glue(
    "{wkdir}/data/interim/foldseek_results/html_output"
  ),
  recursive = TRUE,
  full.names = TRUE,
  pattern = ".html"
)
html_paths %>% length()

# pdb_paths <- list.files(
#   glue(
#     "{wkdir}/data/interim/foldseek_results/pdb_output"
#   ),
#   recursive = TRUE,
#   full.names = TRUE,
#   pattern = ".pdb"
# )










# pdb_df <- data.table("pdb_path" = pdb_paths) %>%
#   mutate(
#     gene_target_id =
#       str_before_nth(pdb_path, "/", -2) %>% str_after_last(., "/"),
#     filename = str_after_nth(pdb_path, "/", -2) %>%
#       str_after_last(., "/") %>%
#       str_after_first("_"),
#     db = case_when(
#       grepl("Alphafold_Proteome", filename) ~
#         "Alphafold_Proteome",
#       grepl("foldseek_Alphafold_SwissProt", filename) ~
#         "foldseek_Alphafold_SwissProt",
#       grepl("foldseek_PDB", filename) ~
#         "foldseek_PDB",
#       grepl("foldseek_ESMAtlas30", filename) ~
#         "
#         ",
#       grepl("foldseek_Alphafold_UniProt50", filename) ~
#         "foldseek_Alphafold_UniProt50",
#       TRUE ~ "ERROR"
#       ),
#     res_pair = strex::str_after_first(filename, db),
#     target = case_when(
#       grepl("MGYP", res_pair) ~ str_after_last(res_pair, "_"),
#       grepl("AF-", res_pair) ~ str_after_nth(res_pair, "_", -2),
#       db == "foldseek_PDB" ~ str_after_nth(res_pair, "_", -2),
#       TRUE ~ "ERROR"
#     ),
#     query = str_remove(res_pair, target),
#     query = case_when(
#       grepl("MODEL", query) ~ gsub("_MODEL_\\d+", "", query),
#       TRUE ~ query
#     ),
#     query = str_before_last(query, "_"),
#     target = str_remove(target, "\\.pdb")
#   )

# pdb_df


# pdb_df %>%
#   slice_sample(n = 50) %>%
#   pull(target)
# pdb_df %>%
#   slice_sample(n = 50) %>%
#   pull(query)





# "SERPINA3_ENSG00000196136"

# Homone regulation
# NR3C1


# Metabolism ---
# "INS_ENSG00000254647"
# "KL_ENSG00000133116"
# "ADM2_ENSG00000128165"
# ADIPOR1

# Immune signaling
# Complement  ===

# Receptors  ===
# "TLR8_ENSG00000101916"

# Chemokines  ===

# Cytokines  ===
# "IL10_ENSG00000136634"
# "TGFB3_ENSG00000119699"




# foldseek_meta_res$gene_target_id %>%
#   unique() %>%
#   keep(grepl("ADM", .))





# # 2. Select target gene from results table and get PDBs of interest to plot
# fsk_hit_df_all <-
#   foldseek_meta_res %>%
#   mutate(
#     query_pdb = gsub("_[A-Za-z]$", "", query)
#   ) %>%
#   filter(gene_target_id == "IL34_ENSG00000157368")

# # getting the top 3 query structures per target
# fsk_query_hits_per_DB_df <- fsk_hit_df_all %>%
#   group_by(query_pdb, foldseek_DB) %>%
#   slice_max(order_by = lddt, n = 1, with_ties = FALSE) %>%
#   group_by(foldseek_DB) %>%
#   slice_max(order_by = lddt, n = 3, with_ties = FALSE) %>%
#   select(foldseek_DB, query)

# db_list <- fsk_query_hits_per_DB_df$foldseek_DB %>% unique()
# fsk_query_hits_per_DB <- db_list %>%
#   purrr::set_names() %>%
#   purrr::map(
#     ~ filter(fsk_query_hits_per_DB_df, foldseek_DB == .x) %>%
#       pull(query) %>%
#       unique()
#   )

#  top_query_hits <- unlist(fsk_query_hits_per_DB) %>%
#    unname() %>%
#    unique()

# # Getting the top 3 hits per query hit per DB (up to 3x3x5 = 45 rows)
# hits_to_plot_res <- db_list %>%
#   purrr::map_dfr(
#     ~ filter(fsk_hit_df_all, foldseek_DB == .x) %>%
#       filter(query %in% fsk_query_hits_per_DB[[.x]]) %>%
#       group_by(query, foldseek_DB) %>%
#       slice_max(order_by = lddt, n = 3, with_ties = FALSE)
#   )

# # 3. Adding in query pdb paths
# # 4. Loading in aligned foldseek result PDBs
# hits_to_plot <- hits_to_plot_res %>%
#   mutate(
#     pdb_output_path =
#       gsub("2023-07-22_easysearch", "pdb_output", id) %>%
#       gsub(".m8", "", .),
#     target_pdb_path =
#       glue("{pdb_output_path}", "{query}_{target}.pdb"),
#     target_pdb_path_exists = file.exists(target_pdb_path),
#     query = toupper(query)
#   ) %>%
#   left_join(query_pdb_paths, by = join_by(query, gene_target_id, query_pdb))

# # 5. Plotting
# model_results <- list()
# for (db in db_list) {
#   db_filt <- hits_to_plot %>%
#     filter(foldseek_DB == db, target_pdb_path_exists) %>%
#     tidyr::drop_na(query_pdb_path)
#   if (nrow(db_filt) == 0) {
#     next
#   }
#   for (query_path in unique(db_filt$query_pdb_path)) {
#     message("Query: ", basename(query_path))
#     query_filt <- db_filt %>%
#       filter(query_pdb_path == query_path)
#     model_results[[db]][[basename(query_path)]] <-
#       view_foldseek_pdb(
#         query_pdb_path = query_path, # should be a single value
#         query_chain = unique(query_filt$query_chain), # should be a single value
#         target_pdbs_path_list = query_filt$target_pdb_path
#       )
#   }
# }







# mod$TGFB1_ENSG00000105329$PDB %>% names
# mod$TGFB1_ENSG00000105329$PDB$`1KLC.pdb`
# mod$TGFB1_ENSG00000105329$PDB[[3]]


# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# taxadump_dir <-
#   glue(
#     "{wkdir}/data/interim/refseq_genomes/genus_sampling/download"
#   )

# genbank_df <- read.delim(
#   glue("{taxadump_dir}/assembly_summary_genbank.txt"),
#   skip = 1
# )
# genbank_df %>% glimpse

# refseq_df <- read.delim(
#   glue("{taxadump_dir}/assembly_summary_refseq.txt"),
#   skip = 1
# )
# refseq_df %>% glimpse

# library(ggrepel)


# foldseek_meta_res %>%
#   filter(gene_target_id %in% gois) %>%
#   filter(foldseek_DB == "Alphafold UniProt50") %>%
#   filter(fident > 0.8) %>%
#   View

# foldseek_meta_res %>%
#   filter(gene_target_id %in% gois) %>%
#   filter(fident < 0.4, L2 == "Bacteria") %>%
#   View



# eid_query_chains %>% filter(grepl(".pdb_C", query))


eid_query_chains %>% View

foldseek_meta_res %>% glimpse
foldseek_meta_res %>%
  filter

foldseek_res_df_all_temp %>%
  filter(grepl(".pdb_a", query)) %>%
  View

foldseek_meta_res %>%
  filter(gene_target_id == "INS_ENSG00000254647") %>%
  View


foldseek_res_df_all_temp %>%
  mutate(query_chain = strex::str_after_last(query, "._")
  # case_when(
  #   grepl(".pdb_", query) ~ strex::str_after_last(query, ".pdb_"),
  #   TRUE ~ ""
  # )
  ) %>%
  pull(query_chain) %>% 
  unique



rcsb_aln_df <- readRDS(
  glue(
    "{wkdir}/data/interim/pdb_search/",
    "2023-08-07_RCSB-chain-target-global-alignment.rds"
    )
)

rcsb_aln_df %>% glimpse
# rcsb_aln_df %>% View

rcsb_aln_df %>%
  filter(grepl("3H44", pdb_path)) %>%
    group_by(pdb_path) %>%
  slice_max(order_by = global_score, with_ties = FALSE, n = 1) %>%
    View



# setA <- rcsb_aln_df %>%
#   group_by(pdb_path) %>%
#   slice_max(order_by = local_score, with_ties = FALSE, n = 1) %>%
#   mutate(
#     pdb_chain =
#       glue("{strex::str_after_last(pdb_path, '/')}_{chain}")
#   ) %>%
#   pull() %>%
#   unique()

# setB <- rcsb_aln_df %>%
#   group_by(pdb_path) %>%
#   slice_max(order_by = global_score, with_ties = FALSE, n = 1) %>%
#   mutate(
#     pdb_chain =
#       glue("{strex::str_after_last(pdb_path, '/')}_{chain}")
#   ) %>%
#   pull() %>%
#   unique()

# # setA
# # setB
# # setdiff(setA, setB)
# # setdiff(setB, setA)


# uniprot_acc <- "P00520?key=AIzaSyCeurAJz7ZGjPQUtEaerUkBZ3TaBkXrY94"


View(goi_df)

goi_df <- readRDS(
  glue(
    "{wkdir}/data/interim/foldseek_results/",
    "2023-08-10_goi-results-with-paths.rds"
  )
)

goi_df_plot <- goi_df %>% 
  slice_max(order_by = lddt, n = 15, with_ties = FALSE)

# goi_df_plot %>% filter(gene_target_id == "TGFB3_ENSG00000119699") %>% View
# goi_df_plot %>%
#   filter(gene_target_id == "TGFB3_ENSG00000119699") %>%
#   pull(target_pdb_path) %>% unique

# 
# for goi in ...
# "NPPA_ENSG00000175206"
# IL6ST_ENSG00000134352
# goi <- "IL6ST_ENSG00000134352"
# goi <- "CD79B_ENSG00000007312"
# goi <- "CD8B_ENSG00000172116"
# goi <- "CD86_ENSG00000114013"
# goi <- "CD28_ENSG00000178562"
# goi <- "IL12B_ENSG00000113302"
# goi <- "CCL8_ENSG00000108700"
# goi <- "TGFA_ENSG00000163235"
# goi <- "CCL3_ENSG00000278567"
# goi <- "CCL3L1_ENSG00000276085"
# goi <- "CCL3L3_ENSG00000277768"
# goi <- "IGF2_ENSG00000167244"
# goi <- "INS_ENSG00000254647"
# goi <- "CDH1_ENSG00000039068" # something funky going on with pdb_mods$CDH1_ENSG00000039068$`1O6S.pdb_A`$`Alphafold UniProt50`$Bacteria
# goi <- "THBS1_ENSG00000137801"
# goi <- "IL10_ENSG00000136634"
# goi <- "CCL4L2_ENSG00000275313"
# goi <- "CCL3_ENSG00000278567"
# goi <- "CCL4_ENSG00000275302"
# goi <- "CCL5_ENSG00000271503"
# goi <- "CCL7_ENSG00000108688"
# goi <- "CCL8_ENSG00000108700"
# goi <- "CCL11_ENSG00000172156"
# goi <- "CCL18_ENSG00000278167"
# goi <- "CCL19_ENSG00000172724"
# goi <- "CCL20_ENSG00000115009"
# goi <- "ANGPTL1_ENSG00000116194"
goi <- "IL34_ENSG00000157368"

# pb <- progress::progress_bar$new(
#   total = length(gois),
#   format = "[:bar] :percent eta: :eta"
# )

# for (goi in gois) {
#   output_file <- glue(
#     "{wkdir}/data/interim/foldseek_results/",
#     "top_pdb_models/2023-08-07_{goi}_.rds"
#   )
#   if (file.exists(output_file)) next
#   message(glue("Processing {goi}..."))
#   pdb_mods <- list()
#   pb$tick()
#   gti_df <- goi_df_plot %>%
#     filter(gene_target_id == goi) %>%
#     filter(foldseek_DB %in% c("Alphafold UniProt50", "ESMAtlas30"))
#   for (qry in unique(gti_df$query)) {
#     gti_query_df <- filter(gti_df, query == qry)
#     for (db in unique(gti_query_df$foldseek_DB)) {
#       gti_query_db_df <- filter(gti_query_df, foldseek_DB == db)
#       for (domain in unique(gti_query_db_df$L2)) {
#         if (is.na(domain)) {
#           domain <- "Unknown"
#           gti_query_db_domain_df <- gti_query_db_df %>%
#             slice_max(order_by = lddt, n = 4)
#         } else {
#           gti_query_db_domain_df <- filter(gti_query_db_df, L2 == domain) %>%
#             slice_max(order_by = lddt, n = 4)
#         }
#         pdb_mods[[goi]][[qry]][[db]][[domain]] <- view_foldseek_pdb(
#           query_pdb_path = unique(gti_query_db_domain_df$query_pdb_path),
#           query_chain = unique(gti_query_db_domain_df$query_chain),
#           target_pdbs_path_list = gti_query_db_domain_df$target_pdb_path,
#           ca_inputs = FALSE,
#           align = TRUE
#         )
#       }
#     }
#   }
#   saveRDS(
#     pdb_mods,
#     glue(
#       "{wkdir}/data/interim/foldseek_results/",
#       "top_pdb_models/{Sys.Date()}_{goi}_.rds"
#     )
#   )
# }





readRDS(pdb_mod_summaries[["CD48_ENSG00000117091"]])

readRDS(pdb_mod_summaries[["CD1A_ENSG00000158477"]])

names(pdb_mod_summaries)



readRDS(
  glue(
    "{wkdir}/data/interim/foldseek_results/top_pdb_models/",
    # "2023-08-07_CD8A_ENSG00000153563_.rds" # great
    # "2023-08-07_CDH1_ENSG00000039068_.rds" # great examples
    # "2023-08-07_CCL4_ENSG00000275302_.rds" # good
    # "2023-08-07_CCL3_ENSG00000278567_.rds" #
    # "2023-08-07_CD86_ENSG00000114013_.rds" #  great
    # "2023-08-07_IL17A_ENSG00000112115_.rds" # great
    # "2023-08-07_CD79B_ENSG00000007312_.rds" # pretty cool, some are messy
    # "2023-08-07_IL6ST_ENSG00000134352_.rds" # very messy alns
    # "2023-08-07_TGFB3_ENSG00000119699_.rds"
    # "2023-08-07_CCL4L2_ENSG00000276070_.rds"
    # "2023-08-07_CCL4_ENSG00000275302_.rds"
    # "2023-08-07_CRHR1_ENSG00000120088_.rds"
    # "2023-08-07_ENG_ENSG00000106991_.rds"
    # "2023-08-07_IL10_ENSG00000136634_.rds"
    # "2023-08-07_C3_ENSG00000125730_.rds"
    # "2023-08-07_C5_ENSG00000106804_.rds"
  )
)

readRDS(pdb_mod_summaries[["C5_ENSG00000106804"]])




# tst <- readRDS("/central/groups/MazmanianLab/joeB/Microbiota-Immunomodulation/Microbiota-Immunomodulation/data/interim/foldseek_results/top_pdb_models/2023-08-07_IL9_ENSG00000145839_.rds")

# pdb_mods
# pdb_mods$CDH1_ENSG00000039068$`1O6S.pdb_A`$`Alphafold UniProt50`$Bacteria
# pdb_mods$CCL3L3_ENSG00000277768$`3H44.pdb_A`$`Alphafold UniProt50`$Eukaryota



# rcsb_aln_df <- readRDS(
#   glue(
#     "{wkdir}/data/interim/pdb_search/",
#     "2023-08-07_RCSB-chain-target-global-alignment.rds"
#   )
# )

# rcsb_aln_df %>%
#   filter(grepl("1O6S", pdb_path)) %>%
#   pull(pdb_path)
#   View






















# hmp_2012_metadata <- readRDS(
#   glue("{wkdir}/data/interim/curatedMetagenomicData/hmp_2012_metadata.rds")
# )
# hmp_2012_relab <- readRDS(
#   glue("{wkdir}/data/interim/curatedMetagenomicData/hmp_2012_relab.rds")
# )
# hmp_2012_relab_ncbi <- readRDS(
#   glue("{wkdir}/data/interim/curatedMetagenomicData/hmp_2012_relab_ncbi.rds")
# )
# # split data by phylogeny
# altExps(hmp_2012_relab) <- splitByRanks(hmp_2012_relab)
# altExps(hmp_2012_relab_ncbi) <- splitByRanks(hmp_2012_relab_ncbi)
# # create phyloseq
# phy_hmp <- makePhyloseqFromTreeSummarizedExperiment(
#   hmp_2012_relab,
#   abund_values = "relative_abundance"
# )
# phy_hmp_ncbi <- makePhyloseqFromTreeSummarizedExperiment(
#   hmp_2012_relab_ncbi,
#   abund_values = "relative_abundance"
# )
# taxids_map <- data.frame(
#   "sciname" = phy_hmp %>% taxa(),
#   "taxid" = phy_hmp_ncbi %>% taxa()
# ) %>%
#   mutate(sciname_formatted = sciname %>%
#     str_replace_all("sp. ", "sp_") %>%
#     str_replace_all(" ", "_") %>%
#     str_replace_all(":", "_") %>%
#     str_replace_all("\\[", "") %>%
#     str_replace_all("\\]", ""))


# ps_species_merged <- phy_hmp_ncbi %>%
#   merge_samples("body_site", fun = mean) %>%
#   microbiome::transform("clr")
# phyloseq::sample_data(ps_species_merged)$body_site <-
#   rownames(phyloseq::sample_data(ps_species_merged))
# hmp_ncbi_df <- psmelt(ps_species_merged) %>%
#   dplyr::rename("taxid" = "OTU") %>%
#   mutate(taxid = as.integer(taxid)) %>%
#   dplyr::select(taxid, Sample, Abundance, body_site)

# hmp_ncbi_df %>% glimpse

# foldseek_res_df_hmp_annot <- foldseek_res_df %>%
#   dplyr::left_join(hmp_ncbi_df, relationship = "many-to-many")

# foldseek_res_df_hmp_annot %>% glimpse

# foldseek_res_df_hmp_annot %>%
#   ggplot(aes(x = bits, y = gene_target_id, color = pdb_source)) +
#   geom_point() +
#   facet_grid(body foldseek_DB)
#   theme_bw()


# top_gene_target_id <- foldseek_res_df_hmp_annot %>%
#   group_by(gene_target_id) %>%
#   summarize(mean_bits = mean(bits)) %>%
#   slice_max(order_by = mean_bits, n = 200, with_ties = FALSE) %>%
#   slice_min(order_by = mean_bits, n = 50, with_ties = FALSE)

# p <- foldseek_res_df_hmp_annot %>%
#   filter(gene_target_id %in% top_gene_target_id$gene_target_id) %>%
#   ggplot(aes(x = Abundance, y = lddt, color = taxname)) +
#   facet_wrap( ~ body_site) +
#   geom_point() +
#   theme(legend.position = "none")

# ggsave("testplot.png", p, width = 10, height = 10)


# Are distributions of scores different between databases?





















# p <- blastp_res %>%
#   ggplot(aes(x = pident, y = -log10(evalue))) +
#   geom_point(aes(color = file), alpha = 0.6) +
#   scale_color_viridis_d() +
#   theme_bw() +
#   theme(legend.position = "none")

# ggsave(
#   glue("{wkdir}/figures/s.aureus-blastp/blastp_res.png"),
#   p,
#   width = 8, height = 8
# )

# p_filt <- blastp_res %>%
#   filter(pident > 30, evalue < 0.05) %>%
#   ggplot(aes(x = pident, y = -log10(evalue))) +
#   geom_point(aes(color = file), alpha = 0.6) +
#   scale_color_viridis_d() +
#   theme_bw() +
#   theme(legend.position = "none")

# ggsave(
#   glue("{wkdir}/figures/s.aureus-blastp/blastp_res-trimed.png"),
#   p_filt,
#   width = 8, height = 8
# )





# ______________________________________________________________________________
# # MGNIFY - HMP 2012 - Staphylococcus aureus ----
# # ______________________________________________________________________________

library(plyr)
library(tidyverse)
library(reshape2)
library(phyloseq)
library(microbiome)
library(httr)
library(urltools)
library(MGnifyR)

# Set up the MGnify client instance
mgclnt <- mgnify_client(
  usecache = T,
  cache_dir = glue("{wkdir}/tmp/MGnify_cache")
)

# Searching for HMP dataset
human_samples <- mgnify_query(
  mgclnt, "samples",
  biome_name = "Human", maxhits = 1000
)

human_samples %>% dim
human_samples %>% glimpse
# human_samples %>% View



# # Set up the MGnify client instance
# mgclnt <- mgnify_client(
#   usecache = T,
#   cache_dir = glue("{wkdir}/tmp/MGnify_cache")
# )
# # Retrieve the list of analyses associated with a study
# hmp_acc_list <- mgnify_analyses_from_studies(
#   mgclnt, "MGYS00001056",
#   usecache = TRUE
# )
# hmp_acc_oral <- mgnify_analyses_from_studies(
#   mgclnt, "MGYS00005569",
#   usecache = TRUE
# )
# hmp_acc_skin <- mgnify_analyses_from_studies(
#   mgclnt, "MGYS00000604",
#   usecache = TRUE
# )

# intersect(hmp_acc_list, hmp_acc_oral)
# intersect(hmp_acc_list, hmp_acc_skin)
# intersect(hmp_acc_oral, hmp_acc_skin)
# # setdiff(hmp1_accession_list, hmp2_accession_list)

# # Download all associated study/sample and analysis metadata
# meta_dataframe <- mgnify_get_analyses_metadata(mgclnt,
#   accession_list,
#   usecache = T
# )
# meta_dataframe %>% glimpse

# # Convert analyses outputs to a single `phyloseq` object
# psobj <- mgnify_get_analyses_phyloseq(
#   mgclnt, meta_dataframe$analysis_accession,
#   usecache = T
# )
# psobj
# psobj %>% taxa %>% keep(grepl("1284", .))


# # Retrieve Interpro assignment counts for these analyses
# ip_df <- mgnify_get_analyses_results(
#   mgclnt, meta_dataframe$analysis_accession,
#   retrievelist = c("interpro-identifiers"), usecache = T
# )
# head(ip_df)
# ip_df %>% str



# names(MGnifyR::mgnify_analyses_from_studies)

#______________________________________________________________________________

library("ggplot2")
theme_set(theme_bw())
library("rnaturalearth")
library("rnaturalearthdata")
library(sf)

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)








# # concatenate individual fasta files into a multi-fasta file for each gene-target-id
# # goi_df %>% glimpse
# goi_df %<>%
#   mutate(
#     target_faa_path =
#       gsub("target_pdbs", "target_fasta", target_pdb_path) %>% gsub(".pdb$", ".faa", .),
#     target_faa_path_exists = file.exists(target_faa_path)
#   )

# # for goi in gois ...
# goi <- "TGFB3_ENSG00000119699"
# concat_dir <-
#   glue("{wkdir}/data/interim/foldseek_results/target_fasta_concatenated")
# clus_dir <-
#   glue("{wkdir}/data/interim/foldseek_results/target_fasta_clustered")
# dir.create(concat_dir, showWarnings = FALSE)
# dir.create(clus_dir, showWarnings = FALSE)
# concat_fasta <- glue("{concat_dir}/{goi}_all.faa")

# faa_paths <- goi_df %>%
#   filter(target_faa_path_exists) %>%
#   filter(gene_target_id == goi) %>%
#   pull(target_faa_path) %>%
#   unique()

# pb <- progress_bar$new(total = length(faa_paths),
#   format = "[:bar] :current/:total (:percent)"
#   )

# for (f in faa_paths) {
#   pb$tick()
#   shell_do(glue("cat {f} >> {concat_fasta}"))
# }
# mmseq_cmd <- glue(
#   "mmseqs easy-linclust",
#   " {concat_fasta}",
#   " {clus_dir}/{goi}",
#   " {wkdir}/tmp",
#   " --min-seq-id 0.9",
#   " -c 0.95",
#   " --cov-mode 1"
# )
# shell_do(mmseq_cmd)


uhgg_metadata %>% glimpse
uhgg_metadata$Genome_type %>% table
uhgg_metadata %>%
  filter(Genome == Species_rep) %>%
  pull(Genome_type) %>% table
  

# of the 4,744 representative species genomes, 909 are cultured isolates
# Isolate     MAG
#     909    3835



# diamond_results <- list.files(
#   glue("{fldsk_dir}/uhgg_diamond_results"),
#   pattern = "*.m8", full.names = TRUE
# ) %>%
#   keep(fs::path_ext_remove(basename(.)) %in% tois)
# diamond_results %>% length

# future::plan(multisession, workers = 32)
# mgnfiy_diamond_res_raw <- diamond_results %>%
#   keep(file.size(.) > 0) %>%
#   keep(fs::path_ext_remove(basename(.)) %in% tois) %>%
#   furrr::future_map_dfr(
#     ~ read.delim(
#       .,
#       header = FALSE,
#       sep = "\t"
#     )
#   )
# future::plan(sequential)

# colnames(mgnfiy_diamond_res_raw) <- diamond_header
# mgnfiy_diamond_res_raw %>% glimpse

# select best target-query score by bit score







# What gene_targe_ids are found in each database?

atree <- ape::read.tree(glue("{wkdir}/data/input/mgnify/ar122_iqtree.nwk"))
btree <- ape::read.tree(glue("{wkdir}/data/input/mgnify/bac120_iqtree.nwk"))

shared_bacteria <- uhgg_results_df %>%
  filter(grepl("d__Bacteria", Lineage)) %>%
  pull(dmnd_genome) %>%
  unique() %>%
  intersect(btree$tip.label)
shared_archaea <- uhgg_results_df %>%
  filter(grepl("d__Archaea", Lineage)) %>%
  pull(dmnd_genome) %>%
  unique() %>%
  intersect(atree$tip.label)

atree <- filter_tree(atree, shared_archaea)
btree <- filter_tree(btree, shared_bacteria)


genome_eid_stats_uhgg <- uhgg_results_df %>%
  group_by(dmnd_genome, gene_target_id) %>%
  summarize(genome_hits = n(), .groups = "drop")

genome_eid_stats_uhgg_mat <- genome_eid_stats_uhgg %>%
  pivot_wider(
    names_from = "gene_target_id",
    values_from = "genome_hits",
    values_fill = 0
  ) %>%
  column_to_rownames(var = "dmnd_genome")

# uhgg_order <- seriate_matrix_rows(genome_eid_stats_uhgg_mat)
uhgg_eid_order <- seriate_matrix_rows(t(genome_eid_stats_uhgg_mat))
uhgg_genome_order <- c(btree$tip.label, atree$tip.label)


uhgg_lineage_df <- uhgg_results_df %>%
  select(dmnd_genome, Lineage) %>%
  distinct() %>%
  tidyr::separate(Lineage,
    into =
      c(
        "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"
      ),
    sep = ";",
    remove = FALSE
  )
uhgg_bact <- uhgg_lineage_df %>%
  filter(Domain == "d__Bacteria") %>%
  pull(dmnd_genome) %>%
  unique()
uhgg_arch <- uhgg_lineage_df %>%
  filter(Domain == "d__Archaea") %>%
  pull(dmnd_genome) %>%
  unique()
phy_cols <- c(
  RColorBrewer::brewer.pal(8,"Set2"),
  RColorBrewer::brewer.pal(9, "Set1"),
  RColorBrewer::brewer.pal(7,"Set3")
)

# names(phy_cols) <- unique(genome_query_df$Phylum)
p_atree <- ggtree(atree, size = 0.2)
p_btree <- ggtree(btree, size = 0.2)

heatmap_uhgg_eid_bact <- genome_eid_stats_uhgg %>%
  filter(dmnd_genome %in% uhgg_bact) %>% 
  mutate(cytchm = case_when(
    gene_target_id %in% sec_cytokines_chemokines ~ "Cytokine/Chemokine",
    TRUE ~ "Other"
  )) %>%
  mutate(dmnd_genome = factor(dmnd_genome, levels = uhgg_genome_order)) %>%
  mutate(gene_target_id = factor(gene_target_id, levels = uhgg_eid_order)) %>%
  ggplot(aes(y = dmnd_genome, x = gene_target_id, fill = genome_hits)) +
  geom_tile() +
  labs(
    x = NULL, y = NULL,
    fill = "Unique Foldseek Targets\n with Homologs in Genome"
  ) +
  scale_fill_viridis_c(option = "F", limits = c(1, 200), oob = scales::squish) +
  theme_bw() +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    title.theme = element_text(size = 12),
    direction = "horizontal",
    ticks = FALSE,
    barwidth = 12,
    barheight = 1,
    label.theme = element_text(size = 10)
  )) +
  facet_grid(cols = vars(cytchm), space = "free", scales = "free") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
  )

heatmap_uhgg_eid_arch <- genome_eid_stats_uhgg %>%
  filter(dmnd_genome %in% uhgg_arch) %>% 
  mutate(cytchm = case_when(
    gene_target_id %in% sec_cytokines_chemokines ~ "Cytokine/Chemokine",
    TRUE ~ "Other"
  )) %>%
  mutate(dmnd_genome = factor(dmnd_genome, levels = uhgg_genome_order)) %>%
  mutate(gene_target_id = factor(gene_target_id, levels = uhgg_eid_order)) %>%
  ggplot(aes(y = dmnd_genome, x = gene_target_id, fill = genome_hits)) +
  geom_tile() +
  labs(
    x = NULL, y = NULL,
    fill = "Unique Foldseek Targets\n with Homologs in Genome"
  ) +
  scale_fill_viridis_c(option = "F", limits = c(1, 200), oob = scales::squish) +
  theme_bw() +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    title.theme = element_text(size = 12),
    direction = "horizontal",
    ticks = FALSE,
    barwidth = 12,
    barheight = 1,
    label.theme = element_text(size = 10)
  )) +
  facet_grid(cols = vars(cytchm), space = "free", scales = "free") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
  )

phylum_bar <- uhgg_lineage_df %>%
  mutate(dmnd_genome = factor(dmnd_genome, levels = uhgg_genome_order)) %>%
  mutate(root = ".") %>%
  ggplot(aes(x = root, y = dmnd_genome, fill = Phylum)) +
  geom_tile() +
  scale_fill_manual(values = phy_cols) +
  theme_void()

heatmap_assembled_bact <- heatmap_uhgg_eid_bact %>%
  insert_left(phylum_bar, width = .02) %>%
  insert_left(p_btree, width = .2)
heatmap_assembled_arch <- heatmap_uhgg_eid_arch %>%
  insert_left(phylum_bar, width = .02) %>%
  insert_left(p_atree, width = .2)

ggsave(
  glue(
    "{wkdir}/figures/foldseek/",
    "{Sys.Date()}_DIAMOND_UHGG-Bacteria_eid-genome-heatmap.png"
  ),
  heatmap_assembled_bact,
  width = 22, height = 18,
  dpi = 600
)
ggsave(
  glue(
    "{wkdir}/figures/foldseek/",
    "{Sys.Date()}_DIAMOND_UHGG-Archaea_eid-genome-heatmap.png"
  ),
  heatmap_assembled_arch,
  width = 14, height = 6,
  dpi = 600
)













##### 1. UHGG Representative Species
##### 2. UHGP-95 Protein Catalog
##### 3. hCom2 Synthetic Human Gut Microbiome

# hcom_results_df <- readRDS(
#   glue("{fldsk_dir}/2023-08-17_DIAMOND-results_hCom2_plus_metadata.rds")
# )
# uhgp95_results_df <- readRDS(
#   glue("{fldsk_dir}/2023-08-17_DIAMOND-results_UHGP95_plus_metadata.rds")
# )
# uhgg_results_df <- readRDS(
#   glue("{fldsk_dir}/2023-08-17_DIAMOND-results_UHGG_plus_metadata.rds")
# )

uhgg_results_df %>% glimpse



# mgnfiy_diamond_res <- mgnfiy_diamond_res_raw %>%
#   as.data.table() %>%
#   mutate(Genome = str_before_first(dmnd_target, "_")) %>%
#    left_join(uhgg_metadata) %>%
#   inner_join(goi_df_trim, by = "target", relationship = "many-to-many")

# mgnfiy_diamond_res_max <- mgnfiy_diamond_res %>%
#   group_by(query_pdb, target, Genome) %>%
#   slice_max(order_by = dmnd_bitscore, n = 1, with_ties = FALSE)
# mgnfiy_diamond_res_max %>% glimpse

# mgnfiy_diamond_res_max_eid <- mgnfiy_diamond_res_max %>%
#   group_by(gene_target_id, Genome) %>%
#   slice_max(order_by = dmnd_bitscore, n = 1, with_ties = FALSE)
# mgnfiy_diamond_res_max_eid %>% glimpse

# saveRDS(
#   mgnfiy_diamond_res,
#   glue(
#     "{wkdir}/data/interim/foldseek_results/",
#     "{Sys.Date()}_UHGG-rep-species_diamond_results.rds"
#   )
# )
# saveRDS(
#   mgnfiy_diamond_res_max,
#   glue(
#     "{wkdir}/data/interim/foldseek_results/",
#     "{Sys.Date()}_UHGG-rep-species_diamond_results_genome-max.rds"
#   )
# )
# saveRDS(
#   mgnfiy_diamond_res_max_eid,
#   glue(
#     "{wkdir}/data/interim/foldseek_results/",
#     "{Sys.Date()}_UHGG-rep-species_diamond_results_genome-eid-max.rds"
#   )
# )


# mgnfiy_diamond_res <- readRDS(
#   glue(
#     "{wkdir}/data/interim/foldseek_results/",
#     "2023-08-11_UHGG-rep-species_diamond_results.rds"
#   )
# )
# mgnfiy_diamond_res_max <- readRDS(
#   glue(
#     "{wkdir}/data/interim/foldseek_results/",
#     "2023-08-13_UHGG-rep-species_diamond_results_genome-max.rds"
#   )
# )
mgnfiy_diamond_res_max_eid <- readRDS(
  glue(
    "{wkdir}/data/interim/foldseek_results/",
    "2023-08-14_UHGG-rep-species_diamond_results_genome-eid-max.rds"
  )
)

bit_score_matrix <- mgnfiy_diamond_res_max_eid %>%
  ungroup() %>%
  dplyr::select(target, Genome, dmnd_bitscore) %>%
  group_by(target, Genome) %>%
  summarise(dmnd_bitscore = mean(dmnd_bitscore), .groups = "drop") %>%
  pivot_wider(
    names_from = "target",
    values_from = "dmnd_bitscore",
    values_fill = 0
  ) %>%
  column_to_rownames(var = "Genome") %>%
  as.matrix()
bit_score_matrix %>% dim



# query_order <-
#   seriate_matrix_rows(t(bit_score_matrix), dist_method = "binary")
# saveRDS(
#   query_order,
#   glue(
#     "{wkdir}/data/interim/foldseek_results/",
#     "{Sys.Date()}_UHGG-rep-species_query_order.rds"
#   )
# )
query_order <- readRDS(
  glue(
    "{wkdir}/data/interim/foldseek_results/",
    "2023-08-13_UHGG-rep-species_query_order.rds"
  )
)

atree <- ape::read.tree(glue("{wkdir}/data/input/mgnify/ar122_iqtree.nwk"))
btree <- ape::read.tree(glue("{wkdir}/data/input/mgnify/bac120_iqtree.nwk"))

shared_bacteria <- genome_query_df %>%
  filter(Domain == "d__Bacteria") %>%
  pull(Genome) %>%
  unique() %>%
  intersect(btree$tip.label)
shared_archaea <- genome_query_df %>%
  filter(Domain == "d__Archaea") %>%
  pull(Genome) %>%
  unique() %>%
  intersect(atree$tip.label)

atree <- filter_tree(atree, shared_archaea)
btree <- filter_tree(btree, shared_bacteria)
# genome_order <- c(btree$tip.label, atree$tip.label)

genome_query_df <- mgnfiy_diamond_res_max_eid %>%
  mutate(detected = TRUE) %>%
  mutate(target = factor(target, levels = query_order)) %>%
  tidyr::separate(Lineage,
    into =
      c(
        "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"
      ),
    sep = ";",
    remove = FALSE
  ) %>%
  mutate(
    # Phylum = case_when(
    #   grepl("p__Firmicutes", Phylum) ~ "p__Firmicutes",
    #   TRUE ~ Phylum
    # ),
    Phylum = gsub("p__", "", Phylum)
  ) %>%
    drop_na(Genome)

genome_query_df %>% glimpse


# genome_query_df$detected %>% unique
# genome_query_df %>% glimpse
genome_query_df$Domain %>% table()
genome_query_df_bact <- genome_query_df %>%
  filter(Domain == "d__Bacteria") %>%
  mutate(Genome = factor(Genome, levels = btree$tip.label), ordered = TRUE)
genome_query_df_arch <- genome_query_df %>%
  filter(Domain == "d__Archaea") %>%
  mutate(Genome = factor(Genome, levels = atree$tip.label), ordered = TRUE)
phy_cols <- c(
  RColorBrewer::brewer.pal(8,"Set2"),
  RColorBrewer::brewer.pal(9, "Set1"),
  RColorBrewer::brewer.pal(7,"Set3")
)
names(phy_cols) <- unique(genome_query_df$Phylum)
p_atree <- ggtree(atree, size = 0.2)
p_btree <- ggtree(btree, size = 0.2)

# Plotting functions
ggpk_heatmap <- function(...) {
  ggpacket(...) %+%
    geom_tile(.id = "tile", color = NA, ...) %+%
    theme_bw() %+%
    labs(y = NULL, x = "FoldSeek PDB Hits") %+%
    scale_fill_manual(values = c("TRUE" = "black")) %+%
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank()
    )
}
ggpk_clean_bars <- function(...) {
  ggpacket(...) %+%
    geom_col(position = "identity", fill = "black") %+%
    labs(x = NULL, y = NULL) %+%
    theme_bw() %+%
    theme(
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank()
    )
}
ggpk_side_bar <- function(...) {
  ggpacket(...) %+%
    geom_tile() %+%
    scale_fill_manual(values = phy_cols) %+%
    theme_void()
}

assemble_heatmap <- function(df, tree_plot) {
  heatmap <- df %>%
    ggplot(aes(x = target, y = Genome, fill = detected)) +
    ggpk_heatmap()
  query_counts <- df %>%
    group_by(target) %>%
    summarize(target_counts = n()) %>%
    ggplot(aes(x = target, y = target_counts)) +
    ggpk_clean_bars() +
    scale_y_log10() +
    theme(axis.text.x = element_blank())
  gen_counts <- df %>%
    group_by(Genome) %>%
    summarize(gen_counts = n()) %>%
    ggplot(aes(y = Genome, x = gen_counts)) +
    ggpk_clean_bars() +
    scale_x_log10() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  phylum_bar <- df %>%
    mutate(root = ".") %>%
    ggplot(aes(x = root, y = Genome, fill = Phylum)) +
    ggpk_side_bar()
  heatmap_assembled <- heatmap %>%
    insert_top(query_counts, height = .05) %>%
    insert_right(gen_counts, width = .1) %>%
    insert_left(phylum_bar, width = .05) %>%
    insert_left(tree_plot, width = .6)
  return(heatmap_assembled)
}


# setdiff(btree$tip.label, unique(genome_query_df_bact$Genome))
pbac_heatmap_assembled <-
  assemble_heatmap(genome_query_df_bact, tree_plot = p_btree)
parch_heatmap_assembled <-
  assemble_heatmap(genome_query_df_arch, tree_plot = p_atree)

ggsave(
  glue(
    "{wkdir}/figures/foldseek/",
    "{Sys.Date()}_DIAMOND_target-UHGG-heatmap_bacteria.png"
  ),
  pbac_heatmap_assembled,
  width = 12, height = 12
)
ggsave(
  glue(
    "{wkdir}/figures/foldseek/",
    "{Sys.Date()}_DIAMOND_target-UHGG-heatmap_archaea.png"
  ),
  parch_heatmap_assembled,
  width = 12, height = 12
)




uhgp95_results_df <- readRDS(
  glue("{fldsk_dir}/2023-08-17_DIAMOND-results_UHGP95_plus_metadata.rds")
)

uhgp95_results_df %>% glimpse

uhgp95_results_df$dmnd_target %>% unique %>% length





goi_df_trim_clean %>% glimpse
gois <- unique(goi_df_trim_clean$gene_target_id)
 

uhgp95_results_df %>%
  group_by(gene_target_id) %>%
  summarize(gtid_n = n()) %>%
  View()

names(pdb_mod_summaries)

# "CD4_ENSG00000010610"
"IL12B_ENSG00000113302"
"IL6ST_ENSG00000134352"
"CD79A_ENSG00000105369"
"CD22_ENSG00000012124"
"CD28_ENSG00000178562"
"CD86_ENSG00000114013"
"CD48_ENSG00000117091"
"CD8A_ENSG00000153563"
"CD8B_ENSG00000172116"
"CDH1_ENSG00000039068"
"CD40_ENSG00000101017"
"CD320_ENSG00000167775"
"CD3E_ENSG00000198851"


tst <- readRDS(pdb_mod_summaries[["IL12B_ENSG00000113302"]])

IL12B_ENSG00000113302$`5MJ3.pdb_A`$`Alphafold UniProt50`$Bacteria


# What is the distribution of results for each target?
top_target_hits_per_genome <-
  genome_query_df %>%
  group_by(target, Genome) %>%
  slice_max(order_by = dmnd_bitscore, n = 1, with_ties = FALSE)

# Which Genomes have the greatest number of hits?
top_target_hits_per_eid_genome <- top_target_hits_per_genome %>%
  group_by(gene_target_id, Genome) %>%
  slice_max(order_by = dmnd_bitscore, n = 1, with_ties = FALSE) %>%
  group_by(gene_target_id) %>%
  summarize(gene_target_hit_n = n())


top_target_hits_per_eid_genome$gene_target_id %>% unique

p_top_uhgg_eid_hits <- top_target_hits_per_eid_genome %>%
   mutate(cythit = gene_target_id %in% sec_cytokines_chemokines) %>%
   ggplot(aes(
     x = gene_target_hit_n,
     y = fct_reorder(gene_target_id, gene_target_hit_n)
   )) +
   scale_x_log10() +
     geom_col(aes(fill = cythit), width = 0.7) +
     labs(
       x = "Representative UHGG Genomes with Detection",
       y = NULL,
       fill = "Cytokine/Chemokine"
     ) +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
    theme_bw() +
    theme(legend.position = "top")

ggsave(
  glue(
    "{wkdir}/figures/foldseek/",
    "{Sys.Date()}_DIAMOND_top-UHGG-hits-by-EID.png"
  ),
  p_top_uhgg_eid_hits,
  width = 8, height = 10
)



p_target_n_dist <- top_target_hits_per_genome %>%
  group_by(target) %>%
  summarize(target_hit_n = n()) %>%
  ggplot(aes(target_hit_n)) +
  geom_histogram(bins = 100, fill = "lightblue") +
  theme_bw() +
  labs(x = "Number of Genomes with detection of target")

# is it often that we see multiple hits within a genome (multiple copy numbers?)
target_hits_per_genome <- genome_query_df %>%
  group_by(target, Genome) %>%
  mutate(target_hit_per_genome = n())
# target_hits_per_genome %>% glimpse

p_target_n_per_genome <- target_hits_per_genome %>%
  select(target_hit_per_genome, target, Genome) %>%
  distinct() %>%
  ggplot(aes(target_hit_per_genome)) +
  geom_histogram(bins = 100, fill = "lightblue") +
  theme_bw() +
  # scale_x_log10() +
  scale_y_log10() +
  labs(x = "Number of target hits within a genome")

p_target_n_dist + p_target_n_per_genome


# target_hits_per_genome %>%
#   select(target_hit_per_genome, target, Genome) %>%
#   distinct() %>%
#   ggplot(aes(target_hit_per_genome)) +
#   stat_ecdf(geom = "point")


# genome_query_df %>% glimpse
# genome_query_df %>% drop_na(Genome) %>% glimpse
# genome_query_df %>% drop_na(Genome)


target_hits_per_genome %>%
  drop_na(Genome) %>% 
  select(target_hit_per_genome, target, Genome) %>%
  distinct() %>%
  arrange(desc(target_hit_per_genome)) %>%
  head(n=1000)


getwd()
source_python('add.py')
# use batchtools to deplot

# query_counts <-
#   genome_query_df_bact %>%
#   # slice_sample(n = 1000) %>%
#   group_by(target) %>%
#   summarize(target_counts = n()) %>%
#   ggplot(aes(x = target, y = target_counts)) +
#   # scale_y_log10() +
#   geom_col(position = "identity", fill = "black") +
#   labs(x = NULL, y = NULL) +
#   theme(
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank()
#   )

# gen_counts <- genome_query_df %>%
#   group_by(Genome) %>%
#   summarize(gen_counts = n()) %>%
#   ggplot(aes(y = Genome, x = gen_counts)) +
#   labs(x = NULL, y = NULL) +
#   geom_col(position = "identity", fill = "black") +
#   theme_void()

# p_heatmap <- genome_query_df %>%
#   ggplot(aes(x = target, y = Genome, fill = detected)) +
#   geom_tile() +
#   theme_bw() +
#   labs(y = NULL, x = "FoldSeek PDB Hits") +
#   scale_fill_manual(values = c("TRUE" = "black")) +
#   theme(
#     panel.grid = element_blank(),
#     axis.text.x = element_blank(),
#     axis.text.y = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.ticks.y = element_blank()
#   )



# Is there any relationship between the number of hits found in a genome and completeness/contamination?
genome_query_df %>% glimpse


genome_summary_stats <- genome_query_df %>%
  group_by(Genome) %>%
  mutate(genome_hits = n()) %>%
  ungroup() %>%
  select(Genome, genome_hits, Completeness, Contamination, Phylum, Lineage) %>%
  distinct()

p_qc1 <- genome_summary_stats %>%
  ggplot(aes(x = genome_hits, y = Completeness)) +
  geom_point(alpha = 0.8, size = 2, aes(color = Phylum, group = Lineage)) +
  theme_bw() +
  scale_x_log10() +
  scale_color_manual(values = phy_cols) +
  geom_smooth(method = "loess", se = TRUE, color = "black") +
  labs(x = "Number of hits", y = "Completeness (%)")

p_qc2 <- genome_summary_stats %>%
  ggplot(aes(x = genome_hits, y = Contamination)) +
  geom_point(alpha = 0.8, size = 2, aes(color = Phylum, group = Lineage)) +
  theme_bw() +
  scale_x_log10() +
  scale_color_manual(values = phy_cols) +
  geom_smooth(method = "loess", se = TRUE, color = "black") +
  labs(x = "Number of hits", y = "Contamination (%)")


p_qc <- p_qc1  %>% insert_top(p_qc2)
ggsave(
  glue("{wkdir}/figures/foldseek/",
  "{Sys.Date()}_DIAMOND_UHGG-genome-hits-contamination-complete.png"),
  p_qc,
  width = 9, height = 9
)

ggplotly(p_qc1)
ggplotly(p_qc2)






# What taxa have the most hits across all proteins of interest?
genome_summary_stats %>%
  slice_max(genome_hits, n = 10) 
  ggplot(aes(x = genome_hits, y = Lineage)) +
  geom_col(aes(fill = Phylum), stat = "identity")


genome_summary_stats %>% View
genome_summary_stats$Genome %>%
  unique()


# mgnfiy_diamond_res %>%
#   ggplot(aes(x = Completeness, y = Contamination)) +
#   geom_point(aes(color = N50), alpha = 0.4) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))


tst <- readRDS(
  glue(
    "{wkdir}/data/interim/foldseek_results/top_pdb_models/",
    "2023-08-07_C3_ENSG00000125730_.rds"
    # "2023-08-07_IL10_ENSG00000136634_.rds"
    # "2023-08-07_TGFB3_ENSG00000119699_.rds"
  )
)

tst$C3_ENSG00000125730$`7ZGJ.pdb_B`$`Alphafold UniProt50`$Bacteria


il10_pdb <- bio3d::read.pdb(file = get.pdb("6X93", URL = TRUE))
# pdbsplit( get.pdb("6X93", URLonly=TRUE) )


pdbsplit(get.pdb("2ICW", URLonly=TRUE) )

unique(il10_pdb$atom$chain)

# IL10RA: BE

glue(
  "trill docktest 1 dock",
  " Smina",
  " split_chain/6X93_B.pdb",
  " split_chain/6X93_A.pdb"
)

# TEST FOR ZACH
# glue(
#   "trill lightdock29 0 dock",
#   " Smina",
#   " SARS2_RBD.pdb",
#   " A0A5C1RDG7.pdb"
# )

#  smina -r lightdock_29_R.pdbqt -l lightdock_29_L.pdbqt --autobox_ligand pocket1_atm.pdb
#  smina -r lightdock_29_R.pdb -l lightdock_29_L.pdb --autobox_ligand

target_gtid_freq <- merged_aln_df %>%
  select(target, gene_target_id) %>%
  distinct() %>%
  group_by(target) %>%
  summarize(count = n()) %>%
  arrange(desc(count))

non_specific_targets <- target_gtid_freq %>%
  filter(count > 1) %>%
  pull(target)

merged_aln_df %>%
  select(target, gene_target_id) %>%
  filter(target %in% non_specific_targets)
  
"CD4_ENSG00000010610"
"TGFB1_ENSG00000105329"
"TGFB2_ENSG00000092969"

"IL6ST_ENSG00000134352"
"IL12B_ENSG00000113302"
"CD79A_ENSG00000105369"
"CD79B_ENSG00000007312"
"CD22_ENSG00000012124"
"CD48_ENSG00000117091"
"CD28_ENSG00000178562"
"CD86_ENSG00000114013"
"CD8A_ENSG00000153563"
"CD8B_ENSG00000172116"
"CD48_ENSG00000117091"
"CD40_ENSG00000101017"
"CDH1_ENSG00000039068"
"C3_ENSG00000125730"
"C5_ENSG00000106804"

fldseek_target_df <- merged_aln_df %>%
  filter(gene_target_id == "TGFB2_ENSG00000092969") %>%
  mutate(tlen = tend - tstart) %>%
  arrange(desc(tlen))
  # arrange(dmnd_evalue)
fldseek_target_df %>% View()
fldseek_target_df$target %>% unique %>%  head

fldseek_target_df %>%
  group_by(target, query) %>%
  summarize(count = n()) %>%
  arrange(desc(count)) %>%
  print(n = 40)

fldseek_target_df %>%
  filter(target == "AF-G8MZC3-F1-model_v4") %>%
  pull(dmnd_genome) %>%
  unique()

# MGYP001134830988
# AF-P01024-F1-model_v4.pdb

fldsk_htmls <- glue("{wkdir}/data/interim/foldseek_results/html_output")

fldseek_target_df %>%
  pull(query_pdb_path) %>%
  unique() %>%
  strex::str_after_nth(., "/", -2) %>%
  fs::path_ext_remove() %>%
  purrr::map_chr(
    ~ glue("{fldsk_htmls}/{.}")
  ) %>%
  list.files(full.names = TRUE) %>%
  keep(grepl("_Alphafold_UniProt50|ESMAtlas30", .))


view_foldseek_pdb(
  query_pdb_path = unique(fldseek_target_df$query_pdb_path)[1],
    query_chain = NA,
    ca_inputs = FALSE,
    align = TRUE,
    target_pdbs_path_list = unique(fldseek_target_df$target_pdb_path)[1]
)




#______________________________________________________________________________
# FOlding proteins of interest from hCOM2

#!/bin/bash
#SBATCH --time=8:00:00   # walltime
#SBATCH --ntasks-per-node=4
#SBATCH --nodes=4 # number of nodes
#SBATCH --gres=gpu:4 # number of GPUs
#SBATCH --mem-per-cpu=60G   # memory per CPU core
#SBATCH -J "tutorial"   # job name
#SBATCH --mail-user="" # change to your email
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=%x-%j.out
# master_addr=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | head -n 1)
# export MASTER_ADDR=$master_addr
# export MASTER_PORT=13579
# srun trill example_3 4 finetune esm2_t33_650M_UR50D trill/data/query.fasta --nodes 4 --strategy deepspeed_stage_2_offload
# trill fold tgfb_test.fa

# RUnning an interactive session with a gpu
srun --job-name "InteractiveGPUJob" --gres=gpu:1 --cpus-per-task 1 --mem 50G --time 30:00 --pty bash
# export TORCH_HOME=/central/groups/MazmanianLab/joeB/cache
export TRANSFORMERS_CACHE=/central/groups/MazmanianLab/joeB/cache
trill tgfb_esmfold_example 1 fold tgfb_test.fa

# trill docktest dock --force_ligand protein --save_visualisation True Smina
# trill docktestCPU 0 dock Smina split_chain/6X93_B.pdb split_chain/6X93_A.pdb

# glue(
#   "srun --gres=gpu:1",
#   " --time=40:00",
#   " --mem=60G",
#   " -J 'tutorial'",
#   " --mail-user='jboktor@caltech.edu'",
#   " --mail-type=BEGIN",
#   " --mail-type=END",
#   " --mail-type=FAIL", 
#   " --output=%x-%j.out",
#   " master_addr=$(scontrol show hostnames '$SLURM_JOB_NODELIST' | head -n 1) &&",
#   " export MASTER_ADDR=$master_addr &&",
#   " export MASTER_PORT=13579 &&",
#   " trill tgfb_esmfold_example 1 fold tgfb_test.fa"
#   # " --strategy deepspeed_stage_2_offload"
# )

pdb_l <- bio3d::read.pdb(glue("{wkdir}/6X93_A.pdbqt"))
pdb_r <- bio3d::read.pdb(glue("{wkdir}/6X93_B.pdbqt"))

r3dmol(
      viewer_spec = m_viewer_spec(
        cartoonQuality = 10,
        lowerZoomLimit = 70,
        upperZoomLimit = 750
      )
    ) %>%
    m_add_model(data = m_bio3d(bio3d::read.pdb(glue("{wkdir}/6X93_A.pdbqt")))) %>%
    m_set_style(
      sel = m_sel(model = 0),
      style = m_style_cartoon(color = "grey")
    ) %>% 
        m_add_model(data = m_bio3d(bio3d::read.pdb(glue("{wkdir}/6X93_B.pdbqt")))) %>%
    m_set_style(
      sel = m_sel(model = 1),
      style = m_style_cartoon(color = "red")
    ) %>% 
    m_zoom_to()

# reticulate::use_condaenv(condaenv = "esmfold", required = TRUE)
# source_python(glue("{src_dir}/python-scripts/helpers_esm_fold.py"))

# Load in results from hCOm2 hits
key_gois <- c(
  sec_cytokines_chemokines,
  "C3_ENSG00000125730",
  "C5_ENSG00000106804"
)

# how many proteins of interest are there? select (cytokines/chemokines/ and complements)
hcom_results_df <- readRDS(
  glue("{fldsk_dir}/2023-08-17_DIAMOND-results_hCom2_plus_metadata.rds")
)

hcom_results_df_keytargets <- hcom_results_df %>%
  filter(gene_target_id %in% key_gois &
    gene_target_id != "CD4_ENSG00000010610") %>%
  group_by(gene_target_id, target) %>%
    slice_max(order_by = dmnd_bitscore, n = 50, with_ties = FALSE)

hcom_results_df_keytargets %>%
  group_by(gene_target_id) %>%
  summarize(fastas = n())

# hcom_results_df_keytargets %>% glimpse
hcom_gid_gene_hits <- hcom_results_df_keytargets %>%
  pull(dmnd_target) %>%
  unique()

# Results inculde a total of: 788 proteins (some overlap across targets)

# Select proteins of interest from full catalog and save as fasta file
hcom2_catalog <- glue("{homedir}/Downloads/RefDBs/hCom2/hCom2.fasta")
hcom2_select_fa <- glue(
  "{wkdir}/data/interim/foldseek_results/",
  "hCom2-top-50-dmnd-hits-per-fldsk-target.fasta"
)

extract_seq_from_catalog(
  fasta_header_list = list(hcom_gid_gene_hits),
  catalog_path = hcom2_catalog,
  output_fasta_path = hcom2_select_fa
)

# read in hcom2 proteins of interest
hcom_seqs <- seqinr::read.fasta(hcom2_select_fa, seqtype = "AA")
hcom2_fa_seqs <- names(hcom_seqs) %>%
  purrr::map(
    ~ seqinr::getSequence(hcom_seqs[[.]], as.string = TRUE) %>%
      unlist() %>%
      gsub("\\*", "", .)
  )
hcom2_pdb_paths <- names(hcom_seqs) %>%
  purrr::map(
    ~ glue(
      "{wkdir}/data/interim/foldseek_results/hCom2_pdbs/{.}.pdb"
    )
  )


batch_esm_fold <- function(fasta_list, pdb_paths, src = src_dir) {
  require(glue)
  require(reticulate)
  use_condaenv(condaenv = "esmfold", required = TRUE)
  source_python(glue("{src}/python-scripts/helpers_esm_fold.py"))
  esm_fold_wrapper_cpu(
    sequence_list = list(fasta_list),
    output_paths = list(pdb_paths)
  )
}

future::plan(
  future.batchtools::batchtools_slurm,
  template = glue("{wkdir}/batchtools_templates/batchtools.slurm.tmpl"),
  resources = list(
    name = glue("{get_time()}_ESMFOLD"),
    memory = "25G",
    ncpus = 8,
    walltime = 140000
  )
)

# Chunk files (5 per job) and download
tic()
n_jobs <- ceiling(length(hcom2_pdb_paths) / 10)
# n_jobs <- 5
hcom2_esmfold_runs <- listenv()
pb <- progress_bar$new(
  format = "[:bar] :percent eta: :eta",
  total = n_jobs
)
for (job in 1:n_jobs) {
  pb$tick()
  input_fas <- chunk_func(hcom2_fa_seqs, n_jobs)[[job]] %>% unlist()
  output_paths <- chunk_func(hcom2_pdb_paths, n_jobs)[[job]] %>% unlist()
  hcom2_esmfold_runs[[job]] %<-% batch_esm_fold(
    fasta_list = input_fas,
    pdb_paths = output_paths
  )
}
toc()

# as.list(hcom2_esmfold_runs)

list.files(glue("{wkdir}/.future/20230823_142239-Mq3tqb"),
  recursive = TRUE, full.names = TRUE
) %>%
  keep(grepl("/logs/", .))



# _________________________________________________________________
# Calculate pairwise similarity to 

label_map <- c(
  "SW_fident_sequences_aa_signalp_trimmed__query_aa" = "HumanGene-to-HumanPDBs",
  "SW_fident_sequences_aa_signalp_trimmed__foldseek_target_aa" = "HumanGene-to-FoldseekPDBs",
  "SW_fident_sequences_aa_signalp_trimmed__hcom2_target_aa" = "HumanGene-to-hCom2Hits",
  "SW_fident_query_aa__foldseek_target_aa" = "HumanPDBs-to-FoldseekPDBs",
  "SW_fident_query_aa__hcom2_target_aa" = "HumanPDBs-to-hCom2Hits",
  "SW_fident_foldseek_target_aa__hcom2_target_aa" = "FoldseekPDBs-to-hCom2Hits"
)

# What is  the distribution of pIdnt scores across layers?

# Distribution plot
merged_aln_df %>% glimpse
# merged_aln_df %>% View

p_dens1 <- merged_aln_df %>%
  select(
    contains("SW_fident"),
    # contains("local_sequence_sim"),
    all_of(unname(unlist(group_dict)))
  ) %>%
  pivot_longer(
    cols = -unname(unlist(group_dict)),
    names_to = "comparison"
  ) %>%
  # pull(comparison)  %>% unique
  ggplot(aes(x = value)) +
    geom_density(aes(color = comparison)) +
    labs(x = "% Identity (Human RCSB/AF2 PDBs to hCom2 Hits)") +
    theme_bw() +
    scale_color_d3(labels = label_map) +
    theme(legend.position = "top")
  

p_dens_facet <- merged_aln_df %>%
  select(
    contains("SW_fident"),
    all_of(unname(unlist(group_dict)))
  ) %>%
  pivot_longer(
    cols = -unname(unlist(group_dict)),
    names_to = "comparison"
  ) %>%
  ggplot(aes(x = value, after_stat(density))) +
    facet_wrap(~gene_target_id, ncol = 2, scales = "free_y") +
    geom_density(aes(color = comparison)) +
    labs(x = "% Identity (Human RCSB/AF2 PDBs to hCom2 Hits)") +
    theme_bw() +
    scale_color_d3(labels = label_map) +
    theme(legend.position = "bottom")

ggsave(
  glue("{wkdir}/figures/foldseek/{Sys.Date()}_SW-fident-density.png"),
  p_dens1,
  width = 10, height = 7
)

ggsave(
  glue("{wkdir}/figures/foldseek/{Sys.Date()}_SW-fident-density-facet-by-gtid.png"),
  p_dens_facet,
  width = 14, height = 16
)





# Initiate future.batchtools backend for parallel processing
# seq_test <- list(
#   "MGATTALTATVSPEDATDKAVSYASSKISVATVNGSGVVTGVSEGSATITATTHDGSKTASTAVTVTAA",
#   "MEVTGITLNQTELNLTAGRTAALKATVLPDNAADKTVTWSSSAPEVAEVDANGTVTAKTAGSATITAQTANGRTVTCTVTVTAAEPAAEEPPKTN",
#   "MNLYATDESGNVADKEIKVIVNEPTLTIEQANVSIAVGETAQLNATVQGANQTITWTSSDESVATVDANGVVTGIKKGTATITASANGIKDTTEIRIQSNQRSSNEKSNTSSSSSKSGSSSSSSSSGSSGSTSTGTHKHTMPTGNIGKWFSSRSELVSYYNSVAEEWNDKWLSGKISNEEYYANCPSGYECWSCSYCGKWTGNFKYN"
# )
# out_test <- list(
#   glue("{wkdir}/GPUTEST_ESMFOLD_MGYG000248927_02411.pdb"),
#   glue("{wkdir}/GPUTEST_ESMFOLD_MGYG000180748_01300.pdb"),
#   glue("{wkdir}/GPUTEST_ESMFOLD_MGYG000177439_00845.pdb")
# )

# job %<-% esm_fold_wrapper(
#   sequence_list = seq_test,
#   output_paths = out_test
# )






#______________________________________________________________________________
# Downloading MGNify Human Sample Metadata

library(plyr)
library(tidyverse)
theme_set(theme_bw())
# library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(reshape2)
library(phyloseq)
library(microbiome)
# library(httr)
library(urltools)
library(MGnifyR)
library(glue)
library(viridis)
library(future)

library(remotes)
install_github("r-spatial/sf")


library(MGnifyR)

# Set up the MGnify client instance
mgclnt <- mgnify_client(
  usecache = T,
  cache_dir = glue("{wkdir}/tmp/MGnify_cache")
)
# microbiome studies
human_studies <- mgnify_query(
  mgclnt, "studies",
  biome_name = "Human", maxhits = 0
)

human_studies %>% glimpse
saveRDS(
  human_studies,
  glue("{wkdir}/data/interim/mgnify/{Sys.Date()}_mgnify_human_studies.rds")
)



human_studies <- readRDS(
  glue("{wkdir}/data/interim/mgnify/2023-09-07_mgnify_human_studies.rds")
)

# future::plan("multisession", workers = 14)
future::plan(
  future.batchtools::batchtools_slurm,
  template = glue("{wkdir}/batchtools_templates/batchtools.slurm.tmpl"),
  resources = list(
    name = "mgnify_sample_meta",
    memory = "1G",
    ncpus = 1,
    walltime = 10800
  ) #, workers = length(human_studies$accession),
)

#' Function to download all sample metadata
#' from MGnify using study accessions
get_mgnify_sample_meta <- function(acc) {
  require(magrittr)
  require(MGnifyR)
  mgclnt <- mgnify_client(
    usecache = TRUE,
    cache_dir = "/central/groups/MazmanianLab/joeB/Microbiota-Immunomodulation/Microbiota-Immunomodulation/tmp/MGnify_cache"
  )
  mgnify_analyses_from_studies(mgclnt, accession = acc, usecache = TRUE) %>%
    mgnify_get_analyses_metadata(mgclnt, ., usecache = TRUE)
}

future.apply::future_lapply(
  human_studies$accession,
  FUN = get_mgnify_sample_meta,
  future.scheduling = FALSE
)

# FUN = function(iris) sum(iris$Sepal.Length)














analyses_ps %>%
  meta() %>% glimpse
  View()


analyses_ps %>%
  meta() %>%
  pull(analysis_experiment.type) %>%
  table

analyses_ps %>%
  meta() %>%
  pull(sample_geographic.location..country.and.or.sea.region.) %>%
  table
analyses_ps %>%
  meta() %>%
  pull(study_attributes.centre.name) %>%
  table


analyses_ps %>%
  meta() %>%
  pull(sample_geographic.location..longitude.) %>%
  table


# alpha diversity analysis
metadf <- meta(analyses_ps)
observed_spec <- microbiome::alpha(analyses_ps, index = "observed")
observed_spec_df <- bind_cols(metadf, observed_spec)
observed_spec_df %>% glimpse

observed_spec_df %>%
  ggplot(aes(x=fct_reorder(study_attributes.accession, observed), y=observed)) +
  geom_point(aes(color = analysis_experiment.type))+
  theme_bw()



library(maps)
library(mapproj)

world_tbl <- map_data("world") %>% 
  as_tibble() %>%
  filter(region != "Antarctica")
world_tbl

activated_sl_meta <- analyses_ps %>%
  meta() %>%
  drop_na(c(sample_longitude, sample_latitude)) %>%
  select(
    analysis_experiment.type, sample_longitude, 
    sample_latitude, study_attributes.accession
  ) %>%
  mutate_at(vars(sample_longitude, sample_latitude), as.numeric) %>%
  group_by(sample_longitude, sample_latitude, analysis_experiment.type) %>% 
  summarize(sample_n = n()) %>% 
  ungroup()
activated_sl_meta %>% glimpse

# world base plot
world_base <- world_tbl %>%
  ggplot() +
  geom_map(
    aes(long, lat, map_id = region),
    map = world_tbl,
    color = "darkgrey", fill = "lightgrey", size = 0.3
  ) +
  theme_minimal() +
    labs(x = "longitude", y = "latitude", size = "Sample N", fill = "Method") +
    scale_size_continuous(range = c(3, 10)) +
    scale_fill_d3() +
    geom_point(
      aes(
        x = sample_longitude, y = sample_latitude,
        size = sample_n, fill = analysis_experiment.type
      ),
      data = activated_sl_meta,
      # fill = "red",
      shape = 21, alpha = 0.5
    ) +
    guides(fill = guide_legend(override.aes = list(size=5))) +
    theme(legend.position = "top")

# generate orthologonal coordinate system for map projection
final_map <- world_base +
  coord_map("harrison", dist = 1, angle = 30, ylim = c(-75, 85))
final_map

ggsave(
  glue("{wkdir}/{Sys.Date()}_mgnify_activated_sludge_sample_map.png"),
  final_map,
  width = 9, height = 6,
  dpi = 600
)


library(RColorBrewer)
library(reshape)

prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(0.01), log10(.2), length = 9), 3)
# Also define gray color palette
gray <- gray(seq(0,1,length=5))

pseq.rel <- analyses_ps %>% 
  subset_samples(analysis_experiment.type == "assembly") %>%
  core(detection = 1, prevalence = 1/95) %>% 
  transform("compositional")

p_ <- pseq.rel %>%
  plot_core(
    plot.type = "heatmap",
    colours = viridis::viridis_pal(option = "")(8),
    prevalences = prevalences,
    detections = detections,
    min.prevalence = prevalence(pseq.rel, sort = TRUE)[100]
  ) +
  labs(x = "Detection Threshold\n(Relative Abundance (%))") +
  #Adjusts axis text size and legend bar height
  theme(
    axis.text.y = element_text(size = 8, face = "italic"),
    axis.text.x.bottom = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )



which(taxa_names(analyses_ps) %nin% names(ncbi_lineage_list))

ncbi_lineage_list[taxa_names(analyses_ps)]



ncbi_lineage_list <- ncbi_lineages$name
names(ncbi_lineage_list) <- ncbi_lineages$taxid




# USeful resources for plotting maps
# http://freerangestats.info/blog/2017/06/04/military-gdp
# https://datavizm20.classes.andrewheiss.com/example/12-example/




# ncbi_lineages <- read.delim(
#   glue(
#     "/central/groups/MazmanianLab/joeB/Downloads/",
#     "2023-08-02_ncbi-taxdump/fullnamelineage.dmp"
#   ),
#   quote = "",
#   header = FALSE
# ) %>%
#   select(-c(V2, V4, V6)) %>%
#     dplyr::rename(taxid = V1, name = V3, lineage = V5) %>%
#     mutate(lineage = case_when(
#       grepl("Viruses", lineage) ~ paste0("Viruses; ", lineage),
#       TRUE ~ lineage
#     )) %>%
#   select(taxid, name) %>%
#   drop_na() %>%
#   distinct()






tax_table(analyses_ps) %>%
  as.data.frame() %>%
  View


# __ DMMS on Assembled genomes
# get a list of samples with fewer than 1000 counts total
sample_sums <- microbiome::abundances(analyses_ps) %>%
  as.data.frame() %>%
  colSums()
low_read_samples <- names(sample_sums[sample_sums < 10])

pseq <- analyses_ps %>% 
  subset_samples(analysis_experiment.type == "assembly") %>%
  subset_samples(analysis_accession %nin% low_read_samples) %>%
  core(detection = 0, prevalence = 0.1)
pseq

# get a list of samples with fewer than 1000 counts total
sample_sums_rnd2 <- microbiome::abundances(pseq) %>%
  as.data.frame() %>%
  colSums()
low_read_samples_rnd2 <- names(sample_sums_rnd2[sample_sums_rnd2 == 0])
pseq %<>% subset_samples(analysis_accession %nin% low_read_samples_rnd2)


# Pick the OTU count matrix
# and convert it into samples x taxa format
dat <- analyses_ps %>%
  subset_samples(analysis_accession %in% sample_names(pseq)) %>%
  prune_taxa(taxa(pseq), .) %>%
  abundances()
count <- as.matrix(t(dat))
min(rowSums(count))


fit <-  lapply(1:4, dmn, count = count, verbose=TRUE)
fit

lplc <- sapply(fit, laplace)
aic <- sapply(fit, AIC)
bic <- sapply(fit, BIC)

# Selecting optimal model
best <- fit[[which.min(unlist(lplc))]]

mixturewt(best)
clust_assignments <- apply(mixture(best), 1, which.max)



# Aitchison distance PCA

# Pick core taxa with with the given prevalence and detection limits
ps_pca <- pseq %>%
  microbiome::transform("clr")

idist <- phyloseq::distance(ps_pca, method = "euclidean")
dist_label <- "Aitchison"
imds <- phyloseq::ordinate(ps_pca, "MDS", distance = idist)

p <- plot_ordination(ps_pca, imds,  axes = c(1, 2, 3))
pcoa_data <- pca_df <- bind_cols(p$data, clust_assignments) %>%
  dplyr::rename(DMM_cluster = ...80) %>%
  mutate(DMM_cluster = factor(DMM_cluster)) 
pca_df %>% glimpse

pco1 <- round((imds$values$Relative_eig[1]) * 100, digits = 2)
pco2 <- round((imds$values$Relative_eig[2]) * 100, digits = 2)
pco3 <- round((imds$values$Relative_eig[3]) * 100, digits = 2)

fig_plotly <-
  pcoa_data %>%
  plot_ly(
    x = ~Axis.1,
    y = ~Axis.2,
    z = ~Axis.3 #,
    # text = ~ paste("Line:", line, "<br>Organism:", organism)
  ) %>%
  layout(
    scene = list(
      xaxis = list(title = glue("PCo1 ({pco1}%)")),
      yaxis = list(title = glue("PCo2 ({pco2}%)")),
      zaxis = list(title = glue("PCo3 ({pco3}%)"))
    )
  )

plotly_pcoa_richness <-
  add_markers(
    fig_plotly,
    color = ~DMM_cluster
  )




library(RColorBrewer)

# # get a list of samples with fewer than 1000 counts total 
# sample_sums <- microbiome::abundances(analyses_ps) %>% as.data.frame() %>% colSums()
# low_read_samples <- names(sample_sums[sample_sums < 1000])
# pseq <- analyses_ps %>% 
#   subset_samples(analysis_experiment.type == "amplicon") %>% 
#   subset_samples(analysis_accession %nin% low_read_samples) %>% 
#   microbiome::transform("compositional") %>%
#   core(detection = 1e-5, prevalence = 0.1)
# # get a list of samples with fewer than 1000 counts total 
# sample_sums_rnd2 <- microbiome::abundances(pseq) %>% as.data.frame() %>% colSums()
# low_read_samples_rnd2 <- names(sample_sums_rnd2[sample_sums_rnd2 == 0])
# ps_trim <- pseq %>% subset_samples(analysis_accession %nin% low_read_samples_rnd2)

prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(0.001), log10(.2), length = 9), 3)

pseq.rel <- ps_pca %>% 
  microbiome::transform("compositional")

p_core_taxa <-
  pseq.rel %>%
  plot_core(
    plot.type = "heatmap",
    colours = viridis::viridis_pal(option = "H")(8),
    prevalences = prevalences,
    detections = detections,
    min.prevalence =min(prevalence(pseq.rel, sort = TRUE))
  ) +
  labs(x = "Detection Threshold\n(Relative Abundance (%))", y= "Taxa") +
  #Adjusts axis text size and legend bar height
  theme(
    axis.text.y = element_blank(), # element_text(size = 8, face = "italic"),
    axis.text.x.bottom = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.key.height = unit(1, 'cm'),
    axis.ticks.y = element_blank()
  )

ggsave(
  glue("{Sys.Date()}_activated_sludge_amplicon_seq_core_taxa.png"),
  p_core_taxa,
  width = 4, height = 4
)




# for (k in seq(ncol(fitted(best)))) {
#   d <- melt(fitted(best))
#   colnames(d) <- c("OTU", "cluster", "value")
#   d <- subset(d, cluster == k) %>%
#      # Arrange OTUs by assignment strength
#      arrange(value) %>%
#      mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
#      # Only show the most important drivers
#      filter(abs(value) > quantile(abs(value), 0.8))     

#   p <- ggplot(d, aes(x = OTU, y = value)) +
#        geom_bar(stat = "identity") +
#        coord_flip() +
#        labs(title = paste("Top drivers: community type", k))
#   print(p)
# }


# analyses_ps %>%
#   taxa() %>%
#   keep(grepl("83333", .))







# Core analyses ---------------

library(RColorBrewer)
library(reshape)
`%nin%` <- Negate(`%in%`)

# get a list of samples with fewer than 1000 counts total 
sample_sums <- microbiome::abundances(analyses_ps) %>% as.data.frame() %>% colSums()
low_read_samples <- names(sample_sums[sample_sums < 1000])
pseq <- analyses_ps %>% 
  subset_samples(analysis_experiment.type == "amplicon") %>% 
  subset_samples(analysis_accession %nin% low_read_samples) %>% 
  microbiome::transform("compositional") %>%
  core(detection = 1e-5, prevalence = 0.1)
# get a list of samples with fewer than 1000 counts total 
sample_sums_rnd2 <- microbiome::abundances(pseq) %>% as.data.frame() %>% colSums()
low_read_samples_rnd2 <- names(sample_sums_rnd2[sample_sums_rnd2 == 0])
ps_trim <- pseq %>% subset_samples(analysis_accession %nin% low_read_samples_rnd2)


prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(0.001), log10(.2), length = 9), 3)

pseq.rel <- ps_trim %>% 
  microbiome::transform("compositional")

p_core_taxa <-
  pseq.rel %>%
  plot_core(
    plot.type = "heatmap",
    colours = viridis::viridis_pal(option = "H")(8),
    prevalences = prevalences,
    detections = detections,
    min.prevalence =min(prevalence(pseq.rel, sort = TRUE))
  ) +
  labs(x = "Detection Threshold\n(Relative Abundance (%))", y= "Taxa") +
  #Adjusts axis text size and legend bar height
  theme(
    axis.text.y = element_blank(), # element_text(size = 8, face = "italic"),
    axis.text.x.bottom = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.key.height = unit(1, 'cm'),
    axis.ticks.y = element_blank()
  )

ggsave(
  glue("{Sys.Date()}_activated_sludge_amplicon_seq_core_taxa.png"),
  p_core_taxa,
  width = 4, height = 4
)











