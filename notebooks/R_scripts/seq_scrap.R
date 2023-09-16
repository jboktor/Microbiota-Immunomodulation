# ______________________________________________________________________________
# NEO2 Variants Sequence Analysis ----

# read in fasta files ----
fasta_paths <- list.files(
  glue("{proj_dir}/fastas_round2"),
  full.names = TRUE
)

# load in fasta files ----

process_fasta_file <- function(path) {
  s1 <- protr::readFASTA(path)
  names(s1) <- paste0("neo2_", 1:length(s1))
  s1 %<>%
    protr::removeGaps(pattern = "X", replacement = "G") %>%
    purrr::set_names(names(s1))
  return(s1)
}

fastas <- fasta_paths %>%
  purrr::set_names(fs::path_ext_remove(basename(.))) %>%
  purrr::map( ~ process_fasta_file(.))

s1 <- protr::readFASTA(fasta_paths[[1]])
control_seq <- s1[["neo2/15"]]

# add first layer of list names to second layer of list names ----
fastas_merged <- fastas %>% unlist

# calculate sequence similarity to Neo2 OG ----
seqalign <-
  purrr::map(
    fastas_merged,
    ~ protr::twoSeqSim(control_seq, .,
      type = "global",
      submat = "BLOSUM62"
    )@score
  )

seq_scores_col <-
  cbind("water_score" = seqalign) %>%
  as.data.frame() %>%
  mutate(water_score = as.numeric(water_score))
seq_scores_col %>% glimpse

descscales <- fastas_merged %>%
  purrr::set_names(names(.)) %>% 
  purrr::map(
    ~ extractDescScales(.,
      propmat = "AAGETAWAY",
      pc = 3, lag = 5, silent = F
    )
  )

# Visualizing PCA components ----
pca_df <- descscales %>%
  bind_rows(.id = "id") %>%
  distinct() %>%
  column_to_rownames("id") %>%
  cbind(seq_scores_col) %>% 
  rownames_to_column("id") %>%
  mutate(dataset_origin = strex::str_before_last(id, "\\."))
pca_df %>% dim
pca_df %>% glimpse()
# View(pca_df)

# Plot ----
library(plotly)

fig_plotly <-
  pca_df %>% 
  filter(water_score > 100) %>% 
  plot_ly(
    x = ~scl1.lag3,
    y = ~scl2.lag3,
    z = ~scl3.lag3,
    # size = 0.8,
    text = ~ paste("ID:", id, "Similarity Score:", water_score)
  )

add_markers(
    fig_plotly, color = ~dataset_origin
  )



#' randomly select 5 sequences across the entire dataset calculate the 
set.seed(42)
rand_samples <- sample(pca_df$id, 25)
remaining_samples <- setdiff(pca_df$id, rand_samples)

pca_df_iter <-
  pca_df %>%
  mutate(
    selected_sample = ifelse(id %in% rand_samples, "Yes", "No")
  )

best_sd <- pca_df_iter %>%
  filter(selected_sample == "Yes") %>%
  select(contains("scl")) %>%
  # select(scl1.lag3, scl2.lag3) %>% 
  dplyr::summarise_all(~ sd(.)) %>%
  sum()

calculate_scl_sd <- function(df) {
  df %>%
    filter(selected_sample == "Yes") %>%
    select(contains("scl")) %>%
    # select(scl1.lag3, scl2.lag3) %>% 
    dplyr::summarise_all(~ sd(.)) %>%
    sum()
}

# pdf(
#   file = glue("test-selection-algo.pdf"),
#   width = 5, height = 5
# )
iter_max <- 2000
set.seed(2023)
pb <- progress::progress_bar$new(
  total = iter_max,
  format = "[:bar] :current/:total :percent eta: :eta"
)
iter_scores <- tibble(iter = 1:iter_max, sd_value = NA)
p_improvements <- list()
selected_hits <- list()
for (iter in 1:iter_max) {
  pb$tick()
  rand_samples <- sample(pca_df$id, 25)
  pca_df_iter <-
    pca_df %>%
    mutate(
      selected_sample = ifelse(id %in% rand_samples, "Yes", "No")
    )

  iter_scores[iter, "sd_value"] <- calculate_scl_sd(pca_df_iter)
  if (iter_scores[iter, "sd_value"] > best_sd) {
    selected_hits <- rand_samples
    best_sd <- iter_scores[iter, "sd_value"]
    p_improvements[[iter]] <- pca_df_iter %>%
      ggplot(aes(scl1.lag3, scl2.lag3)) +
      geom_point(aes(color = selected_sample, size = selected_sample)) +
      scale_color_d3() +
      scale_size_manual(values = c(0.5, 5)) +
      theme_light()
    # plot(p_improvements[[iter]])
  }
}
# dev.off()

# p_iter <- 
iter_scores %>%
  ggplot(aes(iter, sd_value)) +
  geom_point() +
  geom_line()


# p_improvements

fig_plotly_select <-
  pca_df %>%
  mutate(
    selected_sample = ifelse(id %in% selected_hits, "Yes", "No")
  ) %>% 
  plot_ly(
    x = ~scl1.lag3,
    y = ~scl2.lag3,
    z = ~scl3.lag3,
    # size = 0.8,
    text = ~ paste("ID:", id, "Similarity Score:", water_score)
  )

add_markers(
    fig_plotly_select, 
    color = ~selected_sample,
    symbol =  ~ factor(selected_sample),
    symbols = c("circle-open")
  )


# # selection refinement ----
# iter_max_refine <- 10000

# iter_n <- 1
# while( iter_n < iter_max_refine) {
#   iter_n <- iter_n + 1
# }


 
# average distance between all points ----


# ggsave(
#   filename = glue("{proj_dir}/structural_analysis/figures/fastas_round2_seqid_pca.png"),
#   plot = p,
#   width = 5, height = 5
# )

length(descscales)
descscales %>% glimpse

# library("Biostrings")
# library("foreach")
# library("doParallel")
# library(tictoc)

# tic()
# # extractPAAC(fasta_list[1])
#   psimmat <- protr::parSeqSim(
#     fasta_list,  batches  = 2, cores = 6, 
#     type = "local", submat = "BLOSUM62", verbose = TRUE
#   )
# toc()

# tic()
# psim <- protr::parSeqSimDisk(
#   fasta_list,  batches  = 1, cores = 6, 
#   type = "local", submat = "BLOSUM62", verbose = TRUE
# )
# toc()

pdf(
  file = glue("{proj_dir}/fastas_round2_seqid.pdf"),
  width = 7, height = 7
)
for (input_path in fasta_paths) {
  message("Plotting: ", input_path)
  s1 <- protr::readFASTA(input_path)
  names(s1) <- paste0("neo2_", 1:length(s1))
  fasta_list <- sample(s1, 100)
  fasta_list %<>% protr::removeGaps(pattern = "X", replacement = "G")


  psimmat <- protr::parSeqSim(
    fasta_list, cores = 6, type = "local", submat = "BLOSUM62"
  )
  colnames(psimmat) <- names(fasta_list)
  rownames(psimmat) <- names(fasta_list)

  feature_order <-
    dist(x = psimmat, method = "euclidean") %>%
    seriate(method = "R2E") %>%
    get_order()
  ranked_features_seqid <- rownames(psimmat)[feature_order]

  psi_df <- psimmat %>%
    as.data.frame() %>%
    rownames_to_column(var = "seq1") %>%
    pivot_longer(
      -seq1,
      names_to = "seq2",
      values_to = "similarity"
    ) %>%
    mutate(
      seq1 = factor(seq1, ordered = TRUE, levels = ranked_features_seqid),
      seq2 = factor(seq2, ordered = TRUE, levels = ranked_features_seqid)
    )

  p_heatmap <- psi_df  %>%
    ggplot(aes(x = seq1, y = seq2, fill = similarity)) +
    geom_tile() +
    theme_light() +
    labs(fill = "Sequence Identity (%)", x = NULL, y = NULL) +
    scale_fill_viridis_c(option = "F", limits = c(0, 1)) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      legend.key.width = unit(1, 'cm')
    )
    plot(p_heatmap)
}
dev.off()

# input_path <- glue("{proj_dir}/fastas_round2/output_t001_030.fasta")
s1 <- protr::readFASTA(input_path)
names(s1) <- paste0("neo2_", 1:length(s1))
fasta_list <- sample(s1, 100)

psimmat <- protr::parSeqSim(
  fasta_list, cores = 6, type = "local", submat = "BLOSUM62"
)
colnames(psimmat) <- names(fasta_list)
rownames(psimmat) <- names(fasta_list)

feature_order <-
  dist(x = psimmat, method = "euclidean") %>%
  seriate(method = "R2E") %>%
  get_order()
ranked_features_seqid <- rownames(psimmat)[feature_order]

psi_df <- psimmat %>%
  as.data.frame() %>%
  rownames_to_column(var = "seq1") %>%
  pivot_longer(
    -seq1,
    names_to = "seq2",
    values_to = "similarity"
  ) %>%
  mutate(
    seq1 = factor(seq1, ordered = TRUE, levels = ranked_features_seqid),
    seq2 = factor(seq2, ordered = TRUE, levels = ranked_features_seqid)
  )

psi_df  %>%
  ggplot(aes(x = seq1, y = seq2, fill = similarity)) +
  geom_tile() +
  theme_light() +
  labs(fill = "Sequence Identity (%)", x = NULL, y = NULL) +
  scale_fill_viridis_c(option = "F") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.key.width = unit(1, 'cm')
  )

# names(tst) <- paste0("sampled_seq_", 1:l



# amino acid sequences that contain gaps ("-")
# aaseq <- list(
#   "MHGDTPTLHEYMLDLQPETTDLYCYEQLSDSSEXXXXXXXXXXEEDEIDGPAGQAEPDRAHYNIVTFCCKCDSTLRLCVQS",
#   "MHGDTPTLHEYMLDLQPETTDLYCYEQLNDSSEXXXXXXXXXXEEDEIDGPAGQAEPDRAHYNIVTFCCKCDSTLRLCVQS"
# )

# #' # gaps create issues for alignment
# parSeqSim(aaseq)

# # remove the gaps
# nogapseq <- protr::removeGaps(aaseq, pattern = "X", replacement = "G")
# parSeqSim(nogapseq)

# copy all rank 0 PDB files into a new dir and download
for (complex in names(af_pdb_paths)){
  newfilename <- glue("{complex}_rank_0.pdb")
  print(newfilename)
  shell_do(glue("cp {af_pdb_paths[[complex]]} {analysis_dir}/rank_0_PDBs/{newfilename}"))
}


analysis_dir <- glue("{proj_dir}/structural_analysis")
af_pdb_paths <- list.files(
  "/central/scratch/jbok/alphafold-multimer/2023-03-26_neo-IL2",
  full.names = TRUE, recursive = TRUE, pattern = "ranked_0.pdb"
)
names(af_pdb_paths) <- gsub("IL2RBG_", "", basename(dirname(af_pdb_paths)))
# af_models <- readRDS(glue("{analysis_dir}/af_r3dmol_models.rds"))
pdb_list <- readRDS(glue("{analysis_dir}/af_ranked-0-PDBs.rds"))

pdbs <- pdbaln(
  files = pdb_list,
  outfile = glue("{analysis_dir}/aln.fa")
)
glimpse(pdbs)

# Structure Superposition
# core <- core.find(pdbs)
# core.inds <- print(core, vol=1.0)
# col=rep("black", length(core$volume))
# col[core$volume<2]="pink"; col[core$volume<1]="red"
# plot(core, col=col)
# write.pdb(
#   xyz = pdbs$xyz[1, core.inds$xyz],
#   file = glue("{analysis_dir}/core.pdb")
# )



pdbs$id <- names(af_pdb_paths)
xyz <- pdbfit(pdbs)
xyz <- pdbfit(pdbs, outpath = "fitted_IL2")

glimpse(xyz)

# fiting the structures
xyz <- pdbfit(
  pdbs,
  # inds = core.inds, 
  outpath = glue("{analysis_dir}/fitted")
)
# xyz <- bio3d::fit.xyz(
#   fixed = pdbs$xyz[1, ],
#   mobile = pdbs,
#   fixed.inds = gaps$f.inds,
#   mobile.inds = gaps$f.inds,
#   outpath = "rough_fit",
#   full.pdbs = TRUE,
#   verbose=TRUE
# )

gaps <- gap.inspect(pdbs$xyz)



# Root mean square deviation (RMSD):
rd <- rmsd(xyz)
hist(rd, breaks=40, xlab="RMSD (Å)", main="Histogram of RMSD")

# RMSD clustering
hc.rd <- hclust(as.dist(rd))

hclustplot(hc.rd,
  k = 1,
  labels = pdbs$id,
  cex = 0.5,
  ylab = "RMSD (Å)",
  main = "RMSD Cluster Dendrogram", fillbox = FALSE
)

# # Root mean squared fluctuations (RMSF):
# # Ignore gap containing positions
# gaps.res <- gap.inspect(pdbs$ali)
# gaps.pos <- gap.inspect(pdbs$xyz)

# # Tailor the PDB structure to exclude gap positions for SSE annotation
# id <- grep("neo2-15alpha", pdbs$id)
# inds <- atom.select(pdb, resno = pdbs$resno[id, gaps.res$f.inds])
# ref.pdb <- trim.pdb(pdb, inds = inds)

# # Plot RMSF with SSE annotation and labeled with residue numbers (Figure 8.)
# rf <- rmsf(xyz[, gaps.pos$f.inds])
# plot.bio3d(rf, resno=ref.pdb, sse=ref.pdb, ylab="RMSF (Å)",
#            xlab="Residue No.", typ="l")



# Torsion/Dihedral analysis: ----
pdb <- af_pdb_paths %>%
  keep(grepl("neo2-15alpha", .)) %>%
  read.pdb()

tor <- torsion.pdb(pdb)
# Basic Ramachandran plot (Figure 9)
plot(tor$phi, tor$psi, xlab="phi", ylab="psi")

# Locate the two structures in pdbs
ind.a <- grep("neo2-15alpha", pdbs$id)
ind.b <- grep("neoMPNN10", pdbs$id)

# Exclude gaps in the two structures to make them comparable
gaps.xyz2 <- gap.inspect(pdbs$xyz[c(ind.a, ind.b), ])
a.xyz <- pdbs$xyz[ind.a, gaps.xyz2$f.inds]
b.xyz <- pdbs$xyz[ind.b, gaps.xyz2$f.inds]

# Compare CA based pseudo-torsion angles between the two structures
a <- torsion.xyz(a.xyz, atm.inc=1)
b <- torsion.xyz(b.xyz, atm.inc=1)
d.ab <- wrap.tor(a-b)
d.ab[is.na(d.ab)] <- 0

# Plot results with SSE annotation
plot.bio3d(abs(d.ab), resno=pdb, sse=pdb, typ="h", xlab="Residue No.", 
           ylab = "Difference Angle")



# Difference distance matrix analysis (DDM) ---
a <- dm.xyz(a.xyz)
b <- dm.xyz(b.xyz)
plot.dmat( (a - b), nlevels=10, grid.col="gray", xlab="neo2-15alpha", ylab="neoMPNN10")



# Principal Component Analysis (PCA) ----
pc.xray <- pca.xyz(xyz[, gaps.pos$f.inds])
pc.xray
plot(pc.xray)


