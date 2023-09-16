
#' Function to load and trim PDB files to a specific chain
load_trimmed_pdb <- function(filepath, chain_of_interest) {
  pdb <- bio3d::read.pdb(filepath)
  if (!is.na(chain_of_interest)) {
    pdb <- bio3d::trim.pdb(
      pdb,
      bio3d::atom.select(pdb, chain = chain_of_interest)
    )
  }
  return(pdb)
}

# #' Function to extract the amino acid sequence from a PDB file
# get_pdb_seq <- function(pdb) {
#   require(magrittr)
#   ca_inds <- bio3d::atom.select(pdb, "calpha")
#   seq_aa <- bio3d::pdbseq(pdb, inds = ca_inds, aa1 = TRUE) %>%
#       unname() %>%
#       paste0(collapse = "")
#   return(seq_aa)
# }

#' This function converts a PDB file with only c-alpha atoms
#' into a PDB file with all atoms using pulchra (installed in path)
ca_to_all <- function(pdb_path, wk = wkdir) {
    require(glue)
    source(glue("{wk}/notebooks/R-scripts/helpers_general.R"))
    # load in the pdb file
    pdb <- bio3d::read.pdb(pdb_path)
    # save pdb in a temp dir to reformat
    tmp_dir <- glue("{wk}{tempdir()}")
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
    tmp_pdb <- glue("{tmp_dir}/{basename(pdb_path)}")
    bio3d::write.pdb(pdb, file = tmp_pdb)
    # run pulchra command to reconstruct atomisitic
    shell_do(glue("pulchra {tmp_pdb}"))
    tmp_pdb_fixed <- glue(
        "{tmp_dir}/{strex::str_before_first(basename(pdb_path), '.pdb')}.",
        "rebuilt.pdb"
    )
    pdb_all <- bio3d::read.pdb(tmp_pdb_fixed)
    unlink(tmp_pdb)
    unlink(tmp_pdb_fixed)
    return(pdb_all)
}

#' Function to overlay multiple PDB files in a 3D viewer
overlay_foldseek_pdbs <- function(
    pdb_path_list,
    base_pdb = NULL,
    ca_inputs = TRUE,
    align = FALSE, # only works with base_model input
    max = 5) {
  require(r3dmol)
  require(magrittr)
  coln <- if (length(pdb_path_list) > max) {
    max
  } else {
    length(pdb_path_list)
  }
  color_pal <-
    suppressWarnings(RColorBrewer::brewer.pal(coln, "Oranges"))
  base_model <- r3dmol(
    viewer_spec = m_viewer_spec(
      cartoonQuality = 10,
      lowerZoomLimit = 50,
      upperZoomLimit = 750
    )
  )
  if (is.null(base_pdb)) {
    model_cnt <- 1
  } else {
    base_model %<>%
      m_add_model(data = m_bio3d(base_pdb)) %>%
      m_zoom_to() %>%
      m_set_style(
        sel = m_sel(model = 0),
        style = m_style_cartoon(color = "#1F77B4")
      )
    model_cnt <- 0
  }
  for (i in 1:length(pdb_path_list)) {
    if (i <= max) {
      if (ca_inputs) {
        pdb <- ca_to_all(pdb_path_list[i])
      } else {
        pdb <- bio3d::read.pdb(pdb_path_list[i])
      }
      if (align && !is.null(base_pdb)) {
        pdbs_aln <- bio3d::struct.aln(
          base_pdb,
          pdb,
          exefile = "msa",
          max.cycles = 100,
          outpath = paste0("fitlsq", tempdir()),
          write.pdbs = FALSE,
          verbose = FALSE
        )
        pdb$xyz <- pdbs_aln$xyz
      }
      base_model %<>%
        m_zoom_to() %>%
        m_add_model(data = m_bio3d(pdb)) %>%
        m_set_style(
          sel = m_sel(model = i - model_cnt),
          style = m_style_cartoon(color = color_pal[i])
        )
    }
  }
  return(base_model)
}

#' function to plot a base model (a specific chain) and overlay the target hits
view_foldseek_pdb <- function(
    query_pdb_path,
    query_chain = NA,
    ca_inputs = TRUE,
    align = FALSE,
    target_pdbs_path_list) {
  if (is.na(query_chain)) {
    query_pdb <- bio3d::read.pdb(query_pdb_path)
  } else {
    query_pdb <- load_trimmed_pdb(
      query_pdb_path, query_chain
    )
  }
  final_models <- overlay_foldseek_pdbs(
    pdb_path_list = target_pdbs_path_list,
    base_pdb = query_pdb,
    ca_inputs = ca_inputs,
    align = align
  )
  return(final_models)
}

#' This is a function to download PDBs from
#' Deepmind's Alphafold database or
#' MGnify's ESM database
download_pdb <- function(name, outpath) {
  require(glue)
  require(httr)
  # Example names: AF-Q7M0U6-F1-model_v4 or MGYP003670600000
  if (grepl("AF-", name)) {
    response <- GET(
      glue("https://alphafold.ebi.ac.uk/files/{name}.pdb")
    )
  } else if (grepl("MGYP", name)) {
    response <- GET(
      glue("https://api.esmatlas.com/fetchPredictedStructure/{name}")
    )
  } else {
    stop("name must be Alphafold or MGnify ESM PDB")
  }
  # downloading contents
  if (response$status_code == 200) {
    pdb <- content(response, as = "text")
    write.table(
      pdb,
      file = glue("{outpath}/{name}.pdb"),
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE
    )
  } else {
    print(paste0("Error downloading PDB file: ", response$status_code))
  }
}

format_foldseek_pdb_paths <- function(pdb_paths){
  pdb_df <- data.table("pdb_path" = pdb_paths) %>%
    mutate(
      gene_target_id =
        str_before_nth(pdb_path, "/", -2) %>% str_after_last(., "/"),
      filename = str_after_nth(pdb_path, "/", -2) %>%
        str_after_last(., "/") %>%
        str_after_first("_"),
      db = case_when(
        grepl("Alphafold_Proteome", filename) ~
          "Alphafold_Proteome",
        grepl("foldseek_Alphafold_SwissProt", filename) ~
          "foldseek_Alphafold_SwissProt",
        grepl("foldseek_PDB", filename) ~
          "foldseek_PDB",
        grepl("foldseek_ESMAtlas30", filename) ~
          "foldseek_ESMAtlas30",
        grepl("foldseek_Alphafold_UniProt50", filename) ~
          "foldseek_Alphafold_UniProt50",
        TRUE ~ "ERROR"
      ),
      res_pair = strex::str_after_first(string = filename, pattern = db),
      target = case_when(
        grepl("MGYP", res_pair) ~ str_after_last(res_pair, "_"),
        grepl("AF-", res_pair) ~ str_after_nth(res_pair, "_", -2),
        db == "foldseek_PDB" ~ str_after_nth(res_pair, "_", -2),
        TRUE ~ "ERROR"
      ),
      query = str_remove(res_pair, target),
      query = case_when(
        grepl("MODEL", query) ~ gsub("_MODEL_\\d+", "", query),
        TRUE ~ query
      ),
      query = str_before_last(query, "_"),
      target = str_remove(target, "\\.pdb$")
    )
  return(pdb_df)
}

#' Function to extract fasta sequence from bio3d pdb object
extract_fasta_from_pdb <- function(pdb) {
  bio3d::pdbseq(
    pdb,
    inds = bio3d::atom.select(pdb, "calpha"), aa1 = TRUE
  ) %>%
    unname() %>%
    paste0(collapse = "")
}

#' Function to overlay multiple PDB files in a 3D viewer
overlay_pdbs <- function(
    pdb_path_list,
    base_model = NULL,
    max = 5) {
  require(r3dmol)
  coln <- if (length(pdb_path_list) > max) {
    max
  } else {
    length(pdb_path_list)
  }
  color_pal <-
    suppressWarnings(RColorBrewer::brewer.pal(coln, "Oranges"))
  if (is.null(base_model)) {
    base_model <- r3dmol(
      viewer_spec = m_viewer_spec(
        cartoonQuality = 10,
        lowerZoomLimit = 50,
        upperZoomLimit = 750
      )
    )
    model_cnt <- 1
  } else {
    model_cnt <- 0
  }
  for (i in 1:length(pdb_path_list)) {
    if (i <= max) {
      pdb <- ca_to_all(pdb_path_list[i])
      base_model %<>%
        m_zoom_to() %>%
        m_add_model(data = m_bio3d(pdb)) %>%
        m_set_style(
          sel = m_sel(model = i - model_cnt),
          style = m_style_cartoon(color = color_pal[i])
        )
    }
  }
  return(base_model)
}