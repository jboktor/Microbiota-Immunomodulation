
attach(transducin)

pdb <- read.pdb("1tag")
seq <- pdbseq(pdb)
blast <- blast.pdb(seq)

# See Figure 2.
hits <- plot.blast(blast, cutoff=240)
anno <- pdb.annotate(hits)
head(anno[, c("resolution", "ligandId", "citation")])

# Download PDBs and split by chain ID
files <- get.pdb(hits, path="raw_pdbs", split = TRUE)
# Extract and align sequences
pdbs <- pdbaln(files)


core <- core.find(pdbs)
# See Figure 3.
col=rep("black", length(core$volume))
col[core$volume<2]="pink"; col[core$volume<1]="red"
plot(core, col=col)
core.inds <- print(core, vol=1.0)
# write.pdb(xyz=pdbs$xyz[1,core.inds$xyz], file="quick_core.pdb")


# xyz <- pdbfit(pdbs, core.inds)
xyz <- pdbfit(pdbs)
xyz <- pdbfit(pdbs, outpath="fitted_v2")



# #files <- get.pdb(c("4q21","5p21"), URLonly=TRUE)
# files <- get.pdb(c("4q21","5p21"), path=tempdir(), overwrite=TRUE)
# pdbs <- pdbaln(files)
# xyz <- pdbfit(pdbs)
# # Superpose again this time outputing all-atom PDBs to disc
# xyz <- pdbfit( pdbs, outpath="fitted" )

