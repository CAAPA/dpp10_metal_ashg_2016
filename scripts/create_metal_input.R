#Create old bdos METAL input file
###############################################################################
in.frame <- read.delim("../data/raw/old_results/chr2_rsIDs_pos_assoc_file_DPP10_BDOS_mlinfo_rsq_ge_03.txt")

in.frame <- in.frame[,c(
  "Marker", #MARKER
  "Al2",  #REF
  "Al1", #ALT
  "MAF", 
  "Freq1.x", "Freq0", #Used to derive EFFECT
  "PValue",  #PVALUE
  "Chromosome",
  "Start"
)]

in.frame$EFFECT <- "+"
in.frame$EFFECT[in.frame$Freq0 > in.frame$Freq1.x] <- "-"

out.frame <- in.frame[,c(1:4,10,7:9)]
names(out.frame) <- c("MARKER", "REF", "ALT", "MAF", "EFFECT", "PVALUE", "CHR", "BP")
write.table(out.frame, "../data/input/bdos_old.txt",  sep="\t", quote=F, row.names=F, col.names=T)

#Create old graad METAL input file
###############################################################################
in.frame <- read.delim("../data/raw/old_results/chr2_rsIDs_pos_assoc_file_DPP10_GRAAD_mlinfo_rsq_ge_03.txt")

in.frame <- in.frame[,c(
  "Marker", #MARKER
  "Al2",  #REF
  "Al1", #ALT
  "MAF", 
  "t", #Used to derive EFFECT
  "p.val",  #PVALUE
  "Chromosome",
  "Start"
)]

in.frame$EFFECT <- "+"
in.frame$EFFECT[in.frame$t < 0] <- "-"

out.frame <- in.frame[,c(1:4,9,6:8)]
names(out.frame) <- c("MARKER", "REF", "ALT", "MAF", "EFFECT", "PVALUE", "CHR", "BP")
write.table(out.frame, "../data/input/graad_old.txt",  sep="\t", quote=F, row.names=F, col.names=T)

#Create new bdos METAL input file
###############################################################################
in.frame <- read.delim("../data/raw/new_results/2_output_mqlsdose_k_0.20.txt")

in.frame <- in.frame[,c(
  "Marker", #MARKER
  "Freq1", "Freq0", #Used to derive EFFECT
  "PValue"  #PVALUE
)]
in.frame <- in.frame[!is.na(in.frame$PValue),]

info.frame <- read.delim("../data/raw/new_results/2.mach.info")[,c(1,2,3,5,7)]
info.frame <- info.frame[info.frame$MAF >= 0.01,]
info.frame <- info.frame[info.frame$Rsq > 0.5,]
in.frame <- merge(in.frame, info.frame, by.x="Marker", by.y="SNP")

in.frame$Marker <- as.character(in.frame$Marker)
in.frame$CHR <- as.numeric(unlist(strsplit(in.frame$Marker, split=":"))[seq(1,dim(in.frame)[1]*2,2)])
in.frame$BP <- as.numeric(unlist(strsplit(in.frame$Marker, split=":"))[seq(2,dim(in.frame)[1]*2,2)])
in.frame <- in.frame[(in.frame$BP >= 115189899) & (in.frame$BP <= 116613328),]


in.frame$EFFECT <- "+"
in.frame$EFFECT[in.frame$Freq0 > in.frame$Freq1] <- "-"

out.frame <- in.frame[,c(1,5,6,7,10,4,8,9)]
names(out.frame) <- c("MARKER", "REF", "ALT", "MAF", "EFFECT", "PVALUE", "CHR", "BP")
write.table(out.frame, "../data/input/bdos_new.txt",  sep="\t", quote=F, row.names=F, col.names=T)

#Create new graad METAL input file
###############################################################################
in.frame <- read.delim("../data/raw/new_results/chr2_typed_overlap_allchohort.ps", head=F)
names(in.frame) <- c("MARKER", "BETA", "PVALUE")

info.frame <- read.table(gzfile("../data/raw/new_results/chr2.info.gz"), head=T, stringsAsFactors = F)[,c(1:3,5)]
in.frame <- merge(in.frame, info.frame, by.x="MARKER", by.y="SNP")

in.frame$MARKER <- as.character(in.frame$MARKER)
in.frame$CHR <- as.numeric(unlist(strsplit(in.frame$MARKER, split=":"))[seq(1,dim(in.frame)[1]*2,2)])
in.frame$BP <- as.numeric(unlist(strsplit(in.frame$MARKER, split=":"))[seq(2,dim(in.frame)[1]*2,2)])
in.frame <- in.frame[(in.frame$BP >= 115189899) & (in.frame$BP <= 116613328),]

in.frame$EFFECT <- "+"
in.frame$EFFECT[in.frame$BETA < 0] <- "-"

out.frame <- in.frame[,c(1,4,5,6,9,3,7,8)]
names(out.frame) <- c("MARKER", "REF", "ALT", "MAF", "EFFECT", "PVALUE", "CHR", "BP")
write.table(out.frame, "../data/input/graad_new.txt",  sep="\t", quote=F, row.names=F, col.names=T)

