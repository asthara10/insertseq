#!/usr/bin/env Rscript
library(annotatr)
library(dplyr)
library(ggplot2)
library(scales)

args = commandArgs(trailingOnly=TRUE)
bedFile <- args[1]
outFile <- args[2]
outFile2 <- args[3]
outFile3 <- args[4]
pre <- args[5]

# Build the annotations (a single GRanges object)
annots = c("hg38_basicgenes", "hg38_genes_cds", "hg38_cpg_islands",
           "hg38_genes_intergenic")
annotations = build_annotations(genome = 'hg38', annotations = annots)

# Read insertions bed file
df <- read.table(bedFile)
names(df) <- c("chr", "start", "end")
regions <- GRanges(df)

# Intersect the regions we read in with the annotations
annotated = annotate_regions(
  regions = regions,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)
annotated <- as.data.frame(annotated)

# Collapse by symbol
annotated$insertion <- paste(paste(annotated$seqnames, annotated$start, sep=":"), annotated$end, sep="-")
formated <- as.data.frame(annotated %>% dplyr::group_by(insertion) %>% dplyr::distinct(annot.symbol))
formated <- formated[!is.na(formated$annot.symbol),]
write.table(formated, outFile, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

# Collapse by insertion

annotated$insertion <- paste(paste(annotated$seqnames, annotated$start, sep=":"), annotated$end, sep="-")
type.short <- lapply(annotated$annot.type, function(x) strsplit(x, "_")[[1]][3])
select <- type.short %in% c("introns", "exons", "intergenic", "3UTRs", "5UTRs", "promoters", "CpG")
type.filtered <- annotated[select,]
type.filtered$region.type <- type.short[select]
formated.type <- as.data.frame(type.filtered %>% dplyr::group_by(insertion) %>% dplyr::distinct(region.type))
formated.type <- formated.type[!is.na(formated.type$region.type),]
formated.type$region.type <- as.character(formated.type$region.type)
write.table(formated.type, outFile2, quote=FALSE, sep="\t", row.name=FALSE, col.names=FALSE)

types <- as.character(formated.type$region.type)
plot.df <- data.frame("region" = names(table(types)), "count" = as.vector(table(types)))

pie <- ggplot(plot.df, aes(x="", y=count/sum(count)*100, fill=region)) +
  geom_bar(width = 0.5, stat = "identity")  +
  theme_classic() + scale_fill_brewer(palette="Set2") +
  labs(x="", y="% insertions at gene region") +
  theme(panel.border = element_blank(),
        panel.grid=element_blank(), axis.ticks = element_blank(),
        axis.text.x=element_blank(), axis.text.y = element_text(size=rel(1.8)), axis.title=element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(1.2)), legend.title = element_blank(), legend.position = "top",
        strip.text.x = element_text(size = rel(1.8)),
        aspect.ratio = 1.5/1) 

ggsave(outFile3, pie)
