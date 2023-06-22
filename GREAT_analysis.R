library(devtools)
library(rGREAT)
library(wordcloud)
library(dplyr)
library(ggplot2)
library(cowplot)
library(magick)
library(scales)
library(grid)
library(gridExtra)
library(openxlsx)

## Usage:
# Rscript GREAT_analysis.R dir test.bed background.bed out.pdf out.xlsx
#   dir - directory to use in analysis. By default, input and output files should be here
#   test.bed - 3 column bed file with test regions
#   background.bed - 3 column bed file with background regions
#   out.pdf - name of output plot - should be pdf
#   out.xlsx - name of excel file that all data will be saved to

args = commandArgs(trailingOnly=TRUE)
dir=args[1]
infile=args[2]
bkgd_file=args[3]
outplot=args[4]
outfile=args[5]

grid.newpage()
pushViewport(viewport(layout = grid.layout(length(analyses), 3)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

## set up plot - we will be plotting a word cloud of nearby genes and barplot of top 10 biological processes
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 3)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

ontologies=data.frame(full=c("Ensembl Genes",
                             "GO Molecular Function",
                             "GO Biological Process",
                             "GO Cellular Component",
                             "Human Phenotype",
                             "Mouse Phenotype",
                             "Mouse Phenotype Single KO"),
                      short=c("Genes","GO MF","GO BP","GO CC","HP","MP","MP KO"))

setwd(dir)

i=1
wb <- createWorkbook()

print(infile)

test_regions <- read.table(infile, header=F, sep='\t')
format(test_regions, scientific=F)
bkgd_regions <- read.table(bkgd_file, header=F, sep='\t')
format(bkgd_regions, scientific=F)

job = submitGreatJob(test_regions, bg = bkgd_regions, species='hg38', version='4.0.4', 
                     includeCuratedRegDoms = TRUE,
                     rule='basalPlusExt', adv_upstream = 5.0, adv_downstream = 1.0, adv_span = 1000.0)

tbl = getEnrichmentTables(job, ontology=ontologies$full)

for (ontology in ontologies$full) {
  print(ontology)
  df <- as.data.frame(tbl[[ontology]])
  if (ontology %in% c("GO Molecular Function","GO Biological Process","GO Cellular Component")) {
    df <- df[df$Hyper_Fold_Enrichment > 2,]
  }
  if (ontology=='Ensembl Genes') {
    png("tmp.png", width = 1000, height = 1000, unit = "px")
    wordcloud(words = df$name, 
              freq = -log10(df$Hyper_Adjp_BH),
              min.freq = 1,  
              max.words=400,
              random.order=FALSE,
              rot.per=0.35,
              colors=brewer.pal(8, "Dark2"))
    dev.off()
    p1 <- ggdraw() + draw_image("tmp.png") + ggtitle(paste(analysis, ontology))
    print(p1, vp = vplayout(1, 1))
  } else if (ontology=='GO Biological Process') {
    p2 <- ggplot(df[1:10,], aes(x=factor(name, levels=df$name), y=-log10(Hyper_Adjp_BH))) + 
      geom_bar(stat="identity",fill='black') + coord_flip() +
      theme_cowplot() + theme(axis.text.x = element_text(size = 8),
                              axis.text.y = element_text(size = 8),
                              axis.title.x = element_text(size = 10),
                              axis.title.y = element_text(size = 10)) + 
      scale_x_discrete(labels = label_wrap(65)) +
      ggtitle(paste(infile)) + 
      ylab('-log10(FDR Q-value)') + xlab('Enriched Biological Processes')
    print(p2, vp = vplayout(1, 2:3))
  }
  sheetname = infile
  addWorksheet(wb = wb, sheetName = sheetname, gridLines = FALSE)
  writeDataTable(wb = wb, sheet = sheetname, x = df)
}

saveWorkbook(wb, outfile, overwrite = TRUE)
ggsave(outplot, width=11, height=3)


