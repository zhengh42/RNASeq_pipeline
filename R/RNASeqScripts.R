require(DESeq2)
DESeq2Analysis=function(data,coldata,type,normmethod){
  if(type=="tximport"){
    dds <- DESeqDataSetFromTximport(data,colData = coldata,design= ~condition)
  }
  if(type=="count"){
    dds <- DESeqDataSetFromMatrix(data,colData = coldata,design= ~condition)
  }
  
  
  #dds <- dds[ rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)
  if(normmethod=="rld"){countnorm <- rlog(dds, blind=FALSE)}
  if(normmethod=="vst"){countnorm <- vst(dds, blind=FALSE)}

  ### top genes
  #select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:1000]
  #p1<- Heatmap(assay(countnorm)[select,],
  #      clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
  #      clustering_method_rows  = "average",clustering_method_columns  = "average",
  #      show_row_names =F,
  #      name="distance",
  #      column_names_gp = gpar(fontsize = 5),row_names_gp = gpar(fontsize = 5))
  
  #draw(p1,newpage = T)

  ### sample clustering
  #sampleDists <- dist(t(assay(countnorm)))
  #sampleDistMatrix <- as.matrix(sampleDists)
  #colnames(sampleDistMatrix) <- NULL

  #colors <- colorRampPalette( rev(brewer.pal(9, "Oranges")) )(255)
  #p2<- Heatmap(sampleDistMatrix,
  #      clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
  #      clustering_method_rows  = "average",clustering_method_columns  = "average",
  #      col=colors,
  #      name="Distance",
  #      column_names_gp = gpar(fontsize = 5),row_names_gp = gpar(fontsize = 5)
  #           )
  #draw(p2,newpage = T)
  
  ### sample PCA

  pcaData <- plotPCA(countnorm, intgroup=c("condition"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  p3=ggplot(pcaData, aes(PC1, PC2, color=condition)) +
    theme_bw()+
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() +
     theme(text = element_text(color="#525252"),panel.grid = element_blank(),panel.border = element_blank(),
         legend.title=element_blank(), 
          axis.line = element_line(color = 'grey'),axis.ticks=element_line(color = 'grey')) 
    #scale_color_manual(values = c("case"="#ffb7c5","case2"="#b7c5ff","ctrl"="#b7e9ff"))

  res <- results(dds)
  #resLFC <- lfcShrink(dds, coef=2)

  #p4=plotMA(res, ylim=c(-2,2))
  #p5=plotMA(resLFC, ylim=c(-2,2))

  #DESeq2.gene.out <- as.data.frame(res)
  #DESeq2.gene.out$log2FoldChange_lfc <- resLFC$log2FoldChange
  #DESeq2.gene.out$lfcSE_lfc <- resLFC$lfcSE
  return(list("dds"=dds,"countnorm"=countnorm,"pcaplot"=p3))
}

