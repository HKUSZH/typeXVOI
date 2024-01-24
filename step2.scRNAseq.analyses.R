library(Seurat)
library(dplyr)
library(patchwork)
load("../resources/TF.CD.receptors.human.mouse.expanded.RData")
source("../my-func-lib/WRITE_seurate_to_web.R")
source("../my-func-lib/helpers.GETSIG.R")

source("../my-func-lib/gsea.visualization.R")

############################################

load("../data/WNT1.2.3.4.5.RData")

load("../data/LLD2.3.4.5.RData")
load("../data/LLD.RData")

UMAPPlot(WNT5)
pdf("tsne.plot.wnt5.lld5.pdf")
	TSNEPlot(WNT3, label=T, label.size=8, pt.size=4)
	TSNEPlot(LLD5, label=T, label.size=8, pt.size=4)

	TSNEPlot(WNT5, label=T, label.size=8, pt.size=4)
	TSNEPlot(LLD5, label=T, label.size=8, pt.size=4)
	PCAPlot(LLD5, label=T, label.size=8, group.by="newPops", pt.size=4)
	PCAPlot(WNT5, label=T, label.size=8, group.by="newPops", pt.size=4)

dev.off()

#################################
load("../data/Signature.genes.of.the.8clusters.in.P6.RData")
LSIG.mouse<-GETSIGm(MARK.cellPopsName)

sapply(LSIG.mouse, function(x){
	sapply(LSIG.WNT5, function(y){
		xy<-intersect(toupper(x), y)
		lenxy<-length(xy)
		round(lenxy/length(x)/length(y)*2e4, 1)
	})
})

names(LSIG.LLD5)<-c("L2-imOb", "L1-prog", "L3-mOb", "L5-periosteum", "L4-Osteocyte")
sapply(LSIG.mouse, function(x){
	sapply(LSIG.LLD5, function(y){
		xy<-intersect(toupper(x), y)
		lenxy<-length(xy)
		mat1<-rbind(c(2e4, length(x)), c(length(y), lenxy))
		res1<-chisq.test(mat1)
		OR<-round(lenxy/length(x)/length(y)*2e4, 1)
		paste0(OR, ";", signif(res1$p.value, 2))
	})
})
sapply(LSIG.mouse, function(x){
	sapply(LSIG.LLD5, function(y){
		xy<-intersect(toupper(x), y)
		lenxy<-length(xy)
	})
})

vec2<-as.vector(sapply(LSIG.mouse, function(x){
	xy<-intersect(toupper(x), LSIG.LLD5[[3]])
	paste0(xy, collapse=", ")
}))

vec3<-as.vector(sapply(LSIG.mouse, function(x){
	xy<-intersect(toupper(x), LSIG.WNT5[[3]])
	paste0(xy, collapse=", ")
}))

vec4<-as.vector(sapply(LSIG.mouse, function(x){
	xy<-intersect(toupper(x), LSIG.WNT5[[4]])
	paste0(xy, collapse=", ")
}))
data.frame(vec3, names(LSIG.mouse))
data.frame(vec4, names(LSIG.mouse))

#################################
WNT.genes<-read.table(file="WNT.genes", sep="\t")

intersect(LSIG.WNT5[[3]], WNT.genes[,2])
intersect(LSIG.WNT5[[4]], WNT.genes[,2])

#################################
UMAPPlot(LLD5, label=T, label.size=8)
pdf("tsne.W3.W4.pdf", width=8)
	WNT[["inWNT5"]]<- "0"
	WNT[["inWNT5"]][colnames(WNT) %in% colnames(WNT5)[WNT5[["newPops"]][,1]%in%c("W1-imOb/prog")], 1]<-"W1-imOb/prog"
	WNT[["inWNT5"]][colnames(WNT) %in% colnames(WNT5)[WNT5[["newPops"]][,1]%in%c("W2-mOb")], 1]<-"W2-mOb"
	WNT[["inWNT5"]][colnames(WNT) %in% colnames(WNT5)[WNT5[["newPops"]][,1]%in%c("W3")], 1]<-"W3"
	WNT[["inWNT5"]][colnames(WNT) %in% colnames(WNT5)[WNT5[["newPops"]][,1]%in%c("W4")], 1]<-"W4"
	WNT[["inWNT5"]]<-factor(WNT[["inWNT5"]][,1], sort(unique(WNT[["inWNT5"]][,1])))

	TSNEPlot(WNT, label=T, pt.size=2, label.size=8, group.by="inWNT5")
	Idents(WNT)<-"inWNT5"
	FeaturePlot(WNT, features=c("RUNX2", "ADIPOQ", "PPARG", "CXCR4"), cols=c("grey", "brown"), reduction="tsne")
	FeaturePlot(WNT5, features=c("RUNX2", "ADIPOQ", "PPARG", "CXCR4"), cols=c("grey", "brown"), reduction="tsne")
	FeaturePlot(WNT5, features=c("RUNX2", "MCAM"), cols=c("grey", "brown"), reduction="tsne")

	FeaturePlot(WNT5, features=c("RUNX2", "ALPL", "SFRP2", "GDF10"), cols=c("grey", "brown"), reduction="tsne")
	FeaturePlot(WNT5, features=c("COL10A1", "STEAP4", "SFRP2", "GDF10"), cols=c("grey", "brown"), reduction="tsne")

	

dev.off()


pdf("violin.0.W1.W2.W3.W4.pdf", height=6, width=12)
	VlnPlot(WNT, features=c("RUNX2", "COL1A1", "COL10A1", 
		"STEAP4", "SFRP1", "SFRP2","FZD4", "GDF10", "LTBP1", "COL4A1", "COL4A2", 
		"PTH1R", "CD44"), pt.size=0, ncol=6)
dev.off()


	FeaturePlot(LLD, features=c("COL10A1"), cols=c("grey", "brown"), pt.size=2, reduction="tsne")


UMAPPlot(WNT, label=T, label.size=8, group.by="inWNT5")
FeaturePlot(WNT, features=c("RUNX2", "ALPL"))
FeaturePlot(WNT, features=c("RUNX2", as.vector(sapply(LSIG.WNT5[3:4], head, 10))), ncol=8)

FeaturePlot(LLD, features=c("RUNX2", as.vector(sapply(LSIG.WNT5[3:4], head, 10))), ncol=8)


FeaturePlot(WNT, features=c("ADIPOQ", "LPL"))
WNT[["oneB"]]<-factor(WNT@meta.data$oneBased, level=sort(unique(WNT@meta.data$oneBased)))
LLD[["oneB"]]<-factor(LLD@meta.data$oneBased, level=sort(unique(LLD@meta.data$oneBased)))

VlnPlot(WNT5, features=c("RUNX2", "COL1A1", "COL2A1", "SOX9", "ALPL", "IBSP", "SPP1"))

pdf("heatmap.WNT.pdf", height=6, width=12)
	DoHeatmap(WNT, features=as.vector(sapply(LSIG.WNT, head, 10)), group.by="oneB")
	DoHeatmap(LLD, features=as.vector(sapply(LSIG.LLD, head, 10)), group.by="oneB")
dev.off()

FeaturePlot(WNT5, features=c("RUNX2"))

sapply(LSIG.WNT5, function(x){
	sapply(LSIG.LLD5, function(y){
		paste0(c(round(length(intersect(x,y))/length(x)/length(y)*1e4, 2), length(x), length(y), length(intersect(x,y))), collapse=",")
	})
})

sapply(LSIG.WNT5, head, 20)
sapply(LSIG.LLD5, head, 20)


FeaturePlot(WNT, features=c("RUNX2", "SOX9", "TAGLN", "MYH11", "TOP2A", "SP7", "TNFRSF11A"))
FeaturePlot(LLD, features=c("RUNX2", "SP7", "TNFRSF11B", "TNFRSF11A"), reduction="tsne")

LLD5[["newPops"]]<-c("L2-imOb", "L1-prog", "L3-mOb", "L5-periosteum", "L4-Osteocyte")[as.integer(LLD5@meta.data$seurat_clusters)]
table(LLD5@meta.data$newPops, LLD5@meta.data$seurat_clusters)

WNT5[["newPops"]]<-c("W2-mOb", "W1-imOb/prog", "W3-AberrAdi1", "W4-AberrAdi2")[as.integer(WNT5@meta.data$seurat_clusters)]
WNT5[["newPopsNum"]]<-c("W2", "W1", "W3", "W4")[as.integer(WNT5@meta.data$seurat_clusters)]
table(WNT5@meta.data$newPops, WNT5@meta.data$seurat_clusters)


pdf("dot.plot.wnt5.lld5.pdf", height=2.5, width=10)
	DotPlot(LLD5, features=c("RUNX2", "COL1A1", "COL2A1", "SOX9", 
		"COL4A1", "COL4A2", "LEPR", "LPL", "IGFBP4", "CXCL12",
		"ALPL", "FGFR3", "DKK3", "GREM1", "IGFBP5", "WIF1", "SHOX2", "HIF1A",
		"IBSP", "BGLAP", "CREB3L1", "MEPE", "CADM1", "IFITM5", "MMP13", "PTH1R",
		"CD44", "DMP1", "PHEX", "PDPN", "TNFRSF11B", "DKK1", "IRX5", 
		"CILP2", "ASPN", "FAP", "POSTN"), group.by="newPops", cols=c("grey", "brown"))

	DotPlot(WNT5, features=c("RUNX2","COL1A1", "COL1A2", "SOX9", 
		"COL4A1", "COL4A2", "LEPR", "LPL", "IGFBP4", "CXCL12",
		"ALPL", "FGFR3", "DKK3", "GREM1", "IGFBP5", "WIF1", "SHOX2", "HIF1A",
		"IBSP", "BGLAP", "CREB3L1", "MEPE", "CADM1", "IFITM5", "MMP13", "PTH1R",
		"CD44", "DMP1", "PHEX", "PDPN", "TNFRSF11B", "DKK1", "IRX5", 
		"CILP2", "ASPN", "FAP", "POSTN"), group.by="newPops", cols=c("grey", "brown"))

	Idents(WNT5)<-"newPops"
	FeaturePlot(WNT5, features=c("ADIRF", "MEDAG", "COL4A1","COL4A2", "LTBP4", "GDF10", "STEAP4", "LEPR", "CXCL12", "PLIN2", "PPARG", "ADIPOQ",
		"RUNX2", "COL1A1", "ASPN", "POSTN"),
		label=T, reduction="tsne", cols=c("grey", "brown"))

	Idents(WNT5)<-"newPopsNum"
	FeaturePlot(WNT5, features=c("ADIRF", "PPARG", "ADIPOQ", "SFRP2", "SFRP1", "RUNX2", "COL1A1", "SOX9", "PLIN2"),
		label=T, reduction="tsne", cols=c("grey", "brown"))

	FeaturePlot(WNT5, features=c("ADIRF", "CRLF1", "ADIPOQ", "SFRP2", "SFRP1", "RUNX2", "COL1A1", "SOX9", "PLIN2"),
		label=T, reduction="tsne", cols=c("grey", "brown"))

	FeaturePlot(WNT5, features=c("ADIRF", "AHR", "IL6", "SFRP2", "SFRP1", "RUNX2", "COL1A1", "SOX9", "PLIN2"),
		label=T, reduction="tsne", cols=c("grey", "brown"))

	FeaturePlot(WNT5, features=c("ADIRF", "GSK3B", "CTNNB1", "SFRP2", "ACP2", "COX17", "LYPLA2", "MIF", "TOP3B"),
		label=T, reduction="tsne", cols=c("grey", "brown"))

	DotPlot(WNT5, features=as.vector(sapply(LSIG.WNT5[c(2,1,3,4)], head, 9)), group.by="newPops", cols=c("grey", "brown"))
dev.off()

pdf("feature.plot.wnt5.adipo.pdf", height=8, width=10)
	FeaturePlot(WNT5, features=c("ADIRF", "CRLF1", "ADIPOQ", "SFRP2", "SFRP1", "RUNX2", "COL1A1", "SOX9", "PLIN2"),
		label=T, reduction="pca", cols=c("grey", "brown"))

	FeaturePlot(WNT5, features=c("AHR", "IL6", "ADIPOQ", "SFRP4", "LPL", "IGF1", "COL1A1", "SOX9", "PLIN2"),
		label=T, reduction="pca", cols=c("grey", "brown"))
dev.off()

FeaturePlot(WNT5, features=c("RUNX2"), pt.size=4, label=T)
FeaturePlot(LLD5, features=c("RUNX2"), pt.size=4, label=T)

table(as.matrix(GetAssayData(WNT5)["RUNX2", ])>0)
table(as.matrix(GetAssayData(LLD5)["RUNX2", ])>0)

table(as.matrix(GetAssayData(WNT5)["BGLAP", ])>0)
table(as.matrix(GetAssayData(LLD5)["BGLAP", ])>0)

VlnPlot(LLD5, features=c("SOST", "PHEX", "PDPN", "BGLAP"))
VlnPlot(LLD5, features=c("SOST", "PHEX", "PDPN", "BGLAP"))

VlnPlot(LLD5, features=c("PRRX1"))
VlnPlot(WNT5, features=c("PRRX1"))
############################################
WNT6<-subset(WNT5, ident=c(0,1))

WNT6<-NormalizeData(WNT6)
WNT6<-ScaleData(WNT6,features=rownames(WNT6))
WNT6<-FindVariableFeatures(WNT6)
WNT6<-RunPCA(WNT6)
WNT6<-FindNeighbors(WNT6, dims = 1:10)
WNT6<-FindClusters(WNT6, resolution = 0.25)
WNT6<-RunUMAP(WNT6, reduction= "pca", dims= 1:40)
WNT6<-RunTSNE(WNT6, reduction= "pca", dims.use = 1:40, do.fast = T)
UMAPPlot(WNT6, label=T,label.size=6, pt.size = 4) 
UMAPPlot(WNT6, group.by="oneBased", label=T,label.size=6, pt.size = 1) 
VlnPlot(WNT6, features=c("nFeature_RNA", "FRZB", "ASPN", "CXCL12", "LEPR", "GREM1", "STEAP4", 
	"COL1A1", "RUNX2", "HBB", "COL2A1", "SOX9",  "IFITM5", "BGLAP", "PHEX"), pt.size=0)
#LSIG.WNT2<-GETSIG(WNT2)

############################################
WNT_LLD55<-RunCCA(LLD5,WNT5)
WNT_LLD55[["OldClust"]]<-paste0(WNT_LLD55@meta.data$orig.ident,"_",as.integer(as.character(WNT_LLD55@meta.data$seurat_clusters))+1)
table(WNT_LLD55[["OldClust"]][,1], WNT_LLD55@meta.data$orig.ident)

VlnPlot(WNT_LLD55, features=c("nFeature_RNA"), group.by="OldClust")

WNT_LLD<-RunCCA(LLD5,WNT6)

WNT_LLD[["OldClust"]]<-paste0(WNT_LLD@meta.data$orig.ident,"_",as.integer(as.character(WNT_LLD@meta.data$seurat_clusters))+1)
table(WNT_LLD[["OldClust"]][,1], WNT_LLD@meta.data$orig.ident)


WNT_LLD<-NormalizeData(WNT_LLD)
WNT_LLD<-ScaleData(WNT_LLD, fearures=rownames(WNT_LLD))
WNT_LLD<-FindVariableFeatures(WNT_LLD)
WNT_LLD<-RunPCA(WNT_LLD)
WNT_LLD<-RunUMAP(WNT_LLD, reduction= "cca", dims= 1:10)
WNT_LLD<-RunTSNE(WNT_LLD, reduction= "cca", dims.use = 1:40, do.fast = T)
WNT_LLD<-FindNeighbors(WNT_LLD, reduction= "cca", dims = 1:10)
WNT_LLD<-FindClusters(WNT_LLD, resolution = 0.25)
WNT_LLD.5<-FindClusters(WNT_LLD, resolution = 0.5)

levels(WNT_LLD.5@meta.data$seurat_clusters)<-c("2", "1", "5", "3", "4")
WNT_LLD.5[["newClust"]]<-paste0(as.integer(as.character(WNT_LLD.5@meta.data$seurat_clusters)), "_", WNT_LLD.5@meta.data$orig.ident)
WNT_LLD.5[["newClust2"]]<-as.integer(as.character(WNT_LLD.5@meta.data$seurat_clusters))

tab56<-table(as.integer(as.character(WNT_LLD.5@meta.data$seurat_clusters)), WNT_LLD.5@meta.data$orig.ident)

LSIG56<-GETSIG(WNT_LLD.5)[c(2,1,4,5,3)]

WNT_LLD.5[["combinedClust"]]<-factor(as.integer(WNT_LLD.5@meta.data$seurat_clusters), levels=c(2,1,4,5,3))

pdf("heatmap.Wnt.LLD.5.combinedClust.pdf")
	DoHeatmap(WNT_LLD.5, features=as.vector(sapply(LSIG56, head, 10)), group.by="combinedClust")
dev.off()


pdf("piechart.cell.composition.Wnt6.LLD5.pdf")
	ggplot(as.data.frame(tab56)%>%group_by(Var2)%>%
		mutate(prop=100*Freq/sum(Freq),lab=rev(paste0(round(prop,1),"%;n=",Freq)),
		ypos=cumsum(rev(prop))-rev(prop)/2),
			aes(x=1,y=prop,label=lab,fill=Var1)) + 
		geom_col(width=1)  + xlim(c(0,1.5)) + geom_text(aes(y=ypos)) +
		facet_wrap(vars(Var2)) +
		scale_fill_brewer(palette="Set2") +
		coord_polar("y") + ggtitle("deformity distr in the two clusters")

	ggplot(as.data.frame(tab56[1:4, ])%>%group_by(Var2)%>%
		mutate(prop=100*Freq/sum(Freq),lab=rev(paste0(round(prop,1),"%;n=",Freq)),
		ypos=cumsum(rev(prop))-rev(prop)/2),
			aes(x=1,y=prop,label=lab,fill=Var1)) + 
		geom_col(width=1)  + xlim(c(0,1.5)) + geom_text(aes(y=ypos)) +
		facet_wrap(vars(Var2)) +
		scale_fill_brewer(palette="Set2") +
		coord_polar("y") + ggtitle("deformity distr in the two clusters")
dev.off()

table(WNT_LLD.5[["newClust"]][,1], WNT_LLD.5@meta.data$newPops)

pie(apply(tab56[1:4, ], 2, function(x)x/sum(x))[,1])

pdf("dot.plot.WNT6.LLD5.combined.pdf", height=5, width=10)

	DotPlot(WNT_LLD.5, features=c("RUNX2", "COL1A1", "COL2A1", "SOX9", 
		"COL4A1", "COL4A2", "LEPR", "LPL", "IGFBP4", "CXCL12",
		"ALPL", "FGFR3", "DKK3", "GREM1", "IGFBP5", "WIF1", "SHOX2", "HIF1A",
		"IBSP", "BGLAP", "CREB3L1", "MEPE", "CADM1", "IFITM5", "MMP13", "PTH1R",
		"CD44", "DMP1", "PHEX", "PDPN", "TNFRSF11B", "DKK1", "IRX5", 
		"CILP2", "ASPN", "FAP", "POSTN"), group.by="newClust", cols=c("grey", "brown"))

	DotPlot(WNT_LLD.5, features=as.vector(sapply(LSIG56, head, 8)), 
		group.by="newClust", cols=c("grey", "brown"))

dev.off()



table(WNT_LLD[["OldClust"]][,1], Idents(WNT_LLD))
table(WNT_LLD.5[["OldClust"]][,1], Idents(WNT_LLD.5))

TSNEPlot(WNT_LLD, label=T, label.size=6, pt.size=4, group.by="orig.ident")
PCAPlot(WNT_LLD, label=T, label.size=6, pt.size=4, group.by="orig.ident")

#WNT_LLD@meta.data$orig.ident[WNT_LLD@meta.data$orig.ident=="zLLD"]<-"LLD"

pdf("tsne.LLD5.WNT6.pdf", height=5, width=6)
	TSNEPlot(WNT_LLD.5, label=T, label.size=4, pt.size=2, group.by="newClust", shape.by="orig.ident")
	TSNEPlot(WNT_LLD.5, label=T, label.size=8, pt.size=2, group.by="newClust2", shape.by="orig.ident")
	TSNEPlot(WNT_LLD.5, label=T, label.size=8, pt.size=2, group.by="newClust2", split.by="orig.ident", shape.by="orig.ident")
	TSNEPlot(WNT_LLD.5, label=T, label.size=8, pt.size=2, group.by="newClust2", shape.by="orig.ident")
	TSNEPlot(WNT_LLD, label=T, label.size=4, pt.size=2, group.by="newClust", shape.by="orig.ident")
	PCAPlot(WNT_LLD, label=T, label.size=4, pt.size=2, group.by="OldClust", shape.by="orig.ident")
dev.off()

table(WNT_LLD.5@meta.data$seurat_clusters, WNT_LLD.5@meta.data$OldClust)

GETDEG<-function(GRP1="LLD_2", GRP2="WNT_2", GRPBY="OldClust"){
	DEGtmp<-FindMarkers(WNT_LLD.5, ident.1 = GRP1,
		ident.2 = GRP2, group.by = GRPBY,
		only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
	print(table(DEGtmp$p_val_adj<0.05, abs(DEGtmp$avg_log2FC)>1))
	return(DEGtmp)
}
GETDEG2<-function(DEGtmp){
	DEGtmp2<-DEGtmp[which(DEGtmp$p_val_adj<0.05 & abs(DEGtmp$avg_log2FC)>1), ]
	DEGtmp2
}

DEG_prog<-GETDEG(GRP1="1_WNT", GRP2="1_LLD", GRPBY="newClust")
DEG_imOb<-GETDEG(GRP1="2_WNT", GRP2="2_LLD", GRPBY="newClust")
DEG_mOb<-GETDEG(GRP1="3_WNT", GRP2="3_LLD", GRPBY="newClust")
DEG_perio<-GETDEG(GRP1="5_WNT", GRP2="5_LLD", GRPBY="newClust")
DEG_WL<-GETDEG(GRP1="WNT", GRP2="LLD", GRPBY="orig.ident")

DEG_prog2<-GETDEG2(DEG_prog)
DEG_imOb2<-GETDEG2(DEG_imOb)
DEG_mOb2<-GETDEG2(DEG_mOb)
DEG_perio2<-GETDEG2(DEG_perio)
DEG_WL2<-GETDEG2(DEG_WL)

########################################
hBMSC<-read.delim("../data/RNAseq/gene_exp.diff", sep="\t")
DEGh<-hBMSC %>% filter(significant=="yes") %>% arrange(log2.fold_change.)

########################################
library(gplots)
LosteoUp<-list(prog2=rownames(DEG_prog2)[DEG_prog2$avg_log2FC>0],
	imOb2=rownames(DEG_imOb2)[DEG_imOb2$avg_log2FC>0],
	mOb2=rownames(DEG_mOb2)[DEG_mOb2$avg_log2FC>0])
	#WL2=rownames(DEG_WL2)[DEG_WL2$avg_log2FC>0])
LosteoDown<-list(prog2=rownames(DEG_prog2)[DEG_prog2$avg_log2FC<0],
	imOb2=rownames(DEG_imOb2)[DEG_imOb2$avg_log2FC<0],
	mOb2=rownames(DEG_mOb2)[DEG_mOb2$avg_log2FC<0])
	#WL2=rownames(DEG_WL2)[DEG_WL2$avg_log2FC<0])

LosteoUp2<-c(LosteoUp, list(perio2=rownames(DEG_perio2)[DEG_perio2$avg_log2FC>0]))
LosteoDown2<-c(LosteoDown, list(perio2=rownames(DEG_perio2)[DEG_perio2$avg_log2FC<0]))

sapply(split(DEGh$gene, sign(DEGh$log2.fold_change.)), function(x){
	intersect(x, unique(unlist(DEG_WL2$genes)))
})

pdf("venn.wnt.ctrl.up.down.degs.pdf")
	vup<-venn(LosteoUp)
	vdown<-venn(LosteoDown)

	vup2<-venn(LosteoUp2)
	vdown2<-venn(LosteoDown2)
dev.off()

attr(vup, "intersections")
attr(vdown, "intersections")

attr(vup2, "intersections")
attr(vdown2, "intersections")

pdf("Pathways.up.and.down.pdf", height=16)
	Canonical.Up.allWnt<-get.gsea("Canonical-up.tsv", db="", num.shown=70)
	Canonical.Down.allWnt<-get.gsea("Canonical-down.tsv", db="", num.shown=70)

	Canonical.Up.allWnt[[1]][grep("WNT", names(Canonical.Up.allWnt[[1]]))]

	GO.Up.allWnt<-get.gsea("GO-up.tsv", db="", num.shown=100)
	GO.Down.allWnt<-get.gsea("GO-down.tsv", db="", num.shown=100)

	Canonical.Up.allWnt[[1]][grep("WNT", names(Canonical.Up.allWnt[[1]]))]

dev.off()


split(rownames(DEG_WL2), sign(DEG_WL2$avg_log2FC))

VOLCANO<-function(DEGtmp){
	DEGtmp[["genes"]]<-rownames(DEGtmp)
	DEGtmp[["logp"]]<- -log10(DEGtmp$p_val)
	DEGtmp[["status"]]<- c("lower", "mid", "higher")[sign(DEGtmp$avg_log2FC)+2]
	DEGtmp2<-DEGtmp[which(DEGtmp$p_val_adj<0.05 & abs(DEGtmp$avg_log2FC)>1), ]

	g1<-ggplot(DEGtmp, aes(x=avg_log2FC, y=logp, label=genes))+
		geom_point(aes(size=logp), color="grey")+
		geom_point(data=DEGtmp2, aes(size=logp, color=status))+
		geom_text(data=DEGtmp2, aes(size=logp))
	plot(g1)
	return(DEGtmp2)
}
pdf("WNT-vs-LLD.DEG.pdf", width=9)
	DEG_prog2<-VOLCANO(DEG_prog)
	DEG_imOb2<-VOLCANO(DEG_imOb)
	DEG_mOb2<-VOLCANO(DEG_mOb)
	DEG_perio2<-VOLCANO(DEG_perio)
	DEG_WL2<-VOLCANO(DEG_WL)
dev.off()


VlnPlot(WNT_LLD, features=c("IGFBP4"), group.by="OldClust")
VlnPlot(WNT_LLD, features=c("LPL"), group.by="OldClust")
VlnPlot(WNT_LLD, features=c("RUNX2", "ADIPOQ"), group.by="OldClust")


###############
proteomics<-readxl::read_excel("../proteomics_summary.compare.xlsx", sheet = "OI__Control")
DEPs<-proteomics[proteomics[["P-value"]]<0.05, ]
LDEPs<-sapply(split(DEPs$Gene, sign(DEPs[[8]])), sort)
LDEGs<-sapply(split(DEG_WL2$genes, sign(DEG_WL2[[2]])), sort)

sapply(LDEPs, function(x){
	sapply(LDEGs, function(y){
		paste0(intersect(x, y), collapse=", ")
	})
})

sapply(LDEPs, function(x){
	sapply(LSIG.LLD5, function(y){
		paste0(intersect(x, y), collapse=", ")
	})
})

colnames(proteomics)[1:13]<-c("Accession", "Name", "Gene", "Description", "MolkD",
	"Length", "OIControl", "log2OIControl", "Pvalue", "Padjust", "DiffSig", "Count", "OI1IFITM5")
proteomics[["logp"]]<- -log10(proteomics[["Pvalue"]])
proteomics[["status"]]<- c("zneg", "pos")[(sign(proteomics[["log2OIControl"]])+3)/2]
DEPs<-proteomics%>%filter(Pvalue<0.05)

g1<-ggplot(proteomics, aes(x=log2OIControl, y=logp, label=Gene))+
		geom_point(aes(size=logp), color="grey")+
		geom_point(data=DEPs, aes(size=logp, color=status))+
		geom_text(data=DEPs, aes(size=logp))
plot(g1)


