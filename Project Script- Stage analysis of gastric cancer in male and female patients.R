# Load Important Packages
library("BiocManager")
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("EDASeq")
library("gplots")
library("sesameData")
library("SummarizedExperiment")

# Get an overview on  LGG data
getProjectSummary("TCGA-STAD")

# Query the LGG Data
STAD.query <- GDCquery(project = "TCGA-STAD",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       barcode = "STAD.metaData_Gender")


# Download the dataset
GDCdownload(STAD.query)


# Prepare the dataset
STAD.Data <- GDCprepare(STAD.query, summarizedExperiment = TRUE)

# Check the output of the data for STAD Data
head(STAD.Data)
View(STAD.Data)

# Explore the data
STAD.Data$gender
STAD.Data$barcode
STAD.Data$patient
STAD.DAta$
  

STAD.metaData_Gender <- data.frame("Gender" = STAD.Data$gender, 
                            "Stage" = STAD.Data$paper_TNM.Stage,
                            "Barcode" = STAD.Data$barcode)


# Extract the raw data from the prepared dataset
STAD.Raw <- assays(STAD.Data)

#Select unstranded data
dim(STAD.Raw$unstranded)
View(STAD.Raw$unstranded)

# Group the stages into early and late stage
Early_stage <- c("Stage_IA", "Stage_IB", "Stage_IIA", "Stage_IIB")
Late_stage <- c("Stage_III", "Stage_IIIB", "Stage_IIIC", "Stage_IV")

# Subset the data
SelectedBarcodes_Gender <- c(subset(STAD.metaData_Gender, Gender == "female" & Stage %in% Early_stage)$Barcode,
                             subset(STAD.metaData_Gender, Gender == "female" & Stage %in% Late_stage)$Barcode,
                             subset(STAD.metaData_Gender, Gender == "male" & Stage %in% Early_stage)$Barcode,
                             subset(STAD.metaData_Gender, Gender == "male" & Stage %in% Late_stage)$Barcode)
                             

#Retrieve the unstranded data for the selected barcodes
Selected_STAD_GenderData <- STAD.Raw$unstranded[, c(SelectedBarcodes_Gender)]
dim(Selected_STAD_GenderData)
View(Selected_STAD_GenderData)

#Normalize Data
STAD.Gender_normalized <- TCGAanalyze_Normalization(tabDF = Selected_STAD_GenderData, geneInfo = geneInfoHT, method = "geneLength")

# Then Filter
STAD.Gender_filtered <- TCGAanalyze_Filtering(tabDF = STAD.Gender_normalized,
                                              method = "quantile",
                                              qnt.cut = 0.25)

View(STAD.Gender_filtered)
dim(STAD.Gender_filtered)

# Create matrices for Female samples
mat1_Female_Early_stage <- STAD.Gender_filtered[, colnames(STAD.Gender_filtered) %in% subset(STAD.metaData_Gender, Gender == "female" & Stage %in% earlystage)$Barcode]
mat2_Female_Late_stage <- STAD.Gender_filtered[, colnames(STAD.Gender_filtered) %in% subset(STAD.metaData_Gender, Gender == "female" & Stage %in% latestage)$Barcode]

# Create matrices for Male samples
mat1_Male_Early_stage <- STAD.Gender_filtered[, colnames(STAD.Gender_filtered) %in% subset(STAD.metaData_Gender, Gender == "male" & Stage %in% earlystage)$Barcode]
mat2_Male_Late_stage <- STAD.Gender_filtered[, colnames(STAD.Gender_filtered) %in% subset(STAD.metaData_Gender, Gender == "male" & Stage %in% latestage)$Barcode]

# Perform DEA for female Samples (early stages vs late stages)
results_female_Samples <- TCGAanalyze_DEA(mat1 = mat1_Female_Early_stage,
                                          mat2 = mat2_Female_Late_stage,
                                          Cond1type = "female Stage %in% Early_stage",
                                          Cond2type = "female Stage %in% Late_stage",
                                          pipeline = "edgeR",
                                          fdr.cut = 0.01,
                                          logFC.cut = 2)

# Perform DEA for male Samples (early stages vs late stages)
results_male_Samples <- TCGAanalyze_DEA(mat1 = mat1_Male_Early_stage,
                                          mat2 = mat2_Male_Late_stage,
                                          Cond1type = "male Stage %in% Early_stage",
                                          Cond2type = "male Stage %in% Late_stage",
                                          pipeline = "edgeR",
                                          fdr.cut = 0.01,
                                          logFC.cut = 2)

#Volcano Plot for Female samples
FEMALE_volcano_plot <- plot(results_female_Samples$logFC, -log10(
  results_female_Samples$FDR))

#Volcano Plot for Male samples
MALE_volcano_plot <- plot(results_male_Samples$logFC, -log10(
  results_male_Samples$FDR))

#DEA with expression levels for female samples
female_Expression.level <- TCGAanalyze_LevelTab(results_female_Samples,
                                                "female Stage %in% Early_stage",
                                                "female Stage %in% Late_stage"
                                                , mat1_Female_Early_stage, mat2_Female_Late_stage)

#DEA with expression levels for male samples
male_Expression.level <- TCGAanalyze_LevelTab(results_male_Samples,
                                                "male Stage %in% Early_stage",
                                                "male Stage %in% Late_stage"
                                                , mat1_Male_Early_stage, mat2_Male_Late_stage)

head() # To view expression result. 
head(female_Expression.level)
head(male_Expression.level)

# Create vectors for the barcodes of female samples
Female_Early_Stage_barcodes <- subset(STAD.metaData_Gender, Gender == "female" & Stage %in% Early_stage)$Barcode
Female_Late_Stage_barcodes <- subset(STAD.metaData_Gender, Gender == "female" & Stage %in% Late_stage)$Barcode

# Create vectors for the barcodes of male samples
Male_Early_Stage_barcodes <- subset(STAD.metaData_Gender, Gender == "male" & Stage %in% Early_stage)$Barcode
Male_Late_Stage_barcodes <- subset(STAD.metaData_Gender, Gender == "male" & Stage %in% Late_stage)$Barcode

# Combine the barcodes for heatmap data for female samples
selected_barcodes_Female_samples <- c(Female_Early_Stage_barcodes,
                                          Female_Late_Stage_barcodes)

# Combine the barcodes for heatmap data for male samples
selected_barcodes_Male_samples <- c(Male_Early_Stage_barcodes,
                                    Male_Late_Stage_barcodes)

# Subset the heatmap data for female samples
heat.data.Female <- STAD.Gender_filtered[rownames(female_Expression.level), 
                                             selected_barcodes_Female_samples]


# Subset the heatmap data for male samples
heat.data.Male <- STAD.Gender_filtered[rownames(male_Expression.level), 
                                         selected_barcodes_Male_samples]

# Create a color vector for the Female samples
sample_colors_female <- c(rep("midnightblue", 55), rep("red4", 30)) 

heatmap.2(x = as.matrix(heat.data.Female),
          col = hcl.colors(10, palette = 'Blue-Red 3'),  
          Rowv = FALSE, Colv = TRUE,
          scale = 'row',
          sepcolor = 'black',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.5, cexCol = 1,
          main = "Heatmap of Female Samples, Early Stage vs Late Stage",
          na.color = 'black',
          ColSideColors = sample_colors_female)

#Add a Legend
legend("topright", legend = c("Early stage", "Late stage"),
       fill = c("midnightblue", "red4"),
       border = "black",  # Black outline
       bty = "n")  # No box around the legend

dev.off()

# Create a color vector for the Male samples
sample_colors_male <- c(rep("purple3", 87), rep("gold", 52)) 

heatmap.2(x = as.matrix(heat.data.Male),
          col = hcl.colors(10, palette = 'Blue-Red 3'),  
          Rowv = FALSE, Colv = TRUE,
          scale = 'row',
          sepcolor = 'black',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.5, cexCol = 1,
          main = "Heatmap of Male Samples, Early Stage vs Late Stage",
          na.color = 'black',
          ColSideColors = sample_colors_male)

#Add a Legend
legend("topright", legend = c("Early stage", "Late stage"),
       fill = c("purple", "gold"),
       border = "black",  # Black outline
       bty = "n")  # No box around the legend

dev.off()

# Upregulated and downregulated genes for Female samples
Female_upreg_early_stage <- rownames(subset(results_female_Samples, logFC > 2 & "female Stage %in% Early_stage" > 0))
Female_upreg_late_stage <- rownames(subset(results_female_Samples, logFC > 2 & "female Stage %in% Late_stage" > 0))
Female_downreg_early_stage <- rownames(subset(results_female_Samples, logFC < 2 & "female Stage %in% Early_stage" > 0))
Female_downreg_late_stage <- rownames(subset(results_female_Samples, logFC < 2 & "female Stage %in% Late_stage" > 0))


# Upregulated and downregulated genes for Male samples
Male_upreg_early_stage <- rownames(subset(results_male_Samples, logFC > 2 & "male Stage %in% Early_stage" > 0))
Male_upreg_late_stage <- rownames(subset(results_male_Samples, logFC > 2 & "male Stage %in% Late_stage" > 0))
Male_downreg_early_stage <- rownames(subset(results_male_Samples, logFC < 2 & "male Stage %in% Early_stage" > 0))
Male_downreg_late_stage <- rownames(subset(results_male_Samples, logFC < 2 & "male Stage %in% Late_stage" > 0))

library("biomaRt")
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Convert Ensemble IDs to gene symbols for female samples
upregulated_Early_Female_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Female_upreg_early_stage,
  mart = mart
)

upregulated_Late_Female_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Female_upreg_late_stage,
  mart = mart
)

downregulated_Early_Female_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Female_downreg_early_stage,
  mart = mart
)

downregulated_Late_Female_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Female_downreg_late_stage,
  mart = mart
)


# Convert Ensemble IDs to gene symbols for male samples
upregulated_Early_Male_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Male_upreg_early_stage,
  mart = mart
)

upregulated_Late_Male_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Male_upreg_late_stage,
  mart = mart
)

downregulated_Early_Male_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Male_downreg_early_stage,
  mart = mart
)

downregulated_Late_Male_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Male_downreg_late_stage,
  mart = mart
)

# Functional Enrichment Analysis

up.EA.female_earlyStage <- TCGAanalyze_EAcomplete(TFname = "female Stage %in% Early_stage", RegulonList = upregulated_Early_Female_symbols$hgnc_symbol)
up.EA.female_lateStage <- TCGAanalyze_EAcomplete(TFname = "female Stage %in% Late_stage", RegulonList = upregulated_Late_Female_symbols$hgnc_symbol)
up.EA.male_earlyStage <- TCGAanalyze_EAcomplete(TFname = "male Stage %in% Early_stage", RegulonList = upregulated_Early_Male_symbols$hgnc_symbol)
up.EA.male_lateStage <- TCGAanalyze_EAcomplete(TFname = "male Stage %in% Late_stage", RegulonList = upregulated_Late_Male_symbols$hgnc_symbol)
down.EA.female_earlyStage <- TCGAanalyze_EAcomplete(TFname = "female Stage %in% Early_stage (Downregulated)", RegulonList = downregulated_Early_Female_symbols$hgnc_symbol)
down.EA.female_lateStage <- TCGAanalyze_EAcomplete(TFname = "female Stage %in% Late_stage (Downregulated)", RegulonList = downregulated_Late_Female_symbols $hgnc_symbol)
down.EA.male_earlyStage <- TCGAanalyze_EAcomplete(TFname = "male Stage %in% Early_stage (Downregulated)", RegulonList = downregulated_Early_Male_symbols$hgnc_symbol)
down.EA.male_lateStage <- TCGAanalyze_EAcomplete(TFname = "male Stage %in% Late_stage (Downregulated)", RegulonList = downregulated_Late_Male_symbols$hgnc_symbol)

#Visualization of the enrichment analysis
UP.Female_earlyStage <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.female_earlyStage$ResBP),
                                        GOBPTab = up.EA.female_earlyStage$ResBP,
                                        GOCCTab = up.EA.female_earlyStage$ResCC,
                                        GOMFTab = up.EA.female_earlyStage$ResMF,
                                        PathTab = up.EA.female_earlyStage$ResPat,
                                        nRGTab = up.EA.female_earlyStage,
                                        nBar = 10,
                                        text.size = 2,
                                        fig.width = 30,
                                        fig.height = 15)

UP.Female_lateStage <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.female_lateStage$ResBP),
                                                GOBPTab = up.EA.female_lateStage$ResBP,
                                                GOCCTab = up.EA.female_lateStage$ResCC,
                                                GOMFTab = up.EA.female_lateStage$ResMF,
                                                PathTab = up.EA.female_lateStage$ResPat,
                                                nRGTab = up.EA.female_lateStage,
                                                nBar = 10,
                                                text.size = 2,
                                                fig.width = 30,
                                                fig.height = 15)


DOWN.Female_earlyStage <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.female_earlyStage$ResBP),
                                                GOBPTab = down.EA.female_earlyStage$ResBP,
                                                GOCCTab = down.EA.female_earlyStage$ResCC,
                                                GOMFTab = down.EA.female_earlyStage$ResMF,
                                                PathTab = down.EA.female_earlyStage$ResPat,
                                                nRGTab = down.EA.female_earlyStage,
                                                nBar = 10,
                                                text.size = 2,
                                                fig.width = 30,
                                                fig.height = 15)

DOWN.Female_lateStage <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.female_lateStage$ResBP),
                                               GOBPTab = down.EA.female_lateStage$ResBP,
                                               GOCCTab = down.EA.female_lateStage$ResCC,
                                               GOMFTab = down.EA.female_lateStage$ResMF,
                                               PathTab = down.EA.female_lateStage$ResPat,
                                               nRGTab = down.EA.female_lateStage,
                                               nBar = 10,
                                               text.size = 2,
                                               fig.width = 30,
                                               fig.height = 15)



UP.Male_lateStage <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.male_lateStage$ResBP),
                                                GOBPTab = up.EA.male_lateStage$ResBP,
                                                GOCCTab = up.EA.male_lateStage$ResCC,
                                                GOMFTab = up.EA.male_lateStage$ResMF,
                                                PathTab = up.EA.male_lateStage$ResPat,
                                                nRGTab = up.EA.male_lateStage,
                                                nBar = 10,
                                                text.size = 2,
                                                fig.width = 30,
                                                fig.height = 15)


UP.Male_earlyStage <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.male_earlyStage$ResBP),
                                                GOBPTab = up.EA.male_earlyStage$ResBP,
                                                GOCCTab = up.EA.male_earlyStage$ResCC,
                                                GOMFTab = up.EA.male_earlyStage$ResMF,
                                                PathTab = up.EA.male_earlyStage$ResPat,
                                                nRGTab = up.EA.male_earlyStage,
                                                nBar = 10,
                                                text.size = 2,
                                                fig.width = 30,
                                                fig.height = 15)

DOWN.Male_earlyStage <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.male_earlyStage$ResBP),
                                                  GOBPTab = down.EA.male_earlyStage$ResBP,
                                                  GOCCTab = down.EA.male_earlyStage$ResCC,
                                                  GOMFTab = down.EA.male_earlyStage$ResMF,
                                                  PathTab = down.EA.male_earlyStage$ResPat,
                                                  nRGTab = down.EA.male_earlyStage,
                                                  nBar = 10,
                                                  text.size = 2,
                                                  fig.width = 30,
                                                  fig.height = 15)


DOWN.Male_lateStage <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.male_lateStage$ResBP),
                                                 GOBPTab = down.EA.male_lateStage$ResBP,
                                                 GOCCTab = down.EA.male_lateStage$ResCC,
                                                 GOMFTab = down.EA.male_lateStage$ResMF,
                                                 PathTab = down.EA.male_lateStage$ResPat,
                                                 nRGTab = down.EA.male_lateStage,
                                                 nBar = 10,
                                                 text.size = 2,
                                                 fig.width = 30,
                                                 fig.height = 15)
