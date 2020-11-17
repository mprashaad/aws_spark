#Install the packages Seurat, here and feather (optional)

library(Seurat)
#Set the working directory to where this repository is stored on your computer
setwd("C:/Users/mpras/Desktop/GitHub")
library(here)

#Uncomment the line below if feather is to be used
#library(feather)

source(here("tabula-muris", "boilerplate.R"))
organs <- c('Aorta',
            'Bladder',
            'Brain_Myeloid',
            'Brain_Non-Myeloid',
            'Diaphragm',
            'Fat',
            'Heart',
            'Kidney',
            'Large_Intestine',
            'Limb_Muscle',
            'Liver',
            'Lung',
            'Mammary_Gland',
            'Marrow',
            'Pancreas',
            'Skin',
            'Spleen',
            'Thymus',
            'Tongue',
            'Trachea')
#length(organs)
for(i in 1:1){
  organ_of_interest = organs[i]
  # Load the per-plate metadata
  plate_metadata_filename = here('00_facs_raw_data', 'metadata_FACS.csv')
  
  plate_metadata <- read.csv(plate_metadata_filename, sep=",", header = TRUE)
  colnames(plate_metadata)[1] <- "plate.barcode"
  
  
  # Load the gene names and set the metadata columns by opening the first file
  filename = here('00_facs_raw_data', 'FACS', paste0(organ_of_interest, '-counts.csv'))
  
  raw.data = read.csv(filename, sep=",", row.names=1)
  
  plate.barcodes = lapply(colnames(raw.data), function(x) strsplit(strsplit(x, "_")[[1]][1], '.', fixed=TRUE)[[1]][2])
  
  barcode.df = t.data.frame(as.data.frame(plate.barcodes))
  
  rownames(barcode.df) = colnames(raw.data)
  colnames(barcode.df) = c('plate.barcode')
  
  rnames = row.names(barcode.df)
  meta.data <- merge(barcode.df, plate_metadata, by='plate.barcode', sort = F)
  row.names(meta.data) <- rnames
  
  # Sort cells by cell name
  meta.data = meta.data[order(rownames(meta.data)), ]
  raw.data = raw.data[, rownames(meta.data)]
  
  # Find ERCC's, compute the percent ERCC, and drop them from the raw data.
  erccs <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = TRUE)
  percent.ercc <- Matrix::colSums(raw.data[erccs, ])/Matrix::colSums(raw.data)
  ercc.index <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = FALSE)
  raw.data <- raw.data[-ercc.index,]
  
  # Create the Seurat object with all the data
  tiss <- CreateSeuratObject(counts = raw.data, meta.data = meta.data, project = organ_of_interest)
  tiss$percent.ercc <- percent.ercc
  
  # Change default name for sums of counts from nUMI to nReads
  colnames(tiss@meta.data)[colnames(tiss@meta.data) == 'nUMI'] <- 'nReads'
  
  # Create metadata columns for annotations
  tiss@meta.data[,'free_annotation'] <- NA
  tiss@meta.data[,'cell_ontology_class'] <- NA
  tiss <- subset(x = tiss, subset = nFeature_RNA > 500 & nCount_RNA > 50000)
  tiss <- NormalizeData(object = tiss, scale.factor = 1e6)

  tiss <- ScaleData(object = tiss)
  tiss_var <- FindVariableFeatures(object = tiss)
  test <- subset(tiss_var, features = VariableFeatures(tiss_var))
  
  #csv version
  mat <- as.matrix(GetAssayData(object = test, slot = "scale.data"))
  write.csv(mat,paste0("processed_data/",organ_of_interest,"-norm_featselect_scale.csv"))
  
  #feather version - comment csv version chunk and uncomment below if feather is to be output format
  #mat <- as.data.frame(GetAssayData(object = test, slot = "scale.data"))
  #write_feather(mat, paste0("processed_data/",organ_of_interest,"-norm_featselect_scale.feather"))
}

