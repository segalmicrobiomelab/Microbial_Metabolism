#Segal Lab
#Programmer: Kendrew Wong
#Project: COPD functional 

#clear all objects
rm(list=ls(all.names=T))
gc()

##open up the R script with all of the function created for the project
source("KW_functions_v1k12.R")

#packages needed for the script
pkg_list <- c("devtools", "tidyverse", "readr", "readtext", "vegan", "ade4", "biomformat", "cachem", 
              "utf8", "backports", "colorspace", "rhdf5", "DelayedArray","Biobase", "Biostrings", "magick",
              "phyloseq", "ape", "phangorn", "ggpubr", "decontam", "ggplot2", "reshape2",
              "tidyr", "matrixStats", "DESeq2", "edgeR", "limma", "Glimma", "RColorBrewer",
              "pheatmap", "ggrepel", "scales", "data.table", "fBasics", "forcats", "maptools", 
              "lubridate", "boot", "table1", "stringr", "papaja", "flextable",
              "Gmisc", "glue", "htmlTable", "grid", "magrittr", "rmarkdown", "plotly",
              "microbiome", "DT", "webshot", "lubridate", "png", "RCurl", "cowplot", "janitor",
              "optmatch", "MatchIt", "decontam", "qdap", "stringr","openxlsx", "chisq.posthoc.test", 
              "FSA", "cobalt", "ggplotify", "grid", "gridExtra", "combinat", "mixOmics", "gplots", "plyr", 
              "readxl", "DESeq2", "mia", "microbiomeMarker", "jpeg", "openxlsx",
              "mia", "miaViz", "corrplot", "ggcorrplot", "cowplot", "gridGraphics",
              "ade4", "ggthemes", "Hmisc", "rdist", "rstatix", "ggpubr", "pdftools", "convertGraph",
              "qiime2R","metagMisc","YesSiR", "microViz","pulsar","patchwork","SpiecEasi")

#loading packages
for (i in pkg_list){
  eval(bquote(library(.(i))))
}

#function for MT MG data transformation:
MTMG_transform_data<- function (dataframe, 
                                meta_dataframe, 
                                meta_dataframe_var_interest=list(),
                                meta_dataframe_subjectID,
                                output_dataframe, 
                                first_column_name){
  #setting up the files so they are standardized
  #keeping the key variables from the meta file
  meta_dataframe_keyvariable<-meta_dataframe %>% select(meta_dataframe_var_interest)
  colnames_meta_dataframe_keyvariable<-colnames(meta_dataframe_keyvariable)
  meta_dataframe_keyvariable<- meta_dataframe_keyvariable[meta_dataframe_keyvariable[[meta_dataframe_subjectID]] !="", ] #remove the missing value 
  #setting up the dataframe
  rownames_dataframe<-dataframe[,1] #first column to be row names
  dataframe<-dataframe[,colnames(dataframe)!=first_column_name]
  dataframe <- as.data.frame(sapply(dataframe, as.numeric))
  rownames(dataframe)<-rownames_dataframe
  #turn the dataframe with each row being sample ID 
  dataframe_transpose<-data.frame(t(dataframe))
  dataframe_transpose$sampleID<-rownames(dataframe_transpose)
  #extract the unique sample name from the sample ID
  dataframe_transpose$Sample.ID_Segal.Lab<-sapply(strsplit(dataframe_transpose$sampleID, "_"), function(x) x[1])
  rownames(dataframe_transpose) <- NULL
  dataframe_transpose <- dataframe_transpose[,!colnames(dataframe_transpose) == "sampleID"] 
  #combine the row with same sample ID (there are two lanes for each sample ID) 
  dataframe_extract <- dataframe_transpose %>%
    group_by(Sample.ID_Segal.Lab) %>%
    summarize_all(.funs = sum)
  #merge in the key variables from the metafile
  dataframe_extract_key<-merge(dataframe_extract, meta_dataframe_keyvariable, by=meta_dataframe_subjectID, all=TRUE)
  dataframe_extract_key <- dataframe_extract_key %>% select(all_of(colnames_meta_dataframe_keyvariable), everything()) #move the key variables to the front 
  assign(output_dataframe,dataframe_extract_key, envir=.GlobalEnv )
}

MTMG_transform2_data <- function (dataframe_edit, 
                                  sample_type_var,
                                  sample_type_keep=list(),
                                  meta_dataframe_var_remove=list(),
                                  output_dataframe){
  #keeping only tongue samples 
  dataframe_edit2<-dataframe_edit %>% filter(dataframe_edit[[sample_type_var]] %in% sample_type_keep)
  dataframe_edit2_T<-data.frame(t(dataframe_edit2))
  dataframe_edit2_T$sample<-rownames(dataframe_edit2_T)
  dataframe_edit2_T<-dataframe_edit2_T %>% select(sample, everything())
  dataframe_edit2_T <- dataframe_edit2_T[!(rownames(dataframe_edit2_T) %in% meta_dataframe_var_remove), ]
  assign(output_dataframe,dataframe_edit2_T, envir=.GlobalEnv )
}

###############Metabolites based on genome of MOC###############
MOC_metabolites_renamed<-read.csv("MOC_metabolites_metascyc.csv")      
      
#####pathways that are related to metabolism on KEGG#####
microbial_metabolism_kegg<-read.csv("microbial_metabolism_kegg.csv",header = FALSE)
Mouse_map <- read.delim2("Mouse_Map.txt", sep="\t")
###############Ex vivo###############
####exvivo: MG pathway output####
MG_count_ko_pathway_final_exvivo<-read.csv("MG_count_ko_pathway_final_exvivo.csv")
beta_diversity_nonPS_KW(input_table=MG_count_ko_pathway_final_exvivo,
                        sample_type_var_name="Innoculation", 
                        sample_types=c("Prevotella", "Smitis", "Veionella"), 
                        sample_type_color=c("#E4B9A1", "#D4F1AC", "#C9A9DC"), 
                        p_value="yes", 
                        turn_relative_ab="yes",
                        p_value_location="TR",
                        output_name="SFigure5C__PCA_COPD_functional_MG_pathway_exvivo")

####exvivo: MT pathway output####
MT_count_ko_pathway_final_exvivo<-read.csv("MT_count_ko_pathway_final_exvivo.csv")
beta_diversity_nonPS_KW(input_table=MT_count_ko_pathway_final_exvivo,
                        sample_type_var_name="Innoculation", 
                        sample_types=c("Prevotella", "Smitis", "Veionella"), 
                        sample_type_color=c("#E4B9A1", "#D4F1AC", "#C9A9DC"), 
                        p_value="yes", 
                        turn_relative_ab="yes",
                        output_name="SFigure5C__PCA_COPD_functional_MT_pathway_exvivo")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

####Lung ONLY####
#####LUNG: MG ortho output#####
MG_count_ko_ortho_final_lung<-read.csv("MG_count_ko_ortho_final_lung.csv")
Edger_nonPS_KW(input_table=MG_count_ko_ortho_final_lung, 
               sample_type_var_name="Innoculation", 
               sample_types=c("PBS_Lung","MOC_Lung"), 
               sample_type_color=c( "ivory3", "#86A7C7"),  
               FDR_cut_off=0.2,
               number_display=20,
               plot_title="EdgeR: Metagenome ortho- Lung",
               output_name="EdgeR_COPD_functional_MG_ortho_Lung_MOC__PBS",
               legend_onplot="yes",
               width=13)

#####LUNG: MG module output#####
MG_count_ko_module_final_lung<-read.csv("MG_count_ko_module_final_lung.csv")
Edger_nonPS_KW(input_table=MG_count_ko_module_final_lung, 
               sample_type_var_name="Innoculation", 
               sample_types=c("PBS_Lung","MOC_Lung"), 
               sample_type_color=c( "ivory3", "#86A7C7"),  
               FDR_cut_off=0.2,
               number_display=20,
               plot_title="EdgeR: Metagenome module- Lung",
               output_name="EdgeR_COPD_functional_MG_module_Lung_MOC__PBS",
               legend_onplot="yes",
               width=13)

#####LUNG: MG pathway output#####
MG_count_ko_pathway_final_lung<-read.csv("MG_count_ko_pathway_final_lung.csv")
####Figure 4D: MG pathway EdgeR####
#filtering to include only KEGG pathway on metabolism 
MG_count_ko_pathway_final_lung_FILTERED <- MG_count_ko_pathway_final_lung %>%
  filter(sample == "Innoculation" | 
           sub("\\..*", "", sample) %in% microbial_metabolism_kegg$V1)

Edger_nonPS_KW(input_table=MG_count_ko_pathway_final_lung_FILTERED,
               sample_type_var_name="Innoculation", 
               sample_types=c("PBS_Lung","MOC_Lung"), 
               sample_type_color=c( "ivory3", "#86A7C7"),  
               FDR_cut_off=0.2,
               number_display=20,
               graph_option="lollipop",
               plot_title="EdgeR: Metagenome pathway- Lung",
               output_name="Figure_4D___EdgeR_COPD_functional_MG_pathway_Lung_MOC__PBS_all",
               legend_onplot="yes",
               width=12)

              res <- subset(res_export_removeresults, Name != "map00981..Insect.hormone.biosynthesis")
              res <- subset(res, Name != "map00140..Steroid.hormone.biosynthesis")
              p<-ggplot(res, aes(y=reorder(Name,-start), x=logFC,fill=col,size=abundance)) +
                geom_point(color="black",alpha=0.8,shape=21)+
                geom_segment(data=res[res$adj.P.Val<0.2,],aes(yend=reorder(Name,-start)), xend=(-30), color= "black", linetype = "solid",linewidth=1)+ 
                scale_fill_manual(values=c("B"="#86A7C7","A"="ivory3","D"="white"))+ 
                scale_size_continuous(name="Relative Abundance",range=c(5, 20))+                
                ggtitle("EdgeR: Metagenome pathway- Lung")+
                theme(panel.background = element_blank(),
                      panel.border=element_rect(fill=NA),
                      panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
                      panel.grid.minor = element_blank(),
                      strip.background=element_blank(),
                      axis.title=element_text(size=20,face="bold"),
                      axis.text.x=element_text(size=18, face="bold"),
                      axis.text.y=element_text(face="bold",size=10),
                      axis.ticks=element_line(colour="gray70"),
                      legend.background = element_rect(color=NA),
                      legend.key = element_rect(colour = "transparent", fill = "white"),
                      plot.title=element_text(size=23, face="bold"))+
                xlab("") +
                ylab("")+
                geom_vline(xintercept=0, color="red",linetype="dashed")+
                guides(fill="none")
              pdf_output<-paste0("Figure_4D___EdgeR_COPD_functional_MG_pathway_Lung_MOC__PBS_FINAL",".pdf")
              pdf(pdf_output, width=12, height=8)
              show(p)
              dev.off()
                                    
#####LUNG: MT ortho output#####
MT_count_ko_ortho_final_lung<-read.csv("MT_count_ko_ortho_final_lung.csv")              
Edger_nonPS_KW(input_table=MT_count_ko_ortho_final_lung, 
               sample_type_var_name="Innoculation", 
               sample_types=c("PBS_Lung","MOC_Lung"), 
               sample_type_color=c( "ivory3", "#86A7C7"),  
               FDR_cut_off=0.2,
               number_display=20,
               plot_title="EdgeR: Metatranscriptome ortho- Lung",
               output_name="EdgeR_COPD_functional_MT_ortho_Lung_MOC__PBS",
               legend_onplot="yes",
               width=15.5)

#####LUNG: MT module output#####
MT_count_ko_module_final_lung<-read.csv("MT_count_ko_module_final_lung.csv")
Edger_nonPS_KW(input_table=MT_count_ko_module_final_lung, 
               sample_type_var_name="Innoculation", 
               sample_types=c("PBS_Lung","MOC_Lung"), 
               sample_type_color=c( "ivory3", "#86A7C7"),  
               FDR_cut_off=0.2,
               number_display=20,
               plot_title="EdgeR: Metatranscriptome module- Lung",
               output_name="EdgeR_COPD_functional_MT_module_Lung_MOC__PBS",
               legend_onplot="yes",
               width=12)

#####LUNG: MT pathway output#####
MT_count_ko_pathway_final_lung<-read.csv("MT_count_ko_pathway_final_lung.csv")
####Figure 4E: MT pathway EdgeR####
#filtering to include only KEGG pathway on metabolism 
MT_count_ko_pathway_final_lung_FILTERED <- MT_count_ko_pathway_final_lung %>%
  filter(sample == "Innoculation" | 
           sub("\\..*", "", sample) %in% microbial_metabolism_kegg$V1)

Edger_nonPS_KW(input_table=MT_count_ko_pathway_final_lung_FILTERED,
               sample_type_var_name="Innoculation", 
               sample_types=c("PBS_Lung","MOC_Lung"), 
               sample_type_color=c( "ivory3", "#86A7C7"),  
               FDR_cut_off=0.2,
               number_display=20,
               graph_option="lollipop",
               plot_title="EdgeR: Metatranscriptome pathway- Lung",
               output_name="Figure_4E___EdgeR_COPD_functional_MT_pathway_Lung_MOC__PBS_all",
               legend_onplot="yes",
               width=12)

                res <- subset(res_export_removeresults, Name != "map00981..Insect.hormone.biosynthesis")
                res <- subset(res, Name != "map00140..Steroid.hormone.biosynthesis")
                res <- subset(res, Name != "map00100..Steroid.biosynthesis")
                p<-ggplot(res, aes(y=reorder(Name,-start), x=logFC,fill=col,size=abundance)) +
                  geom_point(color="black",alpha=0.8,shape=21)+
                  geom_segment(data=res[res$adj.P.Val<0.2,],aes(yend=reorder(Name,-start)), xend=(-30), color= "black", linetype = "solid",linewidth=1)+ 
                  scale_fill_manual(values=c("B"="#86A7C7","A"="ivory3","D"="white"))+ 
                  scale_size_continuous(name="Relative Abundance",range=c(5, 20))+                
                  ggtitle("EdgeR: Metatranscriptome pathway- Lung")+
                  theme(panel.background = element_blank(),
                        panel.border=element_rect(fill=NA),
                        panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
                        panel.grid.minor = element_blank(),
                        strip.background=element_blank(),
                        axis.title=element_text(size=20,face="bold"),
                        axis.text.x=element_text(size=18, face="bold"),
                        axis.text.y=element_text(face="bold",size=10),
                        axis.ticks=element_line(colour="gray70"),
                        legend.background = element_rect(color=NA),
                        legend.key = element_rect(colour = "transparent", fill = "white"),
                        plot.title=element_text(size=23, face="bold"))+
                  xlab("") +
                  ylab("")+
                  geom_vline(xintercept=0, color="red",linetype="dashed")+
                  guides(fill="none")
                pdf_output<-paste0("Figure_4E___EdgeR_COPD_functional_MT_pathway_Lung_MOC__PBS_FINAL",".pdf")
                pdf(pdf_output, width=12, height=8)
                show(p)
                dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
metabolites_lung_final_filtered<-read.csv("metabolites_lung_final_filtered.csv")
####PLSDA: MOC Lung vs PBS Lung####
PLS_DA_nonPS_KW(input_table=metabolites_lung_final_filtered, 
                sample_type_var_name="Innoculation", 
                sample_types=c("PBS_Lung","MOC_Lung"), 
                sample_type_color=c( "ivory3", "#86A7C7"),   
                top_display_loading=100, 
                top_display_way=c("top_each_side"),
                flip_loading_direction="yes",
                output_name="PLSDA_COPD_functional_metabolite_MOClung_PBSlung",
                output_format="pdf",
                width=10,
                height=12)

#determine the ALL metabolites known to MOC
PLSDA__loading_significant<-data.frame(PLSDA__loading_significant)

PLSDA_MOC_all<- PLSDA__loading_significant %>% select(c("comp1"))
metabolites_PLSDA_MOC_all<- rownames(PLSDA_MOC_all)
metabolites_PLSDA_MOC_all<-data.frame(metabolites_PLSDA_MOC_all)
colnames(metabolites_PLSDA_MOC_all)<-"Metabolite_name"
metabolites_PLSDA_MOC_all$PLSDA_MOCLUNG<-1
colnames(metabolites_PLSDA_MOC_all)<-c("Metabolite_name", "MOC_LUNG")

metabolites_PLSDA_MOC_all_graph<- PLSDA_MOC_all
metabolites_PLSDA_MOC_all_graph$Metabolite_name<-rownames(PLSDA_MOC_all)
metabolites_PLSDA_MOC_all_matchgraph <- merge(MOC_metabolites_renamed,metabolites_PLSDA_MOC_all_graph, by.x="Skeleton_Metabolite", by.y="Metabolite_name", all.x=T)
metabolites_PLSDA_MOC_all_matchgraph <- metabolites_PLSDA_MOC_all_matchgraph[!is.na(metabolites_PLSDA_MOC_all_matchgraph$comp1),]
metabolites_PLSDA_MOC_all_matchgraph <- metabolites_PLSDA_MOC_all_matchgraph[order(metabolites_PLSDA_MOC_all_matchgraph$comp1), ]

write.csv(metabolites_PLSDA_MOC_all_matchgraph,"metabolites_PLSDA_MOC_all_matchgraph.csv")
metabolites_PLSDA_MOC_all_matchgraph <- metabolites_PLSDA_MOC_all_matchgraph %>% select(c("Skeleton_Metabolite","comp1"))

#flipping direction 
metabolites_PLSDA_MOC_all_matchgraph$comp1<-metabolites_PLSDA_MOC_all_matchgraph$comp1*-1

# Create a horizontal bar plot
metabolites_PLSDA_MOC_all_matchgraph<- metabolites_PLSDA_MOC_all_matchgraph %>%   mutate(color_col = ifelse(comp1 > 0, "MOC_Lung", "PBS_Lung"))
write.csv(metabolites_PLSDA_MOC_all_matchgraph, "PLSDA_MOClung_for_known_MOC_metabolites_all.csv")

#remove all the metabolites with 0s
metabolites_PLSDA_MOC_all_matchgraph<-subset(metabolites_PLSDA_MOC_all_matchgraph,comp1!=0 )

metabolites_PLSDA_MOC_all_matchgraph$Skeleton_Metabolite<-gsub("^l-|^d-|^dl-", "", metabolites_PLSDA_MOC_all_matchgraph$Skeleton_Metabolite) #remove "l-", "dl-", "d-" from the label

####Figure 5A: PLSDA####
pp<-ggplot(metabolites_PLSDA_MOC_all_matchgraph, aes(x = comp1, y = reorder(Skeleton_Metabolite, -desc(comp1)))) +
  geom_bar(stat = "identity", aes(fill=color_col), ) +
  scale_fill_manual("sample type",values=c("#86A7C7", "ivory3")) +
  labs(title = "Horizontal Bar Graph Sorted by comp1",
       x = "comp1",
       y = "Skeleton_Metabolite")+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.minor = element_blank(),strip.background=element_blank(),
        plot.title=element_text(face="bold",hjust=0.5, size = 18), 
        plot.subtitle = element_text(hjust=0.5),
        axis.title=element_text(face="bold", size = 16),
        axis.text.x=element_text(face="bold",size = 12),
        axis.text.y=element_text(face="bold",size = 11.5),
        axis.ticks=element_blank(),
        strip.text = element_text(size = 18),
        plot.margin=unit(c(1,1,1,2),"line"),
        legend.position = "bottom",  # Position the legend at the bottom
        legend.direction = "horizontal") + # Display the legend horizontally
  ylab("Metabolites") +  # Change the x-axis label
  xlab("Component 1 of PLSDA") +   # Change the y-axis label
  ggtitle("PLSDA-(known MOC metabolites)")

pdf(file="Figure_5A___PLSDA_MOClung_for_known_MOC_metabolites_all_FILTERED.pdf", width=14, height=14)
    show(pp)
dev.off()    

#~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~#
####relative abundance for bacteria:####
#Prevotella melaninogenica
#Streptococcus mitis
#Veillonella parvula
#~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~#
#####relab- metagenome MG##### 
      MG_count_taxa_bracken <-read.delim2("MG_counts.taxa.bracken.tsv", quote = "", row.names = NULL, stringsAsFactors = FALSE)
      MG_count_taxa_bracken <- MG_count_taxa_bracken %>% dplyr::rename(taxa=row.names)
      MG_count_taxa_bracken<-MG_count_taxa_bracken[MG_count_taxa_bracken$taxa!="Homo sapiens",] #remove the homo sapiens which is clearly a contaminant
      MG_numeric_columns <- MG_count_taxa_bracken[, sapply(MG_count_taxa_bracken, is.numeric)]
      MG_count_taxa_bracken_Relab<-MG_count_taxa_bracken
      MG_count_taxa_bracken_Relab[,-1]<- lapply(MG_numeric_columns, function(x) x / sum(x))
      MTMG_transform_data(dataframe=MG_count_taxa_bracken_Relab, 
                          meta_dataframe=Mouse_map,
                          meta_dataframe_var_interest=c("Sample.ID_Segal.Lab", "Sample_Type_IS", "Innoculation"),
                          meta_dataframe_subjectID="Sample.ID_Segal.Lab",
                          output_dataframe="MG_count_taxa_bracken_Relab_final", 
                          first_column_name="taxa")
      MG_count_taxa_bracken_selectpathogen<- MG_count_taxa_bracken_Relab_final %>% select(c("Sample.ID_Segal.Lab","Sample_Type_IS","Innoculation",
                                                                                         "Prevotella.melaninogenica", "Streptococcus.mitis", "Veillonella.parvula"))
      MG_count_taxa_bracken_selectpathogen<- MG_count_taxa_bracken_selectpathogen[MG_count_taxa_bracken_selectpathogen$Sample_Type_IS!="n.a.",]
      MG_count_taxa_bracken_selectpathogen$sample_name <- paste(MG_count_taxa_bracken_selectpathogen$Sample_Type_IS, MG_count_taxa_bracken_selectpathogen$Sample_type, sep = "__")
      
      MG_count_taxa_bracken_graph_t<- MG_count_taxa_bracken_selectpathogen %>% select(c("Sample_Type_IS","Innoculation","Prevotella.melaninogenica", "Streptococcus.mitis", "Veillonella.parvula"))
      MG_count_taxa_bracken_graph <- MG_count_taxa_bracken_graph_t %>%
        pivot_longer(cols = -c("Sample_Type_IS","Innoculation"), names_to = "Bacteria", values_to = "Relative_Abundance")
      
      #making a boxplot graph for the three bacteria of interest
      bacteria_color_mapping <- c("Prevotella.melaninogenica" = "#E4B9A1", "Streptococcus.mitis" = "#D4F1AC", "Veillonella.parvula" = "#C9A9DC")
      desired_sample_order <- c("PBS_Lung","MOC_Oral","MOC_Lung")
      
      ####Figure 4B: MG by sample types####
      #setting order for the facet
      facet_order <- c('Tongue', 'Lung')  # The desired order of the x-axis factor
      MG_count_taxa_bracken_graph$Sample_Type_IS <- factor(MG_count_taxa_bracken_graph$Sample_Type_IS, levels = facet_order)
      
      q<-ggplot(MG_count_taxa_bracken_graph, aes(x=factor(Innoculation, levels = desired_sample_order), y=Relative_Abundance, fill=Bacteria)) + 
        geom_boxplot()+
        facet_grid(. ~ Sample_Type_IS)+
        ylab("Relative Abundnace") +  # Change the x-axis label
        xlab("Sample type") +  # Change the y-axis label
        labs(title = "Metagenome") + 
        scale_y_continuous(trans = "log10",  breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1), 
                           labels = c(expression(10^-6), expression(10^-5), expression(10^-4), expression(10^-3), expression(10^-2), expression(10^-1), expression(10^0)),
                           limits = c(1e-6, 1e-2)) +
        scale_fill_manual(values = bacteria_color_mapping)+
        theme(panel.background = element_blank(),
              panel.border=element_rect(fill=NA),
              panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
              panel.grid.minor = element_blank(),strip.background=element_blank(),
              plot.title=element_text(face="bold",hjust=0.5, size = 28), 
              plot.subtitle = element_text(hjust=0.5),
              axis.title=element_text(face="bold", size = 14),
              axis.text.x=element_text(face="bold",size = 12),
              axis.text.y=element_text(face="bold",size = 13),
              axis.ticks=element_blank(),
              strip.text = element_text(size = 18),
              plot.margin=unit(c(1,1,1,1),"line"),
              legend.position = "bottom",  # Position the legend at the bottom
              legend.direction = "horizontal",
              legend.key = element_rect(fill = "transparent", colour = "transparent"))  # Display the legend horizontally
      pdf_output="Figure_4B___MG_bacteria_relab_by_sample_types.pdf"
      pdf(file=pdf_output, width=9, height=6)
              show(q)
      dev.off()        
      
      #create a table comparing each bacteria between all the pair to see if significant or not. will use wilcox test given nonparametric
      
      chosen_bacteria<- c("Prevotella.melaninogenica", "Streptococcus.mitis", "Veillonella.parvula")
      Sample_types_compare <- c("PBS_Lung__and__MOC_Oral","PBS_Lung__and__MOC_Lung","MOC_Oral__and__MOC_Lung")
      Sample_types_IS <- c("Lung","Tongue")
      
      #create an empty dataframe which will be the frame for the statistics that will be merge into 
      MG_statistics <- data.frame(
        .y. = character(0),
        Sample_types_IS=character(0),
        bacteria=character(0),
        group1 = character(0),
        group2 = character(0),
        group1_mean=integer(0),
        group2_mean=integer(0),
        n1 = integer(0),
        n2 = integer(0),
        statistic = integer(0),
        p = integer(0),
        p.signif = character(0))
      
      for (i in Sample_types_IS) {
        for (j in chosen_bacteria){
          for (k in Sample_types_compare){
            subsample_select <- MG_count_taxa_bracken_graph[MG_count_taxa_bracken_graph$Sample_Type_IS==i,]      
            subsample_select2 <- subsample_select[subsample_select$Bacteria==j, ]
            
            compare_sublist <- unlist(strsplit(k, "__and__"))
            subsample_select3 <- subsample_select2 %>% filter(Innoculation %in% compare_sublist)
            
            stat.test <- subsample_select3 %>% 
              wilcox_test(Relative_Abundance ~ Innoculation, ref.group=compare_sublist[1]) %>%
              add_significance()
            stat.test<-data.frame(stat.test)
            
            mean_values <- subsample_select3 %>%
              dplyr::group_by(Innoculation) %>%
              dplyr::summarize(mean_Relative_Abundance = mean(Relative_Abundance))
            
            stat.test$group1_mean <- mean_values$mean_Relative_Abundance[mean_values$Innoculation==compare_sublist[1]]
            stat.test$group2_mean <- mean_values$mean_Relative_Abundance[mean_values$Innoculation==compare_sublist[2]]
            stat.test$Sample_type_IS<-i
            stat.test$bacteria<-j
            
            MG_statistics<- rbind(MG_statistics, stat.test)
          }
        }
      }
      MG_statistics<-MG_statistics %>% select(-c(".y.")) #remove the column .y. which is not needed
      MG_statistics <- MG_statistics %>% arrange(Sample_type_IS, bacteria, group1, group2) #sorting based on the sample types
      
      #reorder the export so column order is set
      MG_statistics<- MG_statistics %>% select(c("Sample_type_IS","bacteria","group1","n1","group1_mean",  
                                                 "group2","n2", "group2_mean","statistic","p","p.signif" ))
      write.csv(MG_statistics, "MG_bacteria_relab_by_sample_types.csv",row.names=FALSE)
#####relab- metatranscriptome MT##### 
      MT_count_taxa_bracken <-read.delim2("MT_counts.taxa.bracken.txt", sep="\t")
      MT_count_taxa_bracken<-MT_count_taxa_bracken[MT_count_taxa_bracken$taxa!="Homo sapiens",] #remove the homo sapiens which is clearly a contaminant
      
      MT_numeric_columns <- MT_count_taxa_bracken[, sapply(MT_count_taxa_bracken, is.numeric)]
      MT_count_taxa_bracken_Relab<-MT_count_taxa_bracken
      MT_count_taxa_bracken_Relab[,-1]<- lapply(MT_numeric_columns, function(x) x / sum(x))
      MTMG_transform_data(dataframe=MT_count_taxa_bracken_Relab, 
                          meta_dataframe=Mouse_map,
                          meta_dataframe_var_interest=c("Sample.ID_Segal.Lab", "Sample_Type_IS", "Innoculation"),
                          meta_dataframe_subjectID="Sample.ID_Segal.Lab",
                          output_dataframe="MT_count_taxa_bracken_Relab_final", 
                          first_column_name="taxa")
      MT_count_taxa_bracken_selectpathogen<- MT_count_taxa_bracken_Relab_final %>% select(c("Sample.ID_Segal.Lab","Sample_Type_IS","Innoculation",
                                                                                            "Prevotella.melaninogenica", "Streptococcus.mitis", "Veillonella.parvula"))
      MT_count_taxa_bracken_selectpathogen<- MT_count_taxa_bracken_selectpathogen[MT_count_taxa_bracken_selectpathogen$Sample_Type_IS!="n.a.",]
      MT_count_taxa_bracken_selectpathogen$sample_name <- paste(MT_count_taxa_bracken_selectpathogen$Sample_Type_IS, MT_count_taxa_bracken_selectpathogen$Sample_type, sep = "__")
      
      MT_count_taxa_bracken_graph_t<- MT_count_taxa_bracken_selectpathogen %>% select(c("Sample_Type_IS","Innoculation","Prevotella.melaninogenica", "Streptococcus.mitis", "Veillonella.parvula"))
      MT_count_taxa_bracken_graph <- MT_count_taxa_bracken_graph_t %>%
        pivot_longer(cols = -c("Sample_Type_IS","Innoculation"), names_to = "Bacteria", values_to = "Relative_Abundance")
      
      #making a boxplot graph for the three bacteria of interest
      bacteria_color_mapping <- c("Prevotella.melaninogenica" = "#E4B9A1", "Streptococcus.mitis" = "#D4F1AC", "Veillonella.parvula" = "#C9A9DC")
      desired_sample_order <- c("PBS_Lung","MOC_Oral","MOC_Lung")
      
      ####Figure 4C: MT by sample types####
      #setting order for the facet
      facet_order <- c('Tongue', 'Lung')  # The desired order of the x-axis factor
      MT_count_taxa_bracken_graph$Sample_Type_IS <- factor(MT_count_taxa_bracken_graph$Sample_Type_IS, levels = facet_order)
      
      q<-ggplot(MT_count_taxa_bracken_graph, aes(x=factor(Innoculation, levels = desired_sample_order), y=Relative_Abundance, fill=Bacteria)) + 
        geom_boxplot()+
        facet_grid(. ~ Sample_Type_IS)+
        ylab("Relative Abundnace") +  # Change the x-axis label
        xlab("Sample type") +  # Change the y-axis label
        labs(title = "Metatranscriptome") + 
        scale_y_continuous(trans = "log10",  breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1), 
                           labels = c(expression(10^-6), expression(10^-5), expression(10^-4), expression(10^-3), expression(10^-2), expression(10^-1), expression(10^0)),
                           limits = c(1e-6, 1e-2)) +
        scale_fill_manual(values = bacteria_color_mapping)+
        theme(panel.background = element_blank(),
              panel.border=element_rect(fill=NA),
              panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),
              panel.grid.minor = element_blank(),strip.background=element_blank(),
              plot.title=element_text(face="bold",hjust=0.5, size = 28), 
              plot.subtitle = element_text(hjust=0.5),
              axis.title=element_text(face="bold", size = 14),
              axis.text.x=element_text(face="bold",size = 12),
              axis.text.y=element_text(face="bold",size = 13),
              axis.ticks=element_blank(),
              strip.text = element_text(size = 18),
              plot.margin=unit(c(1,1,1,1),"line"),
              legend.position = "bottom",  # Position the legend at the bottom
              legend.direction = "horizontal",
              legend.key = element_rect(fill = "transparent", colour = "transparent"))  # Display the legend horizontally
      pdf_output="Figure_4C___MT_bacteria_relab_by_sample_types.pdf"
      pdf(file=pdf_output, width=9, height=6)
            show(q)
      dev.off()        
      
      #create an empty dataframe which will be the frame for the statistics that will be merge into 
      MT_statistics <- data.frame(
        .y. = character(0),
        Sample_types_IS=character(0),
        bacteria=character(0),
        group1 = character(0),
        group2 = character(0),
        group1_mean=integer(0),
        group2_mean=integer(0),
        n1 = integer(0),
        n2 = integer(0),
        statistic = integer(0),
        p = integer(0),
        p.signif = character(0))
      
      for (i in Sample_types_IS) {
        for (j in chosen_bacteria){
          for (k in Sample_types_compare){
            subsample_select <- MT_count_taxa_bracken_graph[MT_count_taxa_bracken_graph$Sample_Type_IS==i,]      
            subsample_select2 <- subsample_select[subsample_select$Bacteria==j, ]
            
            compare_sublist <- unlist(strsplit(k, "__and__"))
            subsample_select3 <- subsample_select2 %>% filter(Innoculation %in% compare_sublist)
            
            stat.test <- subsample_select3 %>% 
              wilcox_test(Relative_Abundance ~ Innoculation, ref.group=compare_sublist[1]) %>%
              add_significance()
            stat.test<-data.frame(stat.test)
            
            mean_values <- subsample_select3 %>%
              dplyr::group_by(Innoculation) %>%
              dplyr::summarize(mean_Relative_Abundance = mean(Relative_Abundance))
            
            stat.test$group1_mean <- mean_values$mean_Relative_Abundance[mean_values$Innoculation==compare_sublist[1]]
            stat.test$group2_mean <- mean_values$mean_Relative_Abundance[mean_values$Innoculation==compare_sublist[2]]
            stat.test$Sample_type_IS<-i
            stat.test$bacteria<-j
            
            MT_statistics<- rbind(MT_statistics, stat.test)
          }
        }
      }
      MT_statistics<-MT_statistics %>% select(-c(".y.")) #remove the column .y. which is not needed
      MT_statistics <- MT_statistics %>% arrange(Sample_type_IS, bacteria, group1, group2) #sorting based on the sample types
      #reorder the export so column order is set
      MT_statistics<- MT_statistics %>% select(c("Sample_type_IS","bacteria","group1","n1","group1_mean",  
                                                 "group2","n2", "group2_mean","statistic","p","p.signif" ))
      
      write.csv(MT_statistics, "MT_bacteria_relab_by_sample_types.csv",row.names=FALSE)

##### invivo isotope percent#####
isotope_meta_histogram <- function (input_file,
                                    metabolite,
                                    compare_group1=list(),
                                    compare_group2=list(),
                                    group1_name, group2_name,
                                    sample_types=list(), 
                                    sample_type_color=list(), 
                                    p_value=c("yes","no"), 
                                    output_name, 
                                    width, height, titlelabel=NA, 
                                    xlabel_size=14, ylabel_size=14, axis_title_size=16,
                                    xlabel, ylabel,     
                                    output_format=c("pdf", "png")) {
  
  if(missing(output_format)) { output_format<-"pdf"} #default being output_format is pdf
  if (output_format=="pdf"){
    if (missing(width)) width<-7
    if (missing(height)) height<-10
  } else if (output_format=="png"){
    if (missing(width)) width<-100
    if (missing(height)) height<-150
  } 
  #making sure only yes and no is selected for p_value
  p_value <- match.arg(p_value)
  if(missing(p_value)){p_value<-"yes"}
  
  data_percent<- input_file
  data_percent_mets<-data_percent[data_percent$Metabolite==metabolite,]
  data_percent_mets<-data_percent_mets %>% select(-c("Metabolite"))

  data_percent_mets_group1<-t(subset(data_percent_mets, select = compare_group1))
  data_percent_mets_group1<-data.frame(data_percent_mets_group1)
  colnames(data_percent_mets_group1)<-"percent"
  data_percent_mets_group1$group <- group1_name
  data_percent_mets_group2<-t(subset(data_percent_mets, select = compare_group2))
  data_percent_mets_group2<-data.frame(data_percent_mets_group2)
  colnames(data_percent_mets_group2)<-"percent"
  data_percent_mets_group2$group <- group2_name
  
  stacked_df <- rbind(data_percent_mets_group1, data_percent_mets_group2)
  sample_type_color_t<-sample_type_color
  names(sample_type_color_t) <-sample_types
  stacked_df$group<-factor(stacked_df$group, levels = sample_types)
  
  
  stacked_df<<-stacked_df
  
  #Setting up the graphs
  p<-ggplot(stacked_df, aes(group, percent, fill=group)) +
    stat_boxplot(geom = "errorbar",width = 0.3, color = "gray70")+
    geom_boxplot() +
    scale_fill_manual(values=sample_type_color_t) +
    scale_x_discrete(labels=sample_types) +
    scale_y_continuous() + 
    theme_classic() +
    labs(title=titlelabel,
         x = xlabel, 
         y= ylabel) +
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          axis.title=element_text(face="bold", size=axis_title_size),
          axis.text.x =element_text(size=xlabel_size),
          axis.text.y =element_text(size=ylabel_size),
          legend.position="none")
  p <- p + geom_jitter(width=0.3, size=1, shape= 1, color="azure4") #if no paired line, then will have jitter 
  if (p_value=="yes") { ##if want p_value included in the diagram
    stat.test <- stacked_df %>% wilcox_test(percent ~ group, alternative="less") %>% add_significance()
    stat.test<- stat.test %>% add_xy_position(x = "group", y.trans = function(x){x*1.1})
    if (!exists("step.increase")) {step.increase<-0.05} 
    #limiting p value to 3 digits after decimal
    stat.test$p<-pvalr(stat.test$p, digits = 3)
    p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.02, step.increase = step.increase,  inherit.aes = F, size=5)
  }
  ###determine sample_types string length
  string_length <- max(nchar(as.character(stacked_df$group)))
  if (string_length>20) { #if the sample type string is too long then will title the x axis text 
    p<-p+ theme(axis.text.x = element_text(angle = 8, hjust = 0.5, vjust=0.5))
  }
  if (output_format=="pdf"){
    pdf_output<-paste0(output_name,".pdf")
    pdf(file=pdf_output, width=width, height=height)
    show(p)
    dev.off()
  } else if (output_format=="png"){
    png_output<-paste0(output_name,".png")
    png(file=png_output, res = 300, width=width, height=height , units='mm')
    show(p)
    dev.off()
  } 
  csv_output<-paste0(output_name,".csv")
  write.csv(percent, csv_output)
}

isotope_meta_histogram_multiple <- function(input_file,
                                            metabolites = list(),
                                            metabolite_final_name = NULL,
                                            compare_group1 = list(),
                                            compare_group2 = list(),
                                            group1_name, group2_name,
                                            group1_color, group2_color, 
                                            sample_types = list(), 
                                            output_name, 
                                            width = NULL, height = NULL, titlelabel = NA, 
                                            xlabel_size = 14, ylabel_size = 14, axis_title_size = 16,
                                            xlabel, ylabel,     
                                            output_format = c("pdf", "png"),
                                            show_pvalue = FALSE, significance_level = 0.05,
                                            log_transform = FALSE,
                                            angle_labels = TRUE) {
  
  # Default argument values
  output_format <- match.arg(output_format, choices = c("pdf", "png"), several.ok = FALSE)
  if (is.null(width))  width <- ifelse(output_format == "pdf", 7, 10)
  if (is.null(height)) height <- ifelse(output_format == "pdf", 10, 7)
  
  # Create mapping between metabolites and metabolite_final_name
  if (!is.null(metabolite_final_name)) {
    names(metabolite_final_name) <- metabolites
  } else {
    metabolite_final_name <- metabolites  # Default to original metabolite names
  }
  
  # Setting up the colors
  fill_colors_t <- c(group1_color, group2_color)
  names(fill_colors_t) <- c(group1_name, group2_name)
  all_color <- fill_colors_t
  
  
  fill_colors_t <- c("white", group2_color)
  names(fill_colors_t) <- c(group1_name, group2_name)
  all_fill_color <- fill_colors_t
  
  # Prepare combined data frame for all metabolites
  combined_df <- data.frame()
  p_values <- data.frame()
  
  for (metabolite in metabolites) {
    data_percent_mets <- input_file %>% filter(Metabolite == metabolite) %>% select(-Metabolite)
    
    # Prepare data for group 1
    data_percent_mets_group1 <- data_percent_mets %>% select(all_of(compare_group1)) %>%
      t() %>% as.data.frame() %>% mutate(group = group1_name, metabolite = metabolite)
    colnames(data_percent_mets_group1)[1] <- "percent"
    
    # Prepare data for group 2
    data_percent_mets_group2 <- data_percent_mets %>% select(all_of(compare_group2)) %>%
      t() %>% as.data.frame() %>% mutate(group = group2_name, metabolite = metabolite)
    colnames(data_percent_mets_group2)[1] <- "percent"
    
    # Combine data frames
    stacked_df <- bind_rows(data_percent_mets_group1, data_percent_mets_group2)
    combined_df <- bind_rows(combined_df, stacked_df)
    
    if (show_pvalue) {
      # Perform Wilcoxon test
      p_value <- wilcox.test(data_percent_mets_group1$percent, data_percent_mets_group2$percent, alternative = "less", exact = FALSE)$p.value
      p_values <- rbind(p_values, data.frame(metabolite = metabolite, p_value = p_value))
    }
  }
  
  combined_df$group <- factor(combined_df$group, levels = sample_types)
  combined_df$metabolite <- factor(combined_df$metabolite, levels = metabolites)
  combined_df <<- combined_df
  
  combined_df <- combined_df %>% filter(!is.na(percent))
  if(log_transform==TRUE) {
    combined_df$percent <- log10(combined_df$percent + 1)
  }
  
  # Calculate means and standard errors for each metabolite and group
  summary_df <- combined_df %>%
    group_by(metabolite, group) %>%
    dplyr::summarise(mean_percent = mean(percent), se_percent = sd(percent) / sqrt(length(percent)), .groups = 'drop')
  
  summary_df <- summary_df %>% mutate(fill_color = case_when(group == group2_name ~ group2_color, group == group1_name ~ group1_color))
  
  summary_df$fill_color <- factor(summary_df$fill_color, levels = c(group1_color, group2_color))
  
  summary_df <<- summary_df
  
  if (show_pvalue) {
    summary_df <- summary_df %>%
      dplyr::left_join(p_values, by = "metabolite") %>%
      mutate(significant = ifelse(p_value < significance_level, "*", ""))
  }
  
  # Calculate maximum y value for all bars
  max_y <- max(summary_df$mean_percent + 1.96 * summary_df$se_percent)
  
  # Create plot
  p <- ggplot(summary_df, aes(x = factor(metabolite, levels = metabolites), y = mean_percent, fill = group, color = group)) +
    geom_errorbar(aes(ymin = mean_percent - 1.96 * se_percent, ymax = mean_percent + 1.96 * se_percent),
                  width = 0.5,  # Adjust width if necessary
                  position = position_dodge(width = 0.7),  # Match dodge width of bars
                  color = "black") +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7, size = 1.5) +  # Ensure width is consistent
    scale_fill_manual(values = all_fill_color) +  
    scale_color_manual(values = all_color) +  
    theme_classic() +
    labs(title = titlelabel, x = xlabel, y = ylabel, fill = "Group") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_text(face = "bold", size = axis_title_size),
          axis.text.x = element_text(size = xlabel_size, angle = ifelse(angle_labels, 45, 0), hjust = ifelse(angle_labels, 1, 0.5)),
          axis.text.y = element_text(size = ylabel_size)) +
    scale_x_discrete(labels = metabolite_final_name) +
    guides(fill = "none")
  
  if (show_pvalue) {
    p <- p + geom_text(aes(x = factor(metabolite, levels = metabolites), y = max_y + 0.05 * max_y, label = significant),
                       color = "black", vjust = 0, size = 8)  
  }
  
  # Save plot
  if (output_format == "pdf") {
    pdf(file = paste0(output_name, ".pdf"), width = width, height = height)
    print(p)
    dev.off()
  } else if (output_format == "png") {
    png(file = paste0(output_name, ".png"), res = 300, width = width, height = height, units = 'in')
    print(p)
    dev.off()
  }
  
  # Save data
  write.csv(summary_df, file = paste0(output_name, ".csv"), row.names = FALSE)
}

#~~~~~~~~~~~~~~~~~~~~~~#
#####isotope percent#####
#~~~~~~~~~~~~~~~~~~~~~~#
Metabolites_isotope_percent2 <-read.csv("Metabolites_isotope_percent2_export.csv")
####Figure 6E: invivo isotope relative intensity####
isotope_meta_histogram_multiple(Metabolites_isotope_percent2,
                                metabolites = c("adenosine-13c-5", "inosine-13c-5", "l-glutamate-13c-2_1", "l-methionine-13c-5_1"), 
                                metabolite_final_name = c("adenosine", "inosine", "glutamate", "methionine"),
                                compare_group1 = c("M11_1h_C12_lung", "M12_1h_C12_lung", "M13_1h_C12_lung", "M14_1h_C12_lung"),
                                compare_group2 = c("M1_1h_C13_lung", "M2_1h_C13_lung", "M3_1h_C13_lung", "M4_1h_C13_lung"),
                                group1_name = "lung_12", group2_name = "lung_13",
                                group1_color = "#7EAA73", 
                                group2_color = "#78C1CC",
                                sample_types = c("lung_12", "lung_13"), 
                                xlabel_size = 10, ylabel_size = 20, axis_title_size = 23,
                                width = 10, 
                                xlabel = "Metabolite", ylabel = "Relative Intensity (log 10 of percent)", titlelabel = "C13 Labeling - 1 Hour",
                                output_name = "Figure_6E___all_metabolites_1hr_invivo", show_pvalue = TRUE, log_transform = TRUE,
                                angle_labels = FALSE)

isotope_meta_histogram_multiple(Metabolites_isotope_percent2,
                                metabolites=c("adenosine-13c-5", "inosine-13c-5", "l-glutamate-13c-2_1", "l-methionine-13c-5_1"), 
                                metabolite_final_name=c("adenosine", "inosine", "glutamate",  "methionine"),
                                compare_group1=c("M15_6h_C12_lung","M16_6h_C12_lung","M17_6h_C12_lung","M18_6h_C12_lung","M19_6h_C12_lung","M20_6h_C12_lung"),
                                compare_group2=c("M5_6h_C13_lung","M6_6h_C13_lung","M7_6h_C13_lung","M8_6h_C13_lung","M9_6h_C13_lung"),
                                group1_name="lung_12", group2_name="lung_13",
                                group1_color="#709666",
                                group2_color="#6FB3BD",
                                sample_types=c("lung_12", "lung_13"), 
                                xlabel_size = 10, ylabel_size = 20, axis_title_size = 23,
                                width=10,
                                xlabel="Metabolite", ylabel="Relative Intensity (log 10 of percent)", titlelabel="C13 Labeling - 6 Hour",
                                output_name="Figure_6F___all_metabolites_6hr_invivo", show_pvalue=TRUE,   log_transform = TRUE,
                                angle_labels = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
barplot_KO_metabolites<-function(dataframe,plot_title,output_name){        
  df<-dataframe
  print(df)
  # Calculate the total counts for each column
  df <- mutate(df, total_MOC_Lung = sum(MOC_Lung),
               total_PBS_Lung = sum(PBS_Lung))
  
  # Calculate the percentage of each microbe for MOC_Lung and PBS_Lung
  df <- mutate(df, MOC_Lung = (MOC_Lung / total_MOC_Lung) * 100,
               PBS_Lung = (PBS_Lung / total_PBS_Lung) * 100)
  
  # Filter out rows where both percentages are 0
  df <- filter(df, MOC_Lung != 0 | PBS_Lung != 0)
  
  df$Microbe<-rownames(df)
  # Reorder the Microbe factor based on MOC_Lung_Percentage
  df$Microbe <- factor(df$Microbe, levels = df$Microbe[order(-df$MOC_Lung)])
  
  # Reshape the data into long format for ggplot
  df_long <- tidyr::pivot_longer(df, cols = c(MOC_Lung, PBS_Lung),
                                 names_to = "Condition", values_to = "Percentage")
  
  # Create the bar plot
  vv<-ggplot(df_long, aes(x = Microbe, y = Percentage, fill = Condition)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Taxa (Genus)", y = "Percentage") +
    scale_fill_manual(values = c("#86A7C7", "ivory3"), name = "Condition") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    ggtitle(plot_title)
  output_name_final<-paste0(paste0("barplot_",output_name),".pdf")
  pdf(output_name_final)
  show(vv)
  dev.off()
}


#####Species level MG/MT KO#####
#######MT KO match#######
setwd("~/MT_KO_species") #set the directory to the raw file directory 
metabolite_ko_files <- list.files(getwd())

ii<-1
all_metabolite_percentage_MT <- data.frame(
  row = character(),        # Character column
  category = character(),   # Character column
  value = numeric(),        # Numeric column
  metabolites = character() # Character column
)

for (file in metabolite_ko_files){
  print(ii)
  metabolite_name <- sub("\\.tsv$", "", file)
  print(metabolite_name)
  metabolite_KO <- read.table(file, header = TRUE, sep = "\t")
  metabolite_KO$taxa<-rownames(metabolite_KO)
  metabolite_KO <- metabolite_KO %>% select(c("taxa", everything()))
  if (nrow(metabolite_KO) == 0 ) {
    print(paste0(file," is empty"))
  } else {
    print(file)
    MTMG_transform_data(dataframe=metabolite_KO, 
                        meta_dataframe=Mouse_map,
                        meta_dataframe_var_interest=c("Sample.ID_Segal.Lab", "Sample_Type_IS", "Innoculation"),
                        meta_dataframe_subjectID="Sample.ID_Segal.Lab",
                        output_dataframe="metabolite_KO_updated", 
                        first_column_name="taxa")
    MTMG_transform2_data (dataframe_edit=metabolite_KO_updated, 
                          sample_type_var="Sample_Type_IS",
                          sample_type_keep="Lung",
                          meta_dataframe_var_remove=c("Sample.ID_Segal.Lab", "Sample_Type_IS"),
                          output_dataframe="metabolite_KO_updated2")
    metabolite_KO_updated2[is.na(metabolite_KO_updated2)] <- 0
    metabolite_KO_updated3<-metabolite_KO_updated2[-1,-1]
    metabolite_KO_updated3 <- apply(metabolite_KO_updated3, 2, as.numeric)
    rownames(metabolite_KO_updated3)<-rownames(metabolite_KO_updated2[-1,-1])
    sample_type_KO<- metabolite_KO_updated2[1,- 1]
    # Create a new dataframe with the sums
    metabolite_KO_updated4 <- data.frame(
      MOC_Lung = rowSums(metabolite_KO_updated3[,which(sample_type_KO == "MOC_Lung")]),
      PBS_Lung = rowSums(metabolite_KO_updated3[,which(sample_type_KO == "PBS_Lung")]),
      MOC_Oral = rowSums(metabolite_KO_updated3[,which(sample_type_KO == "MOC_Oral")])
    )
    metabolite_KO_updated4 <- metabolite_KO_updated4 %>% select(c("MOC_Lung","PBS_Lung"))
    
    metabolite_KO_updated4<<-metabolite_KO_updated4
    
    #adding the barplots 
    simplified_metabolite_name <- gsub("\\.species", "", metabolite_name)
    name_output<-paste0("MT KEGG ortho percentage: ",simplified_metabolite_name)
    MT_simplified_metabolite_name<-paste0("MT_",simplified_metabolite_name)
    barplot_KO_metabolites(metabolite_KO_updated4,name_output,MT_simplified_metabolite_name)
    
    if (all(metabolite_KO_updated4 == 0)){ #only doing the graph for ones that are not all zeros
      print(paste0(file," is all zeros"))
    } else {
      # Rows to exclude
      Prevotella <- c("Prevotella.melaninogenica")
      Streptococcus <- c("Streptococcus.mitis")
      Veillonella <- c("Veillonella.parvula")
      exclude_rows<-c("Prevotella.melaninogenica", "Streptococcus.mitis", "Veillonella.parvula")
      # Sum rows of MOC
      sum_rows_Prevotella <- data.frame(t(colSums(metabolite_KO_updated4[(rownames(metabolite_KO_updated4) %in% Prevotella), , drop = FALSE])))
      rownames(sum_rows_Prevotella)<-"Prevotella.melaninogenica"
      sum_rows_Streptococcus <- data.frame(t(colSums(metabolite_KO_updated4[(rownames(metabolite_KO_updated4) %in% Streptococcus), , drop = FALSE])))
      rownames(sum_rows_Streptococcus)<-"Streptococcus.mitis"
      sum_rows_Veillonella <- data.frame(t(colSums(metabolite_KO_updated4[(rownames(metabolite_KO_updated4) %in% Veillonella), , drop = FALSE])))
      rownames(sum_rows_Veillonella)<-"Veillonella.parvula"
      # Sum rows excluding specified ones
      sum_rows_other <- data.frame(t(colSums(metabolite_KO_updated4[!(rownames(metabolite_KO_updated4) %in% exclude_rows), , drop = FALSE])))
      rownames(sum_rows_other)<-"others"
      combined_other<- rbind(metabolite_KO_updated4[(rownames(metabolite_KO_updated4) %in% exclude_rows), , drop = FALSE], sum_rows_other)
      combined_MOC<- rbind(sum_rows_Prevotella, sum_rows_Streptococcus, sum_rows_Veillonella, sum_rows_other)
      combined_MOC_t<-t(combined_MOC)
      
      ########relative abundance
      combined_MOC_t_relab<- as.data.frame(t(apply(combined_MOC_t, 1, function(row) row/sum(row))))
      
      combined_MOC_t_relab$row <- rownames(combined_MOC_t)
      combined_MOC_t_relab<<-combined_MOC_t_relab
      
      df_long <- pivot_longer(combined_MOC_t_relab, cols = -row, names_to = "Category", values_to = "Value")
      print(metabolite_KO_updated4)
      print(combined_MOC_t)
      print(df_long)
      
      df_long$metabolite<- simplified_metabolite_name
      
      all_metabolite_percentage_MT<- rbind(all_metabolite_percentage_MT, df_long)
      
      ii<-ii+1
    } 
  }
}

all_metabolite_percentage_MT <- all_metabolite_percentage_MT %>%
  mutate(Value = ifelse(is.na(Value), 0, Value))

PBS_Lung_isoleucine <- data.frame(
  row = "PBS_Lung",        # Character column
  Category = "NA",   # Character column
  Value = 1,        # Numeric column
  metabolite = "l-isoleucine" # Character column
)

all_metabolite_percentage_MT<-rbind(all_metabolite_percentage_MT, PBS_Lung_isoleucine)

all_metabolite_percentage_MT$Value_updated<-round(all_metabolite_percentage_MT$Value * 100, 2)

custom_colors <- c("Prevotella.melaninogenica" = "#E4B9A1", 
                   "Streptococcus.mitis" = "#D4F1AC", 
                   "Veillonella.parvula" = "#C9A9DC", 
                   "others" = "grey95",
                   "NA"="white")


legend_order <- c("Prevotella.melaninogenica", "Streptococcus.mitis", "Veillonella.parvula", "others", "NA")

#setting PBS_Lung on the left side and then MOC_Lung on the right side 
all_metabolite_percentage_MT$row<-factor(all_metabolite_percentage_MT$row, levels=c("PBS_Lung", "MOC_Lung"))

all_metabolite_percentage_MT<- all_metabolite_percentage_MT[all_metabolite_percentage_MT$metabolite %in% 
                                                              c("l-glutamate","l-tyrosine",  
                                                                "l-methionine", "adenosine", "niacinamide",
                                                                "adenosine_5-monophosphate", "inosine_5-monophosphate"),]

all_metabolite_percentage_MT$metabolite <-gsub("^l-|^d-|^dl-", "", all_metabolite_percentage_MT$metabolite) #remove "l-", "dl-", "d-" from the label

desired_order <- rev(c("niacinamide","methionine", "inosine_5-monophosphate", "adenosine_5-monophosphate", "glutamate", "adenosine","tyrosine"))
all_metabolite_percentage_MT$metabolite <- factor(all_metabolite_percentage_MT$metabolite, levels = desired_order)


####Figure 5B: MT KO stacked bar####
# Create the stacked bar chart
MT_percent_graph <- ggplot(all_metabolite_percentage_MT, aes(x = metabolite, y = Value_updated, fill = Category)) +
  geom_bar(stat = "identity", color = "black", size = 1.1) +
  facet_wrap(~ row) +
  scale_fill_manual(values = custom_colors, breaks = legend_order) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +  # Ensure y-axis starts from 0
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove background color
    axis.line = element_line(colour = "black"),  # Add black axis lines
    text = element_text(size = 12),
    axis.title = element_text(size = 22), 
    plot.title = element_text(size = 24, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),  # Adjust x-axis text size 
    axis.text.x = element_text(size = 16),  # Adjust x-axis text size 
    strip.text = element_text(size = 20, face = "bold"),  # Adjust facet label text size
    panel.spacing = unit(2, "lines")  
  ) +
  labs(
    title = "Metatranscriptome",
    x = "Metabolite",
    y = "Microbial meta-transcriptome (MT) contribution (%)",
    fill = "Category"
  ) +
  coord_flip()  # Flip the axes

pdf("Figure_5B___metabolite_MT_KO_stacked_bargraph_species.pdf", width=18, height=9 )
      show(MT_percent_graph)
dev.off()

#######MG KO match#######
setwd("~/MG_KO_species") #set the directory to the raw file directory 
metabolite_ko_files <- list.files(getwd())

ii<-1
all_metabolite_percentage_MG <- data.frame(
  row = character(),        # Character column
  category = character(),   # Character column
  value = numeric(),        # Numeric column
  metabolites = character() # Character column
)

for (file in metabolite_ko_files){
  print(ii)
  metabolite_name <- sub("\\.tsv$", "", file)
  print(metabolite_name)
  metabolite_KO <- read.table(file, header = TRUE, sep = "\t")
  metabolite_KO$taxa<-rownames(metabolite_KO)
  metabolite_KO <- metabolite_KO %>% select(c("taxa", everything()))
  if (nrow(metabolite_KO) == 0 ) {
    print(paste0(file," is empty"))
  } else {
    print(file)
    MTMG_transform_data(dataframe=metabolite_KO, 
                        meta_dataframe=Mouse_map,
                        meta_dataframe_var_interest=c("Sample.ID_Segal.Lab", "Sample_Type_IS", "Innoculation"),
                        meta_dataframe_subjectID="Sample.ID_Segal.Lab",
                        output_dataframe="metabolite_KO_updated", 
                        first_column_name="taxa")
    MTMG_transform2_data (dataframe_edit=metabolite_KO_updated, 
                          sample_type_var="Sample_Type_IS",
                          sample_type_keep="Lung",
                          meta_dataframe_var_remove=c("Sample.ID_Segal.Lab", "Sample_Type_IS"),
                          output_dataframe="metabolite_KO_updated2")
    metabolite_KO_updated2[is.na(metabolite_KO_updated2)] <- 0
    metabolite_KO_updated3<-metabolite_KO_updated2[-1,-1]
    metabolite_KO_updated3 <- apply(metabolite_KO_updated3, 2, as.numeric)
    rownames(metabolite_KO_updated3)<-rownames(metabolite_KO_updated2[-1,-1])
    sample_type_KO<- metabolite_KO_updated2[1,- 1]
    # Create a new dataframe with the sums
    metabolite_KO_updated4 <- data.frame(
      MOC_Lung = rowSums(metabolite_KO_updated3[,which(sample_type_KO == "MOC_Lung")]),
      PBS_Lung = rowSums(metabolite_KO_updated3[,which(sample_type_KO == "PBS_Lung")]),
      MOC_Oral = rowSums(metabolite_KO_updated3[,which(sample_type_KO == "MOC_Oral")])
    )
    metabolite_KO_updated4 <- metabolite_KO_updated4 %>% select(c("MOC_Lung","PBS_Lung"))
    
    metabolite_KO_updated4<<-metabolite_KO_updated4
    
    #adding the barplots 
    simplified_metabolite_name <- gsub("\\.species", "", metabolite_name)
    name_output<-paste0("MG KEGG ortho percentage: ",simplified_metabolite_name)
    MG_simplified_metabolite_name<-paste0("MG_",simplified_metabolite_name)
    barplot_KO_metabolites(metabolite_KO_updated4,name_output,MG_simplified_metabolite_name)
    
    if (all(metabolite_KO_updated4 == 0)){ #only doing the graph for ones that are not all zeros
      print(paste0(file," is all zeros"))
    } else {
      # Rows to exclude
      Prevotella <- c("Prevotella.melaninogenica")
      Streptococcus <- c("Streptococcus.mitis")
      Veillonella <- c("Veillonella.parvula")
      exclude_rows<-c("Prevotella.melaninogenica", "Streptococcus.mitis", "Veillonella.parvula")
      # Sum rows of MOC
      sum_rows_Prevotella <- data.frame(t(colSums(metabolite_KO_updated4[(rownames(metabolite_KO_updated4) %in% Prevotella), , drop = FALSE])))
      rownames(sum_rows_Prevotella)<-"Prevotella.melaninogenica"
      sum_rows_Streptococcus <- data.frame(t(colSums(metabolite_KO_updated4[(rownames(metabolite_KO_updated4) %in% Streptococcus), , drop = FALSE])))
      rownames(sum_rows_Streptococcus)<-"Streptococcus.mitis"
      sum_rows_Veillonella <- data.frame(t(colSums(metabolite_KO_updated4[(rownames(metabolite_KO_updated4) %in% Veillonella), , drop = FALSE])))
      rownames(sum_rows_Veillonella)<-"Veillonella.parvula"
      # Sum rows excluding specified ones
      sum_rows_other <- data.frame(t(colSums(metabolite_KO_updated4[!(rownames(metabolite_KO_updated4) %in% exclude_rows), , drop = FALSE])))
      rownames(sum_rows_other)<-"others"
      combined_other<- rbind(metabolite_KO_updated4[(rownames(metabolite_KO_updated4) %in% exclude_rows), , drop = FALSE], sum_rows_other)
      combined_MOC<- rbind(sum_rows_Prevotella, sum_rows_Streptococcus, sum_rows_Veillonella, sum_rows_other)
      combined_MOC_t<-t(combined_MOC)
      
      ########relative abundance
      combined_MOC_t_relab<- as.data.frame(t(apply(combined_MOC_t, 1, function(row) row/sum(row))))
      
      combined_MOC_t_relab$row <- rownames(combined_MOC_t)
      combined_MOC_t_relab<<-combined_MOC_t_relab
      
      df_long <- pivot_longer(combined_MOC_t_relab, cols = -row, names_to = "Category", values_to = "Value")
      print(metabolite_KO_updated4)
      print(combined_MOC_t)
      print(df_long)
      
      df_long$metabolite<- simplified_metabolite_name
      
      all_metabolite_percentage_MG<- rbind(all_metabolite_percentage_MG, df_long)
      
      ii<-ii+1
    } 
  }
}

all_metabolite_percentage_MG <- all_metabolite_percentage_MG %>%
  mutate(Value = ifelse(is.na(Value), 0, Value))

all_metabolite_percentage_MG$Value_updated<-round(all_metabolite_percentage_MG$Value * 100, 2)

custom_colors <- c("Prevotella.melaninogenica" = "#E4B9A1", 
                   "Streptococcus.mitis" = "#D4F1AC", 
                   "Veillonella.parvula" = "#C9A9DC", 
                   "others" = "grey95",
                   "NA"="white")


legend_order <- c("Prevotella.melaninogenica", "Streptococcus.mitis", "Veillonella.parvula", "others", "NA")

all_metabolite_percentage_MG<- all_metabolite_percentage_MG[all_metabolite_percentage_MG$metabolite %in% 
                                                              c("l-glutamate","l-tyrosine",  
                                                                "l-methionine", "adenosine", "niacinamide",
                                                                "adenosine_5-monophosphate", "inosine_5-monophosphate"),]

#setting PBS_Lung on the left side and then MOC_Lung on the right side 
all_metabolite_percentage_MG$row<-factor(all_metabolite_percentage_MG$row, levels=c("PBS_Lung", "MOC_Lung"))

all_metabolite_percentage_MG$metabolite <-gsub("^l-|^d-|^dl-", "", all_metabolite_percentage_MG$metabolite) #remove "l-", "dl-", "d-" from the label

desired_order <- rev(c("niacinamide","methionine", "inosine_5-monophosphate", "adenosine_5-monophosphate", "glutamate", "adenosine","tyrosine"))
all_metabolite_percentage_MG$metabolite <- factor(all_metabolite_percentage_MG$metabolite, levels = desired_order)

####Figure 5C: MG KO stacked bar####
# Create the stacked bar chart
MG_percent_graph <- ggplot(all_metabolite_percentage_MG, aes(x = metabolite, y = Value_updated, fill = Category)) +
  geom_bar(stat = "identity", color = "black", size = 1.1) +
  facet_wrap(~ row) +
  scale_fill_manual(values = custom_colors, breaks = legend_order) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +  # Ensure y-axis starts from 0
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove background color
    axis.line = element_line(colour = "black"),  # Add black axis lines
    text = element_text(size = 12),
    axis.title = element_text(size = 22), 
    plot.title = element_text(size = 24, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),  # Adjust x-axis text size 
    axis.text.x = element_text(size = 16),  # Adjust x-axis text size 
    strip.text = element_text(size = 20, face = "bold"),  # Adjust facet label text size
    panel.spacing = unit(2, "lines")  
  ) +
  labs(
    title = "Metagenome",
    x = "Metabolite",
    y = "Microbial metagenome (MG) contribution (%)",
    fill = "Category"
  ) +
  coord_flip()  # Flip the axes

pdf("Figure_5C___metabolite_MG_KO_stacked_bargraph_species.pdf", width=18, height=9)
          show(MG_percent_graph)
dev.off()
   
#~~~~~~~~~~~~~~~~~~~~~~#
##### exvivo isotope percent#####
#~~~~~~~~~~~~~~~~~~~~~~#
setwd("C:/Users/kendr/Dropbox (NYU Langone Health)/Functional.V2/KW_files/COPD_Functional_output_files") 

Metabolites_exvivo_isotope_percent2 <-read.csv("Metabolites_exvivo_isotope_percent2_export.csv")
####Figure 6B: exvivo isotope relative intensity####
isotope_meta_histogram_multiple(Metabolites_exvivo_isotope_percent2,
                                metabolites=c("d-fructose 1,6-bisphosphate-13c-5","adenosine-13c-5", "inosine-13c-5", "l-glutamate-13c-3","l-methionine-13c-5_1"), 
                                metabolite_final_name=c("d-fructose 1,6-bisphosphate","adenosine", "inosine", "glutamate", "methionine"),
                                compare_group1=c("C12.Bacteria.1", "C12.Bacteria.2","C12.Bacteria.3"),
                                compare_group2=c("C13.Bacteria.1","C13.Bacteria.2","C13.Bacteria.3"),
                                group1_name="C12 bacteria", group2_name="C13 bacteria",
                                group1_color="#709666",
                                group2_color="#6FB3BD",
                                sample_types=c("C12 bacteria", "C13 bacteria"), 
                                xlabel_size = 10, ylabel_size = 20, axis_title_size = 23,
                                width=15,
                                xlabel="Metabolite", ylabel="Relative Intensity (log 10 of percent)", titlelabel="C13 Labeling",
                                output_name="Figure_6B___all_metabolites_exvivo_all", show_pvalue=TRUE,   log_transform = TRUE,
                                angle_labels = FALSE)


