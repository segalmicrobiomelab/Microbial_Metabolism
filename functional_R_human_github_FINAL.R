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

unique_metabolite_data_filtered_final<-read.csv("unique_metabolite_data_filtered_final.csv")

unique_metabolite_data_filtered_noBKG<- unique_metabolite_data_filtered_final[,!colnames(unique_metabolite_data_filtered_final) %in% c("COPD.0001.BKG.05.12.2015",
                                                                                                "SmNV.0004.BKG.06.15.2015",
                                                                                                "SmNV.0005.BKG.06.17.2015",
                                                                                                "COPD.0003.BKG.06.18.2015",
                                                                                                "SmNV.0008.BKG.06.22.2015",
                                                                                                "COPD.0007.BKG.06.29.2015",
                                                                                                "SmNV.0013.BKG.01.13.2016",
                                                                                                "COPD.0021.BKG.04.14.2016",
                                                                                                "COPD.0026.BKG.06.23.2016",
                                                                                                "COPD.0030.BKG.09.15.2016")]

####Figure 1A####
beta_diversity_nonPS_KW(input_table=unique_metabolite_data_filtered_final,
                        sample_type_var_name="Label", 
                        sample_types=c("BAL", "UA", "BKG"), 
                        sample_type_color=c("blue", "orange", "ivory3"), 
                        p_value="yes", 
                        turn_relative_ab="yes",
                        output_name="Figure_1A___PCA_COPD_functional_human_metabolites_FILTERED",
                        p_value_location="BR",
                        output_format="pdf")

beta_diversity_nonPS_KW(input_table=unique_metabolite_data_filtered_final,
                        sample_type_var_name="Label", 
                        sample_types=c("BAL", "BKG"), 
                        sample_type_color=c("blue", "ivory3"), 
                        p_value="yes", 
                        turn_relative_ab="yes",
                        output_name="PCA_COPD_functional_human_metabolites_lower_BKG_FILTERED",
                        p_value_location="BL",
                        output_format="pdf")

beta_diversity_nonPS_KW(input_table=unique_metabolite_data_filtered_final,
                        sample_type_var_name="Label", 
                        sample_types=c("UA", "BKG"), 
                        sample_type_color=c("orange", "ivory3"), 
                        p_value="yes", 
                        turn_relative_ab="yes",
                        output_name="PCA_COPD_functional_human_metabolites_upper_BKG_FILTERED",
                        p_value_location="BL",
                        output_format="pdf")

beta_diversity_nonPS_KW(input_table=unique_metabolite_data_filtered_final,
                        sample_type_var_name="Label", 
                        sample_types=c("BAL", "UA"), 
                        sample_type_color=c("blue", "orange"), 
                        p_value="yes", 
                        turn_relative_ab="yes",
                        output_name="PCA_COPD_functional_human_metabolites_lower_upper_FILTERED",
                        p_value_location="BL",
                        output_format="pdf")

####Supp Figure 1####
unique_metabolite_data_filtered_final$Skeleton_Metabolite<-gsub("^l-|^d-|^dl-", "", unique_metabolite_data_filtered_final$Skeleton_Metabolite) #remove "l-", "dl-", "d-" from the label
boxplot_nonPS_KW(input_table=unique_metabolite_data_filtered_final, 
                 sample_type_var_name="Label", 
                 sample_types=c("BAL", "UA", "BKG"),
                 sample_type_color=c("blue", "orange", "ivory3"), 
                 sample_type_color_2nd=c("darkblue","darkorange", "darkgrey"), 
                 compare_stat_option="mean",
                 log_scale="yes",
                 graph_option="boxplot",
                 display_rank="yes",
                 output_prefix="SFigure_1___COPD_functional_metabolite_human_FILTERED")

####Figure 1B####
unique_metabolite_data_filtered_noBKG$Skeleton_Metabolite<-gsub("^l-|^d-|^dl-", "", unique_metabolite_data_filtered_noBKG$Skeleton_Metabolite) #remove "l-", "dl-", "d-" from the label
PLS_DA_nonPS_KW(input_table=unique_metabolite_data_filtered_noBKG, 
                sample_type_var_name="Label", 
                sample_types=c("UA","BAL"), 
                sample_type_color=c( "orange","blue"),   
                top_display_loading=25, 
                flip_loading_direction="no",
                output_name="Figure_1B___PLSDA_COPD_functional_metabolite_lower_upper_FILTERED",
                output_format="pdf",
                width=10,
                height=15)

###MT human data

####MG####
####MG- ortho####
MG_count_ko_ortho__human_final2<-read.csv("MG_count_ko_ortho__human_final2.csv")
Edger_nonPS_KW(input_table=MG_count_ko_ortho__human_final2, 
               sample_type_var_name="Sample_Type_IS", 
               sample_types=c("UA", "BAL"), 
               sample_type_color=c("orange","blue"), 
               FDR_cut_off=0.2,
               number_display=20,
               plot_title="EdgeR: Metagenome KO- BAL vs UA",
               output_name="EdgeR_COPD_functional_MG_KO_ortho_human",
               legend_onplot="yes",
               width=12)

####MG- module####
MG_count_ko_module__human_final2<-read.csv("MG_count_ko_module__human_final2.csv")
Edger_nonPS_KW(input_table=MG_count_ko_module__human_final2, 
               sample_type_var_name="Sample_Type_IS", 
               sample_types=c("UA", "BAL"), 
               sample_type_color=c("orange","blue"), 
               FDR_cut_off=0.2,
               number_display=20,
               plot_title="EdgeR: Metagenome module- BAL vs UA",
               output_name="EdgeR_COPD_functional_MG_module_human",
               legend_onplot="yes",
               width=15)

####MG- pathway####
MG_count_ko_pathway__human_final_FILTERED<-read.csv("MG_count_ko_pathway__human_final_FILTERED.csv")
####Figure 1C####
Edger_nonPS_KW(input_table=MG_count_ko_pathway__human_final_FILTERED, 
               sample_type_var_name="Sample_Type_IS", 
               sample_types=c("UA", "BAL"), 
               sample_type_color=c("orange","blue"), 
               FDR_cut_off=0.2,
               number_display=20,
               graph_option="lollipop",
               plot_title="EdgeR: Metagenome pathway- BAL vs UA",
               output_name="Figure_1C___EdgeR_COPD_functional_MG_pathway_human_all",
               legend_onplot="yes",
               width=12)

              #removing photosythnesis from the results
              res <- subset(res_export_removeresults, Name != "map00196.Photosynthesis...antenna.proteins")
              res <- subset(res, Name != "map00100.Steroid.biosynthesis")
              res <- subset(res, Name != "map00140.Steroid.hormone.biosynthesis")
              res <- subset(res, Name != "map00121: Secondary bile acid biosynthesis")
              res <- subset(res, Name != "map00908.Zeatin.biosynthesis")
              p<-ggplot(res, aes(y=reorder(Name,-start), x=logFC,fill=col,size=abundance)) +
                geom_point(color="black",alpha=0.8,shape=21)+
                geom_segment(data=res[res$adj.P.Val<0.2,],aes(yend=reorder(Name,-start)), xend=(-30), color= "black", linetype = "solid",linewidth=1)+ 
                scale_fill_manual(values=c("B"="blue","A"="orange","D"="white"))+ 
                scale_size_continuous(name="Relative Abundance",range=c(5, 20))+                
                ggtitle("EdgeR: Metagenome pathway- BAL vs UA")+
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
              pdf_output<-paste0("Figure_1C___EdgeR_COPD_functional_MG_pathway_human_FINAL",".pdf")
              pdf(pdf_output, width=12, height=8)
              show(p)
              dev.off()

#MG pathway alpha and beta 
MG_count_ko_pathway__human_final2_all<-read.csv("MG_count_ko_pathway__human_final2_all.csv")
####Supp Figure 2A####
alpha_diversity_nonPS_KW(input_table=MG_count_ko_pathway__human_final2_all,
                         sample_type_var_name="Sample_Type_IS", 
                         sample_types=c("BKG","BAL", "UA"), 
                         sample_type_color=c("ivory3", "blue", "orange"), 
                         p_value="yes", width=5, height=8, 
                         output_name="SFigure_2A___alpha_COPD_functional_human_MG_all",
                         output_format="pdf")
####Supp Figure 2B####
beta_diversity_nonPS_KW(input_table=MG_count_ko_pathway__human_final2_all,
                        sample_type_var_name="Sample_Type_IS", 
                        sample_types=c("BAL", "UA", "BKG"), 
                        sample_type_color=c("blue", "orange", "ivory3"), 
                        p_value="yes", 
                        turn_relative_ab="yes",
                        output_name="SFigure_2B___PCA_COPD_functional_human_MG_all",
                        p_value_location="TR", plot_title="MG pathway",
                        output_format="pdf")
beta_diversity_nonPS_KW(input_table=MG_count_ko_pathway__human_final2_all,
                        sample_type_var_name="Sample_Type_IS", 
                        sample_types=c("BAL", "UA"), 
                        sample_type_color=c("blue", "orange"), 
                        p_value="yes", 
                        turn_relative_ab="yes",
                        output_name="PCA_COPD_functional_human_MG_lower_upper",
                        p_value_location="TR",
                        output_format="pdf")
beta_diversity_nonPS_KW(input_table=MG_count_ko_pathway__human_final2_all,
                        sample_type_var_name="Sample_Type_IS", 
                        sample_types=c("BAL", "BKG"), 
                        sample_type_color=c("blue", "ivory3"), 
                        p_value="yes", 
                        turn_relative_ab="yes",
                        output_name="PCA_COPD_functional_human_MG_lower_BKG",
                        p_value_location="TR",
                        output_format="pdf")
beta_diversity_nonPS_KW(input_table=MG_count_ko_pathway__human_final2_all,
                        sample_type_var_name="Sample_Type_IS", 
                        sample_types=c("UA", "BKG"), 
                        sample_type_color=c("orange", "ivory3"), 
                        p_value="yes", 
                        turn_relative_ab="yes",
                        output_name="PCA_COPD_functional_human_MG_upper_BKG",
                        p_value_location="BL",
                        output_format="pdf")
####Supp Figure 5A####
heatmap_nonPS_KW(input_table=MG_count_ko_pathway__human_final2_all,
                 sample_type_var_name="Sample_Type_IS", 
                 sample_types=c("BAL", "UA", "BKG"), 
                 sample_type_color=c("blue", "orange","ivory3"), 
                 plot_title="MG: BAL, UA, BKG ",
                 top_N_rows=100,
                 width=20, height=15,
                 xaxis_font=2,yaxis_font=1,
                 xaxis_margin=12,yaxis_margin=40,
                 color_palette_set=c("white","darkblue"),
                 log_value="yes",
                 color_palette_num=100,
                 legend_position="topright",
                 output_name="SFigure_5A___COPD_functional_human_MG_heatmap")

####MT####
####MT- ortho####
MT_count_ko_ortho__human_final2<-read.csv("MT_count_ko_ortho__human_final2.csv")
Edger_nonPS_KW(input_table=MT_count_ko_ortho__human_final2, 
               sample_type_var_name="Sample_Type_IS", 
               sample_types=c("UA", "BAL"), 
               sample_type_color=c("orange","blue"), 
               FDR_cut_off=0.2,
               number_display=20,
               plot_title="EdgeR: Metatranscriptome KO- BAL vs UA",
               output_name="EdgeR_COPD_functional_MT_ko_ortho_human",
               legend_onplot="yes",
               width=18.5)
####MT- module####
MT_count_ko_module__human_final2<-read.csv("MT_count_ko_module__human_final2.csv")
Edger_nonPS_KW(input_table=MT_count_ko_module__human_final2, 
               sample_type_var_name="Sample_Type_IS", 
               sample_types=c("UA", "BAL"), 
               sample_type_color=c("orange","blue"), 
               FDR_cut_off=0.2,
               number_display=20,
               plot_title="EdgeR: Metatranscriptome module- BAL vs UA",
               output_name="EdgeR_COPD_functional_MT_module_human",
               legend_onplot="yes",
               width=15)

####MT- pathway####
MT_count_ko_pathway__human_final_FILTERED<-read.csv("MT_count_ko_pathway__human_final_FILTERED.csv")
####Figure 1D####
Edger_nonPS_KW(input_table=MT_count_ko_pathway__human_final_FILTERED, 
               sample_type_var_name="Sample_Type_IS", 
               sample_types=c("UA", "BAL"), 
               sample_type_color=c("orange","blue"), 
               FDR_cut_off=0.2,
               number_display=20,
               graph_option="lollipop",
               plot_title="EdgeR: Metatranscriptome pathway- BAL vs UA",
               output_name="Figure_1D___EdgeR_COPD_functional_MT_pathway_human_all",
               legend_onplot="yes",
               width=12)

                #removing photosythnesis from the results
                res <- subset(res_export_removeresults, Name != "map00196.Photosynthesis...antenna.proteins")
                res <- subset(res, Name != "map00195.Photosynthesis")
                res <- subset(res, Name != "map00590.Arachidonic.acid.metabolism")
                p<-ggplot(res, aes(y=reorder(Name,-start), x=logFC,fill=col,size=abundance)) +
                  geom_point(color="black",alpha=0.8,shape=21)+
                  geom_segment(data=res[res$adj.P.Val<0.2,],aes(yend=reorder(Name,-start)), xend=(-30), color= "black", linetype = "solid",linewidth=1)+ 
                  scale_fill_manual(values=c("B"="blue","A"="orange","D"="white"))+ 
                  scale_size_continuous(name="Relative Abundance",range=c(5, 20))+                
                  ggtitle("EdgeR: Metatranscriptome pathway- BAL vs UA")+
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
                pdf_output<-paste0("Figure_1D___EdgeR_COPD_functional_MT_pathway_human_FINAL",".pdf")
                    pdf(pdf_output, width=12, height=8)
                        show(p)
                    dev.off()

#MT pathway alpha and beta 
MT_count_ko_pathway__human_final2_all<-read.csv("MT_count_ko_pathway__human_final2_all.csv")
####Supp Figure 3A####
alpha_diversity_nonPS_KW(input_table=MT_count_ko_pathway__human_final2_all,
                         sample_type_var_name="Sample_Type_IS", 
                         sample_types=c("BKG","BAL", "UA"), 
                         sample_type_color=c("ivory3", "blue", "orange"), 
                         p_value="yes", width=4, height=6, 
                         output_name="SFigure_3A___alpha_COPD_functional_human_MT_all",
                         output_format="pdf")
####Supp Figure 3B####
beta_diversity_nonPS_KW(input_table=MT_count_ko_pathway__human_final2_all,
                        sample_type_var_name="Sample_Type_IS", 
                        sample_types=c("BAL", "UA", "BKG"), 
                        sample_type_color=c("blue", "orange", "ivory3"), 
                        p_value="yes", 
                        turn_relative_ab="yes",
                        output_name="SFigure_3B___PCA_COPD_functional_human_MT_all",
                        p_value_location="TR", plot_title="MT pathway",
                        output_format="pdf")
beta_diversity_nonPS_KW(input_table=MT_count_ko_pathway__human_final2_all,
                        sample_type_var_name="Sample_Type_IS", 
                        sample_types=c("BAL", "UA"), 
                        sample_type_color=c("blue", "orange"), 
                        p_value="yes", 
                        turn_relative_ab="yes",
                        output_name="PCA_COPD_functional_human_MT_lower_upper",
                        p_value_location="TR",
                        output_format="pdf")
beta_diversity_nonPS_KW(input_table=MT_count_ko_pathway__human_final2_all,
                        sample_type_var_name="Sample_Type_IS", 
                        sample_types=c("BAL", "BKG"), 
                        sample_type_color=c("blue", "ivory3"), 
                        p_value="yes", 
                        turn_relative_ab="yes",
                        output_name="PCA_COPD_functional_human_MT_lower_BKG",
                        p_value_location="TR",
                        output_format="pdf")
beta_diversity_nonPS_KW(input_table=MT_count_ko_pathway__human_final2_all,
                        sample_type_var_name="Sample_Type_IS", 
                        sample_types=c("UA", "BKG"), 
                        sample_type_color=c("orange", "ivory3"), 
                        p_value="yes", 
                        turn_relative_ab="yes",
                        output_name="PCA_COPD_functional_human_MT_upper_BKG",
                        p_value_location="BL",
                        output_format="pdf")

####Supp Figure 5B####
heatmap_nonPS_KW(input_table=MT_count_ko_pathway__human_final2_all,
                 sample_type_var_name="Sample_Type_IS", 
                 sample_types=c("BAL", "UA", "BKG"), 
                 sample_type_color=c("blue", "orange","ivory3"), 
                 plot_title="MT: BAL, UA, BKG ",
                 top_N_rows=100,
                 width=20, height=15,
                 xaxis_font=2,yaxis_font=1,
                 xaxis_margin=12,yaxis_margin=40,
                 color_palette_set=c("white","darkblue"),
                 log_value="yes",
                 color_palette_num=100,
                 legend_position="topright",
                 output_name="SFigure_5B___COPD_functional_human_MT_heatmap")


####TABLE 1####
COPD_human_map_unique<-read.csv("COPD_human_map_unique.csv")

overall_table1_COPD<-table1(~ age_new + age_cat + male+  bmi + bmi_cat +
                              Race_cat + Smoking_status + Pack_Years +
                              Inhaled_Steroids + Inhaled_Beta_Agonist + 
                              FEV1 + FEV1_pre_perc + FEV1_post + FEV1_post_perc +
                              FVC + FVC_pre_perc + FVC_post + FVC_post_perc + 
                              R5_pre + R5_post ,
                                            data=COPD_human_map_unique, overall= "Total", render=render.NEW, caption="Table 1 by highest PGD_v3")
output_table1("overall_table1_COPD")

      



