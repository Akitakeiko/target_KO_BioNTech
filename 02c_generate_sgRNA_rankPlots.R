

rm(list=ls())

library(seqinr)
library(readxl)
library(data.table)
library(magrittr)
library(janitor)
library(dplyr)
library(tibble)
library(tidyr)
library(e1071)
library(corrplot)
library(cowplot) 
library(pheatmap)
library(plotly)
library(ggrepel)
library(MAGeCKFlute)
library(pathview)
library(ggplot2)
library(xtable)  # threw error at this line (when running in terminal)

dir_rra = "/mnt/efs/dmanandhar/research/2022-04-Tcellscreen/tcellscreen/outputs/01_PAN_SCREEN_ANALYSES/mageck_outputs/rra"
wave_12_genes = c("PRDM1", "TNFAIP3", "ZC3H12A", "SOCS1", "PTPN2", "CISH")

csv_contrasts = "/mnt/efs/dmanandhar/research/2022-04-Tcellscreen/tcellscreen/scripts/01_PAN_SCREEN_ANALYSES/intermediate_files/df_contrasts_all.wSCCs.Nov24_2024.csv"
outputdir = "/mnt/efs/cgu"


# sgrna rank plot
plot_sgrna_ranks = function(sdata_df, 
                            title="", 
                            subtitle = "", 
                            gene = c(), #c("PRDM1"),
                            xlim = NULL, 
                            topN = 3, 
                            bottomN=3,
                            bg_gene="Intergenic"){
  if (is.null(xlim)){
    xlim_min = sdata_df$LFC %>% min() %>% floor()
    xlim_max = sdata_df$LFC %>% max() %>% ceiling()
  }else{
    xlim_min = xlim[1]
    xlim_max = xlim[2]
  }
  
  neg_ctrl_gene_only = bg_gene %>% gsub("_.*", "", .)
  
  s_hist = sdata_df %>% 
    filter(grepl(neg_ctrl_gene_only, sgrna)) %>% 
    ggplot(aes(x=LFC)) + 
    geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "lightgray", bins=100) + 
    geom_density(lwd = 0.7, linetype = 1, colour = "#252525") +
    xlim(c(xlim_min, xlim_max)) + 
    labs(title = title,
         subtitle = subtitle, 
         y = "Density") + 
    theme(plot.subtitle = element_text(size=6), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.margin=margin(b=-4,unit="mm"))
  
  
  
  
  s_ranks = sgRankView(sdata_df, top=topN, bottom = bottomN, #top = 7, bottom = 7,
                       neg_ctrl = bg_gene,
                       interval = 0.06, gene=gene) +
    xlim(c(xlim_min, xlim_max))
  
  final_plot = plot_grid(s_hist, s_ranks, ncol=1, align='v', rel_heights = c(2,5), vjust = 1.5)
  return (final_plot)
}



get_title = function(df_contrasts, acontrast_id="contrast_0"){
  df_ = df_contrasts %>% 
    dplyr::filter(contrast_id == acontrast_id)
  title = paste0(df_$screen_condB,"||", df_$condB, " (n=", df_$n_condB_samples, ")\n-vs-\n", df_$screen_condA, "||", df_$condA, " (n=", df_$n_condA_samples, ")")
  return(title)
}

get_sample_level_title = function(df_contrasts, acontrast_id="contrast_0", rep_id = "r0"){
  rep_index = gsub("r", "", rep_id) %>% as.integer() 
  rep_index = rep_index + 1
  
  df_ = df_contrasts %>% dplyr::filter(contrast_id == acontrast_id)
  condB_samples = gsub(" \\| ", ",", df_$condB_samples)
  condA_samples = gsub(" \\| ", ",", df_$condA_samples)
  
  condB_samples = unlist(strsplit(condB_samples, ",")) %>% as.character()
  condA_samples = unlist(strsplit(condA_samples, ",")) %>% as.character()
  
  return (paste0(condB_samples[rep_index], "\n-vs-\n", condA_samples[rep_index]))
}



make_sgrna_rankPlots = function(df_contrasts, f="contrast_30", topN=10, bottomN=10, include_w12genes = FALSE){
  sf = paste0(dir_rra, "/", f, "/", f, ".sgrna_summary.txt") # for eg. f == "contrast_75" ; sf = sgRNA summary file
  gf = paste0(dir_rra, "/", f, "/", f, ".gene_summary.txt")
  df_ = df_contrasts %>% dplyr::filter(contrast_id == f)
  
  sdata_df = read.csv(sf, sep = "\t") # %>% dplyr::select(sgrna, Gene, LFC)
  gdata_df = read.csv(gf, sep = "\t") %>% 
    rename(Gene = "id")
  if (include_w12genes){
    gdata_df = gdata_df %>% filter((`neg.rank` <= bottomN) | (`pos.rank` <= topN)) # note: gf does not have Intergenic gene
  }else{
    gdata_df = gdata_df %>% filter((`neg.rank` <= bottomN) | (`pos.rank` <= topN))  # note: gf does not have Intergenic gene
  }
  gdata_df = gdata_df %>% 
    mutate(rra_rank = `pos.rank`) %>%
    mutate(Gene_wRanks = paste0(Gene, " (", `pos.rank`, "p)") )

  # gdata_df = gdata_df %>%
  #   #mutate(rra_rank = ifelse(`neg.rank` <= bottomN, `neg.rank`, `pos.rank`)) %>% 
  #   mutate(Gene_wRanks = case_when(`neg.rank` <= bottomN ~ paste0(Gene, " (", `neg.rank`, "n)"), 
  #                                  `pos.rank` <= topN ~ paste0(Gene, " (", `pos.rank`, "p)"), 
  #                                  TRUE ~ paste0(Gene, " (", `pos.rank`, "p)") ))
  
  
  # update sdata_df with gene rank
  sdata_df_ = sdata_df %>% 
    right_join( bind_rows(gdata_df %>% dplyr::select("Gene", "Gene_wRanks"), 
                          data.frame("Gene" = c("Intergenic"), "Gene_wRanks" = c("Intergenic")) )) %>% # note: gf does not have Intergenic gene
    dplyr::select(sgrna, Gene_wRanks, LFC) %>% 
    rename(Gene = "Gene_wRanks")
  
  # # put wave12 genes at the bottom
  # gdata_df1 = gdata_df %>% filter(Gene %in% wave_12_genes)
  # gdata_df2 = gdata_df %>% filter(!(Gene %in% wave_12_genes))
  # plot_genes = c(gdata_df1 %>% arrange(`neg.rank`) %>% dplyr::select(Gene_wRanks) %>% deframe(),
  #                gdata_df2 %>% arrange(`neg.rank`) %>% dplyr::select(Gene_wRanks) %>% deframe())
  
  plot_genes = gdata_df %>% arrange(`pos.rank`) %>% dplyr::select(Gene_wRanks) %>% deframe()
  plot_genes = plot_genes[! plot_genes %in% c("Intergenic", "INTERGENIC")]

  
  main_plot = plot_sgrna_ranks(sdata_df_ %>% arrange(Gene), 
                               gene = plot_genes,
                               topN = topN, bottomN = bottomN, subtitle = get_title(df_contrasts, acontrast_id=f))
  
  if (df_$is_paired %in% c(FALSE, "False", "FALSE")){
    ggsave(paste0(outputdir, "/", f, "/", f, ".sgrna_rankplot.png"), width = 5, height = 6.5,units = "in", device='png') 
  }else{
    
    # identify the total number of replicates
    replicates = c()
    for (x in sdata_df$sgrna){
      arep = unlist(strsplit(x, split = "_r"))[-1]
      replicates = c(replicates, paste0("r", arep)) 
    }
    replicates = replicates %>% unique() %>% sort()
    
    
    # make plots for each replicate
    rep_plots = list()
    rep_plots[['main']] = main_plot
    for (arep in replicates){
      sdata_df_r = sdata_df_ %>% filter(grepl(arep, sgrna)) %>% arrange(desc(Gene))
      arep_plot = plot_sgrna_ranks(sdata_df_r, 
                                   gene = plot_genes,
                                   subtitle = get_sample_level_title(df_contrasts, acontrast_id=f, rep_id = arep))
      rep_plots[[arep]] = arep_plot
    }
    
    # put them all together
    if (length(replicates) > 5){
      combined_plot = plot_grid(plotlist = rep_plots, nrow=as.integer(length(rep_plots)/5)+1, align='hv') 
      
      png(paste0(dir_rra, "/", f, "/", f, ".sgrna_rankplot.png"), width = 5*5, height = 6.5* (as.integer(length(rep_plots)/5)+1), units = "in", res = 200)
      print(combined_plot)
      dev.off()
      
    }else{
      combined_plot = plot_grid(plotlist = rep_plots, ncol=rep_plots %>% length(), align='h')    
      
      png(paste0(dir_rra, "/", f, "/", f, ".sgrna_rankplot.png"), width = 6.5*length(replicates), height = 6, units = "in", res = 200)
      print(combined_plot)
      dev.off()
      
    }
  }
}


df_contrasts = read.csv(csv_contrasts, sep = ",")

for (f in c("contrast_164", # poly ras 
            # "contrast_203", # prame skmel ras
            # "contrast_219", # prame 624mel ras
          
            
            "contrast_30", # poly stemness 
            # "contrast_136",  # 624mel stemness
            # "contrast_112", # skmel stemness
            # "contrast_84", # bnt211 stemness
            # "contrast_158", # neostim stemness
            
            "contrast_15" # poly prolif
            # "contrast_128", # 624mel prolif
            # "contrast_108", # skmel prolif
            # "contrast_72", # bnt221 prolif
            # "contrast_155" # neostim prolif
            )){
  make_sgrna_rankPlots(df_contrasts, f=f, topN=20, bottomN=0)  
}



#df_contrasts %>% head()
#df_ = df_contrasts %>% dplyr::filter(contrast_id == f)

# ==================================================================================================================================


options(warn=-1)

for (f in list.files(dir_rra)){
  tryCatch(
    expr = {
      sf = paste0(dir_rra, "/", f, "/", f, ".sgrna_summary.txt") # for eg. f == "contrast_75" ; sf = sgRNA summary file
      sdata_df = read.csv(sf, sep = "\t") %>% dplyr::select(sgrna, Gene, LFC)
      sgrna_rp = plot_sgrna_ranks(sdata_df, gene = wave_12_genes, topN = topN, bottomN = bottomN, subtitle = get_title(df_contrasts, acontrast_id=f))
      ggsave(paste0(outputdir, "/", f, "/", f, ".sgrna_rankplot.png"), width = 7, height = 11, dpi = 180, units = "in", device='png') 
    },
    error = function(e){ 
      print(e)
    },
    finally = {
      print(paste0("done with ", f))
    }
  )
}
