options(warn=-1)
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(reshape2))

ariba_summary_csv <- "summary.fig5.csv"
data_for_plots_MICs_csv <- "fig5.mic_data.csv"

### aRiba_functions_Ngono_plots.R ###

### FUNCTION TO GET %HET IN TARGET MUTATIONS (Figure 5) ###

copyNumber <- function(ariba.table, mic.table, mutation, sec.mut=NULL, antibiotic=NULL, use.log=TRUE, nb.copies=4){

  # Get data #
  ariba <- ariba.table
  colnames(ariba) <- gsub("^X", "", colnames(ariba))
  metadata <- mic.table
  if (is.null(antibiotic)){ stop("Antibiotic name needs to be specified\n") }

  # Get columns corresponding to mutation and %het #
  g <- which(colnames(ariba)==mutation)
  ghet <- which(colnames(ariba)==paste0(mutation, "..."))
  if (length(g)>0 & length(ghet)>0){

    # Get %s #
    het <- ariba[,ghet]
    het[which(is.na(het) & ariba[,g]=="yes")] <- 100
    het[which(ariba[,g]=="no")] <- 0
    names(het) <- rownames(ariba)

    # Create data frame #
    df <- data.frame(het=het, mic=metadata[names(het),antibiotic])
    if (use.log==TRUE){ df$mic <- log2(df$mic) }

    # Check if other 23S mutations present in het==0 #
    if (!is.null(sec.mut)){
      second_mutation <- sec.mut
      second_mut <- c()
      for (x in 1:nrow(df)){
        nam <- rownames(df)[x]
        sec <- as.vector(ariba[nam,second_mutation])
        second_mut <- c(second_mut, sec)
      }
      df$second <- second_mut

      rm <- which(df$het==0 & df$second=="yes")
      if(length(rm)>0){ df <- df[-rm,] }
    }

    # Calculate correlation #
    ct <- lm(df$het~df$mic)
    cor <- cor.test(df$het, df$mic, method="spearman")
    if (cor$p.value<0.01){ pval = 0.01 }else if (cor$p.value<0.05){ pval = 0.05 }else{ pval = "ns" }
    ct.pred <- predict(ct, df[,1:2], interval="confidence")
    ct.pred <- as.data.frame(ct.pred)
    df <- cbind(df, ct.pred)

    breaks <- seq(0,100,by=10)
    bins <- findInterval(df$het, breaks)
    df$groups <- breaks[bins]

    ### PLOT ###
    accent <- brewer.pal(8, "Accent")
    paired <- brewer.pal(12, "Paired")
    azm_mics <- c(0,0.001,0.0025,0.0075, 0.015, 0.03, 0.06,0.125,0.25,0.5,1,2,4,8,16,32,64,128,256,512,1024)
    by <- 100/nb.copies
    breaks <- seq(0, 100, by=by)

    v <- round(cor$estimate,2)
    plabel <- bquote("Spearman\'s "~ rho ~"="~.(v) ~ " ") # Make label for correlation value
    plabel_text <- as.character(as.expression(plabel))

    mut <- gsub("\\.", " ", mutation)
    p3 <- ggplot(df, aes(het, mic))
    for (x in 1:length(breaks)){ p3 <- p3+geom_vline(xintercept=breaks[x], linetype="dotted", colour=accent[2]) }
    p3 <- p3+geom_count(colour=paired[10], alpha=.8)+scale_size(range=c(2,15), breaks = c(5,10,50,100,250))+
      geom_smooth(method="lm", col="#330f53")+
      xlab(paste("% reads with 23S 2597T mutation"))+ylab(bquote(.(antibiotic) ~ "log(MIC) " ~ mu * "g/mL"))+
      scale_x_continuous(breaks=seq(0,100,by=25), limits = c(0,100))+
      scale_y_continuous(breaks=log2(azm_mics), labels=azm_mics)+
      theme_bw()+
      theme(axis.text=element_text(size=14),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line.x = element_line(color="black"),
            axis.line.y = element_line(color="black"),
            axis.title = element_text(size=14),
            axis.title.x = element_text(margin = margin(t = 10)),
            legend.title = element_text(size=12, face="bold", hjust=.5),
            legend.text = element_text(size = 12))+
      labs(size="Number of\n   isolates")+
      annotate("text", x = 25, y = log2(256),
               label = paste("R^2 == ", round(summary(ct)$r.squared,2)), size=4, parse = TRUE)+
      annotate("text", x = 30, y = log2(64),
               label = paste0(plabel_text, " ", paste("(p-value < ", pval, ")")), parse=TRUE)

    ggsave(paste0("summary.", antibiotic, "_", mutation, "_copy_number.pdf"), width=7, height=5)

  }else{ stop("Mutation not found\n") }
}

## Figure 5 N. gonorrhoeae ## 23S C2611T copy number mutation ##

ariba <- read.csv(ariba_summary_csv, header=T, row.names=1)
rownames(ariba) <- gsub(".tsv", "", rownames(ariba))
rownames(ariba) <- gsub("ARIBA_reports/", "", rownames(ariba))
mics <- read.csv(data_for_plots_MICs_csv, header=T, row.names=1) # Table with MIC data

copyNumber(ariba.table=ariba, mic.table=mics, mutation="23S.23S.2597T", sec.mut="23S.23S.2045G", antibiotic="Azithromycin", nb.copies=4)

