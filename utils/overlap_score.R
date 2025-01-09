if (!require("pacman"))
  install.packages("pacman")
pacman::p_load("data.table",
               "tidyverse",
               "reshape2",
               "plyr"
)

# t1: table with 2 columns: coembed labels, raw labels
# t2: table with 2 columns: coembed labels, raw labels
cal_ovlpScore <- function(t1, t2){
  t1.table <- table(t1)
  t2.table <- table(t2)
  t1.pct <- apply(t1.table, 2, function(x){x/sum(x)})
  t2.pct <- apply(t2.table, 2, function(x){x/sum(x)})
  t1.labels <- colnames(t1.pct)
  t2.labels <- colnames(t2.pct)
  ovlpScore.df <- data.frame(anno1=as.character(), anno2=as.character(), ovlpScore=as.numeric())
  for(t1.label in t1.labels){
    for(t2.label in t2.labels){
      t1.pct.df <- data.frame(t1.pct[,t1.label])
      colnames(t1.pct.df) <- "t1"
      t1.pct.df$ident <- rownames(t1.pct.df)
      t2.pct.df <- data.frame(t2.pct[,t2.label])
      colnames(t2.pct.df) <- "t2"
      t2.pct.df$ident <- rownames(t2.pct.df)
      comp.df <- join(t1.pct.df, t2.pct.df, by="ident", type="full")
      comp.df[is.na(comp.df)] <- 0
      comp.df$ident <- NULL
      comp.df <- t(comp.df)
      ovlpScore <- sum(apply(comp.df, 2, min))
      out <- data.frame(anno1=t1.label, anno2=t2.label, ovlpScore=ovlpScore)
      ovlpScore.df <- rbind(ovlpScore.df, out)
    }
  }
  return(ovlpScore.df)
}
