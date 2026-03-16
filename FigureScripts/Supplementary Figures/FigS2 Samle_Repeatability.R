library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(grid)

multiplot <- function(..., plotlist = NULL, cols = 1, layout = NULL) {
  plots <- c(list(...), plotlist)
  numPlots <- length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
                     ncol = cols, nrow = ceiling(numPlots / cols))
  }
  if (numPlots == 1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

Corfigure <- function(dataFrame) {
  ggscatter(dataFrame, x = colnames(dataFrame)[1], y = colnames(dataFrame)[2],
            add = "reg.line", conf.int = TRUE, color = "#FF7C88",
            add.params = list(fill = "lightgray")) +
    stat_cor(method = "pearson", label.x = 0, label.y = 12)
}

log2p1 <- function(x) log2(x + 1)

data  <- read.table("01.Axm_TPM.txt", header = TRUE, sep = "\t")
data1 <- data[, 3:49]

make_pair <- function(L, R) as.data.frame(cbind(log2p1(L), log2p1(R)))

eFC1 <- make_pair(data1$X52_FLD1, data1$X52_FRD1)
eFC2 <- make_pair(data1$X52_FLD2, data1$X52_FRD2)
eFC3 <- make_pair(data1$X52_FLD3, data1$X52_FRD3)
eFC4 <- make_pair(data1$X52_FLD4, data1$X52_FRD4)

eHC1 <- make_pair(data1$X55_HLD1, data1$X55_HRD1)
eHC2 <- make_pair(data1$X55_HLD2, data1$X55_HRD2)
eHC3 <- make_pair(data1$X55_HLD3, data1$X55_HRD3)
eHC4 <- make_pair(data1$X55_HLD4, data1$X55_HRD4)
eHC5 <- make_pair(data1$X55_HLD5, data1$X55_HRD5)

multiplot(Corfigure(eFC1), Corfigure(eFC2), Corfigure(eFC3), Corfigure(eFC4),
          Corfigure(eHC1), Corfigure(eHC2), Corfigure(eHC3), Corfigure(eHC4), Corfigure(eHC5),
          cols = 9)

lFC1 <- make_pair(data1$X55_FRD1, data1$X55_FRD2)
lFC2 <- make_pair(data1$X57_FLD2, data1$X57_FRD2)
lFC3 <- make_pair(data1$X55_FLD3, data1$X55_FRD3)
lFC4 <- make_pair(data1$X55_FLD4, data1$X55_FRD4)

lHC1 <- make_pair(data1$X57_HLD1, data1$X57_HRD1)
lHC2 <- make_pair(data1$X57_HLD2, data1$X57_HRD2)
lHC3 <- make_pair(data1$X57_HLD3, data1$X57_HRD3)
lHC4 <- make_pair(data1$X57_HLD4, data1$X57_HRD4)
lHC5 <- make_pair(data1$X57_HLD5, data1$X57_HRD5)

multiplot(Corfigure(lFC1), Corfigure(lFC2), Corfigure(lFC3), Corfigure(lFC4),
          Corfigure(lHC1), Corfigure(lHC2), Corfigure(lHC3), Corfigure(lHC4), Corfigure(lHC5),
          cols = 9)