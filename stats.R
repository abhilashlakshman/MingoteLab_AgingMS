rm(list=ls())

library(plotly)

da.3mos <- readRDS("~/OneDrive - Indiana University/ShaferLab/Others/Jackie/INTRSECT Ant.-Post._Medial-Lateral/Medial-Lateral/HeatmapData/IndData/Normalized_DA_3mos.rds")
da.12mos <- readRDS("~/OneDrive - Indiana University/ShaferLab/Others/Jackie/INTRSECT Ant.-Post._Medial-Lateral/Medial-Lateral/HeatmapData/IndData/Normalized_DA_12mos.rds")
da.24mos <- readRDS("~/OneDrive - Indiana University/ShaferLab/Others/Jackie/INTRSECT Ant.-Post._Medial-Lateral/Medial-Lateral/HeatmapData/IndData/Normalized_DA_24mos.rds")

daglu.3mos <- readRDS("~/OneDrive - Indiana University/ShaferLab/Others/Jackie/INTRSECT Ant.-Post._Medial-Lateral/Medial-Lateral/HeatmapData/IndData/Normalized_DAGLU_3mos.rds")
daglu.12mos <- readRDS("~/OneDrive - Indiana University/ShaferLab/Others/Jackie/INTRSECT Ant.-Post._Medial-Lateral/Medial-Lateral/HeatmapData/IndData/Normalized_DAGLU_12mos.rds")
daglu.24mos <- readRDS("~/OneDrive - Indiana University/ShaferLab/Others/Jackie/INTRSECT Ant.-Post._Medial-Lateral/Medial-Lateral/HeatmapData/IndData/Normalized_DAGLU_24mos.rds")

th.3mos <- readRDS("~/OneDrive - Indiana University/ShaferLab/Others/Jackie/INTRSECT Ant.-Post._Medial-Lateral/Medial-Lateral/HeatmapData/IndData/Normalized_allTH_3mos.rds")
th.12mos <- readRDS("~/OneDrive - Indiana University/ShaferLab/Others/Jackie/INTRSECT Ant.-Post._Medial-Lateral/Medial-Lateral/HeatmapData/IndData/Normalized_allTH_12mos.rds")
th.24mos <- readRDS("~/OneDrive - Indiana University/ShaferLab/Others/Jackie/INTRSECT Ant.-Post._Medial-Lateral/Medial-Lateral/HeatmapData/IndData/Normalized_allTH_24mos.rds")

da.3mos.array <- do.call(cbind, da.3mos)
da.3mos.array <- array(da.3mos.array, dim = c(dim(da.3mos[[1]]), length(da.3mos)))
da.avg.3mos <- apply(da.3mos.array, c(1,2), mean, na.rm = T)

da.12mos.array <- do.call(cbind, da.12mos)
da.12mos.array <- array(da.12mos.array, dim = c(dim(da.12mos[[1]]), length(da.12mos)))
da.avg.12mos <- apply(da.12mos.array, c(1,2), mean, na.rm = T)

da.24mos.array <- do.call(cbind, da.24mos)
da.24mos.array <- array(da.24mos.array, dim = c(dim(da.24mos[[1]]), length(da.24mos)))
da.avg.24mos <- apply(da.24mos.array, c(1,2), mean, na.rm = T)

daglu.3mos.array <- do.call(cbind, daglu.3mos)
daglu.3mos.array <- array(daglu.3mos.array, dim = c(dim(daglu.3mos[[1]]), length(daglu.3mos)))
daglu.avg.3mos <- apply(daglu.3mos.array, c(1,2), mean, na.rm = T)

daglu.12mos.array <- do.call(cbind, daglu.12mos)
daglu.12mos.array <- array(daglu.12mos.array, dim = c(dim(daglu.12mos[[1]]), length(daglu.12mos)))
daglu.avg.12mos <- apply(daglu.12mos.array, c(1,2), mean, na.rm = T)

daglu.24mos.array <- do.call(cbind, daglu.24mos)
daglu.24mos.array <- array(daglu.24mos.array, dim = c(dim(daglu.24mos[[1]]), length(daglu.24mos)))
daglu.avg.24mos <- apply(daglu.24mos.array, c(1,2), mean, na.rm = T)

th.3mos.array <- do.call(cbind, th.3mos)
th.3mos.array <- array(th.3mos.array, dim = c(dim(th.3mos[[1]]), length(th.3mos)))
th.avg.3mos <- apply(th.3mos.array, c(1,2), mean, na.rm = T)

th.12mos.array <- do.call(cbind, th.12mos)
th.12mos.array <- array(th.12mos.array, dim = c(dim(th.12mos[[1]]), length(th.12mos)))
th.avg.12mos <- apply(th.12mos.array, c(1,2), mean, na.rm = T)

th.24mos.array <- do.call(cbind, th.24mos)
th.24mos.array <- array(th.24mos.array, dim = c(dim(th.24mos[[1]]), length(th.24mos)))
th.avg.24mos <- apply(th.24mos.array, c(1,2), mean, na.rm = T)


da.diff.3.12.pval <- matrix(NA, nrow = 10, ncol = 6)
da.diff.3.24.pval <- matrix(NA, nrow = 10, ncol = 6)

daglu.diff.3.12.pval <- matrix(NA, nrow = 10, ncol = 6)
daglu.diff.3.24.pval <- matrix(NA, nrow = 10, ncol = 6)

th.diff.3.12.pval <- matrix(NA, nrow = 10, ncol = 6)
th.diff.3.24.pval <- matrix(NA, nrow = 10, ncol = 6)

for (i in 1:length(da.3mos.array[,1,1])) {
  for (j in 1:length(da.3mos.array[1,,1])) {
    x = da.3mos.array[i,j,]
    y = da.12mos.array[i,j,]
    z = da.24mos.array[i,j,]
    
    # t12 = wilcox.test(x = x, y = y)
    # t24 = wilcox.test(x = x, y = z)
    
    t12 = t.test(x = y, mu = mean(x, na.rm = T))
    t24 = t.test(x = z, mu = mean(x, na.rm = T))
    
    da.diff.3.12.pval[i,j] <- t12$p.value
    da.diff.3.24.pval[i,j] <- t24$p.value
  }
}

for (i in 1:length(daglu.3mos.array[,1,1])) {
  for (j in 1:length(daglu.3mos.array[1,,1])) {
    x = daglu.3mos.array[i,j,]
    y = daglu.12mos.array[i,j,]
    z = daglu.24mos.array[i,j,]
    
    # t12 = wilcox.test(x = x, y = y)
    # t24 = wilcox.test(x = x, y = z)
    
    t12 = t.test(x = y, mu = mean(x, na.rm = T))
    t24 = t.test(x = z, mu = mean(x, na.rm = T))
    
    daglu.diff.3.12.pval[i,j] <- t12$p.value
    daglu.diff.3.24.pval[i,j] <- t24$p.value
  }
}

for (i in 1:length(th.3mos.array[,1,1])) {
  for (j in 1:length(th.3mos.array[1,,1])) {
    x = th.3mos.array[i,j,]
    y = th.12mos.array[i,j,]
    z = th.24mos.array[i,j,]
    
    # t12 = wilcox.test(x = x, y = y)
    # t24 = wilcox.test(x = x, y = z)
    
    t12 = t.test(x = y, mu = mean(x, na.rm = T))
    t24 = t.test(x = z, mu = mean(x, na.rm = T))
    
    th.diff.3.12.pval[i,j] <- t12$p.value
    th.diff.3.24.pval[i,j] <- t24$p.value
  }
}


pre.da.diff.3.12.pval.adj <- matrix(p.adjust(da.diff.3.12.pval, method = "BH"), nrow = nrow(da.diff.3.12.pval))
pre.da.diff.3.24.pval.adj <- matrix(p.adjust(da.diff.3.24.pval, method = "BH"), nrow = nrow(da.diff.3.24.pval))

pre.daglu.diff.3.12.pval.adj <- matrix(p.adjust(daglu.diff.3.12.pval, method = "BH"), nrow = nrow(daglu.diff.3.12.pval))
pre.daglu.diff.3.24.pval.adj <- matrix(p.adjust(daglu.diff.3.24.pval, method = "BH"), nrow = nrow(daglu.diff.3.24.pval))

pre.th.diff.3.12.pval.adj <- matrix(p.adjust(th.diff.3.12.pval, method = "BH"), nrow = nrow(th.diff.3.12.pval))
pre.th.diff.3.24.pval.adj <- matrix(p.adjust(th.diff.3.24.pval, method = "BH"), nrow = nrow(th.diff.3.24.pval))

da.diff.3.12.pval.adj <- matrix(NA, nrow = nrow(pre.da.diff.3.12.pval.adj), ncol = ncol(pre.da.diff.3.12.pval.adj))
da.diff.3.24.pval.adj <- matrix(NA, nrow = nrow(pre.da.diff.3.24.pval.adj), ncol = ncol(pre.da.diff.3.24.pval.adj))
daglu.diff.3.12.pval.adj <- matrix(NA, nrow = nrow(pre.daglu.diff.3.12.pval.adj), ncol = ncol(pre.daglu.diff.3.12.pval.adj))
daglu.diff.3.24.pval.adj <- matrix(NA, nrow = nrow(pre.daglu.diff.3.24.pval.adj), ncol = ncol(pre.daglu.diff.3.24.pval.adj))
th.diff.3.12.pval.adj <- matrix(NA, nrow = nrow(pre.th.diff.3.12.pval.adj), ncol = ncol(pre.th.diff.3.12.pval.adj))
th.diff.3.24.pval.adj <- matrix(NA, nrow = nrow(pre.th.diff.3.24.pval.adj), ncol = ncol(pre.th.diff.3.24.pval.adj))

for (i in 1:length(pre.da.diff.3.12.pval.adj[,1])) {
  for (j in 1:length(pre.da.diff.3.12.pval.adj[1,])) {
    if (!is.nan(pre.da.diff.3.12.pval.adj[i,j]) & pre.da.diff.3.12.pval.adj[i,j] < 0.05) {
      da.diff.3.12.pval.adj[i,j] <- 1
    } else {
      da.diff.3.12.pval.adj[i,j] <- 0
    }
  }
}
for (i in 1:length(pre.da.diff.3.24.pval.adj[,1])) {
  for (j in 1:length(pre.da.diff.3.24.pval.adj[1,])) {
    if (!is.nan(pre.da.diff.3.24.pval.adj[i,j]) & pre.da.diff.3.24.pval.adj[i,j] < 0.05) {
      da.diff.3.24.pval.adj[i,j] <- 1
    } else {
      da.diff.3.24.pval.adj[i,j] <- 0
    }
  }
}

for (i in 1:length(pre.daglu.diff.3.12.pval.adj[,1])) {
  for (j in 1:length(pre.daglu.diff.3.12.pval.adj[1,])) {
    if (!is.nan(pre.daglu.diff.3.12.pval.adj[i,j]) & pre.daglu.diff.3.12.pval.adj[i,j] < 0.05) {
      daglu.diff.3.12.pval.adj[i,j] <- 1
    } else {
      daglu.diff.3.12.pval.adj[i,j] <- 0
    }
  }
}
for (i in 1:length(pre.daglu.diff.3.24.pval.adj[,1])) {
  for (j in 1:length(pre.daglu.diff.3.24.pval.adj[1,])) {
    if (!is.nan(pre.daglu.diff.3.24.pval.adj[i,j]) & pre.daglu.diff.3.24.pval.adj[i,j] < 0.05) {
      daglu.diff.3.24.pval.adj[i,j] <- 1
    } else {
      daglu.diff.3.24.pval.adj[i,j] <- 0
    }
  }
}

for (i in 1:length(pre.th.diff.3.12.pval.adj[,1])) {
  for (j in 1:length(pre.th.diff.3.12.pval.adj[1,])) {
    if (!is.nan(pre.th.diff.3.12.pval.adj[i,j]) & pre.th.diff.3.12.pval.adj[i,j] < 0.05) {
      th.diff.3.12.pval.adj[i,j] <- 1
    } else {
      th.diff.3.12.pval.adj[i,j] <- 0
    }
  }
}
for (i in 1:length(pre.th.diff.3.24.pval.adj[,1])) {
  for (j in 1:length(pre.th.diff.3.24.pval.adj[1,])) {
    if (!is.nan(pre.th.diff.3.24.pval.adj[i,j]) & pre.th.diff.3.24.pval.adj[i,j] < 0.05) {
      th.diff.3.24.pval.adj[i,j] <- 1
    } else {
      th.diff.3.24.pval.adj[i,j] <- 0
    }
  }
}

# CREATE FONTS
f1 <- list( # define font details
  family = "Arial, sans-serif",
  size = 16,
  color = "black"
)
f2 <- list( # define font details
  family = "Arial, sans-serif",
  size = 12,
  color = "black"
)
# CREATE x- and y-axes
x.ax <- list(
  showgrid = T,
  showline = T,
  titlefont = f1,
  tickfont = f2,
  title = "Normalized M-L",
  linecolor = "black",
  ticks = "outside",
  ticklen = 7,
  tickcolor = "black",
  range = c(0, 1)
)
y.ax <- list(
  showgrid = T,
  showline = T,
  titlefont = f1,
  tickfont = f2,
  title = "A-P (mm)",
  linecolor = "black",
  ticks = "outside",
  ticklen = 7,
  tickcolor = "black"
)

empty.plot <- plot_ly(
  x = seq(0.05, 1, by = 0.1),
  y = colnames(daglu.3mos[[1]]),
  z = t(matrix(NA, nrow = nrow(daglu.avg.3mos), ncol = ncol(daglu.avg.3mos))),
  type = "heatmap",
  zauto = F,
  zmin = 0,
  zmax = 15,
  # colors = colorRamp(c(rgb(0,0,0,1), rgb(0,1,0,1))),
  colors = colorRamp(c(rgb(1,1,1,1), rgb(1,0,1,1))),
  xgap = 0.4,
  ygap = 0.4
)%>%
  layout(
    title = "EMPTY",
    showlegend = T,
    yaxis = y.ax,
    xaxis = x.ax
  )

p.da.3mos <- plot_ly(
  x = seq(0.05, 1, by = 0.1),
  y = colnames(da.3mos[[1]]),
  z = t(da.avg.3mos),
  type = "heatmap",
  zauto = F,
  zmin = 0,
  zmax = 20,
  # colors = colorRamp(c(rgb(0,0,0,1), rgb(1,0,1,1))),
  colors = colorRamp(c(rgb(1,1,1,1), rgb(1,0,1,1))),
  xgap = 0.4,
  ygap = 0.4
)%>%
  layout(
    title = "DA-3mos",
    showlegend = T,
    yaxis = y.ax,
    xaxis = x.ax
  )

p.da.12mos <- plot_ly(
  x = seq(0.05, 1, by = 0.1),
  y = colnames(da.12mos[[1]]),
  z = t(da.avg.12mos),
  type = "heatmap",
  zauto = F,
  zmin = 0,
  zmax = 20,
  # colors = colorRamp(c(rgb(0,0,0,1), rgb(1,0,1,1))),
  colors = colorRamp(c(rgb(1,1,1,1), rgb(1,0,1,1))),
  xgap = 0.4,
  ygap = 0.4
)%>%
  layout(
    title = "DA-12mos",
    showlegend = T,
    yaxis = y.ax,
    xaxis = x.ax
  )

p.da.24mos <- plot_ly(
  x = seq(0.05, 1, by = 0.1),
  y = colnames(da.24mos[[1]]),
  z = t(da.avg.24mos),
  type = "heatmap",
  zauto = F,
  zmin = 0,
  zmax = 20,
  # colors = colorRamp(c(rgb(0,0,0,1), rgb(1,0,1,1))),
  colors = colorRamp(c(rgb(1,1,1,1), rgb(1,0,1,1))),
  xgap = 0.4,
  ygap = 0.4
)%>%
  layout(
    title = "DA-24mos",
    showlegend = T,
    yaxis = y.ax,
    xaxis = x.ax
  )


p.daglu.3mos <- plot_ly(
  x = seq(0.05, 1, by = 0.1),
  y = colnames(daglu.3mos[[1]]),
  z = t(daglu.avg.3mos),
  type = "heatmap",
  zauto = F,
  zmin = 0,
  zmax = 10,
  # colors = colorRamp(c(rgb(0,0,0,1), rgb(0,1,0,1))),
  colors = colorRamp(c(rgb(1,1,1,1), rgb(0,1,0,1))),
  xgap = 0.4,
  ygap = 0.4
)%>%
  layout(
    title = "DAGLU-3mos",
    showlegend = T,
    yaxis = y.ax,
    xaxis = x.ax
  )

p.daglu.12mos <- plot_ly(
  x = seq(0.05, 1, by = 0.1),
  y = colnames(daglu.12mos[[1]]),
  z = t(daglu.avg.12mos),
  type = "heatmap",
  zauto = F,
  zmin = 0,
  zmax = 10,
  # colors = colorRamp(c(rgb(0,0,0,1), rgb(0,1,0,1))),
  colors = colorRamp(c(rgb(1,1,1,1), rgb(0,1,0,1))),
  xgap = 0.4,
  ygap = 0.4
)%>%
  layout(
    title = "DAGLU-12mos",
    showlegend = T,
    yaxis = y.ax,
    xaxis = x.ax
  )

p.daglu.24mos <- plot_ly(
  x = seq(0.05, 1, by = 0.1),
  y = colnames(daglu.24mos[[1]]),
  z = t(daglu.avg.24mos),
  type = "heatmap",
  zauto = F,
  zmin = 0,
  zmax = 10,
  # colors = colorRamp(c(rgb(0,0,0,1), rgb(0,1,0,1))),
  colors = colorRamp(c(rgb(1,1,1,1), rgb(0,1,0,1))),
  xgap = 0.4,
  ygap = 0.4
)%>%
  layout(
    title = "DAGLU-24mos",
    showlegend = T,
    yaxis = y.ax,
    xaxis = x.ax
  )

p.th.3mos <- plot_ly(
  x = seq(0.05, 1, by = 0.1),
  y = colnames(th.3mos[[1]]),
  z = t(th.avg.3mos),
  type = "heatmap",
  zauto = F,
  zmin = 0,
  zmax = 30,
  # colors = colorRamp(c(rgb(0,0,0,1), rgb(0,1,1,1))),
  colors = colorRamp(c(rgb(1,1,1,1), rgb(0,1,1,1))),
  xgap = 0.4,
  ygap = 0.4
)%>%
  layout(
    title = "TH-3mos",
    showlegend = T,
    yaxis = y.ax,
    xaxis = x.ax
  )

p.th.12mos <- plot_ly(
  x = seq(0.05, 1, by = 0.1),
  y = colnames(th.12mos[[1]]),
  z = t(th.avg.12mos),
  type = "heatmap",
  zauto = F,
  zmin = 0,
  zmax = 30,
  # colors = colorRamp(c(rgb(0,0,0,1), rgb(0,1,1,1))),
  colors = colorRamp(c(rgb(1,1,1,1), rgb(0,1,1,1))),
  xgap = 0.4,
  ygap = 0.4
)%>%
  layout(
    title = "TH-12mos",
    showlegend = T,
    yaxis = y.ax,
    xaxis = x.ax
  )

p.th.24mos <- plot_ly(
  x = seq(0.05, 1, by = 0.1),
  y = colnames(th.24mos[[1]]),
  z = t(th.avg.24mos),
  type = "heatmap",
  zauto = F,
  zmin = 0,
  zmax = 30,
  # colors = colorRamp(c(rgb(0,0,0,1), rgb(0,1,1,1))),
  colors = colorRamp(c(rgb(1,1,1,1), rgb(0,1,1,1))),
  xgap = 0.4,
  ygap = 0.4
)%>%
  layout(
    title = "TH-24mos",
    showlegend = T,
    yaxis = y.ax,
    xaxis = x.ax
  )

p.da.3.12 <- plot_ly(
  x = seq(0.05, 1, by = 0.1),
  y = colnames(da.3mos[[1]]),
  z = t(da.diff.3.12.pval.adj),
  type = "heatmap",
  zauto = F,
  zmin = 0,
  zmax = 1,
  colors = colorRamp(c("white", "red")),
  xgap = 0.4,
  ygap = 0.4
)%>%
  layout(
    title = "",
    showlegend = T,
    yaxis = y.ax,
    xaxis = x.ax
  )

p.da.3.24 <- plot_ly(
  x = seq(0.05, 1, by = 0.1),
  y = colnames(da.3mos[[1]]),
  z = t(da.diff.3.24.pval.adj),
  type = "heatmap",
  zauto = F,
  zmin = 0,
  zmax = 1,
  colors = colorRamp(c("white", "red")),
  xgap = 0.4,
  ygap = 0.4
)%>%
  layout(
    title = "",
    showlegend = T,
    yaxis = y.ax,
    xaxis = x.ax
  )

p.daglu.3.12 <- plot_ly(
  x = seq(0.05, 1, by = 0.1),
  y = colnames(daglu.3mos[[1]]),
  z = t(daglu.diff.3.12.pval.adj),
  type = "heatmap",
  zauto = F,
  zmin = 0,
  zmax = 1,
  colors = colorRamp(c("white", "red")),
  xgap = 0.4,
  ygap = 0.4
)%>%
  layout(
    title = "",
    showlegend = T,
    yaxis = y.ax,
    xaxis = x.ax
  )

p.daglu.3.24 <- plot_ly(
  x = seq(0.05, 1, by = 0.1),
  y = colnames(daglu.3mos[[1]]),
  z = t(daglu.diff.3.24.pval.adj),
  type = "heatmap",
  zauto = F,
  zmin = 0,
  zmax = 1,
  colors = colorRamp(c("white", "red")),
  xgap = 0.4,
  ygap = 0.4
)%>%
  layout(
    title = "",
    showlegend = T,
    yaxis = y.ax,
    xaxis = x.ax
  )

p.th.3.12 <- plot_ly(
  x = seq(0.05, 1, by = 0.1),
  y = colnames(th.3mos[[1]]),
  z = t(th.diff.3.12.pval.adj),
  type = "heatmap",
  zauto = F,
  zmin = 0,
  zmax = 1,
  colors = colorRamp(c("white", "red")),
  xgap = 0.4,
  ygap = 0.4
)%>%
  layout(
    title = "",
    showlegend = T,
    yaxis = y.ax,
    xaxis = x.ax
  )

p.th.3.24 <- plot_ly(
  x = seq(0.05, 1, by = 0.1),
  y = colnames(th.3mos[[1]]),
  z = t(th.diff.3.24.pval.adj),
  type = "heatmap",
  zauto = F,
  zmin = 0,
  zmax = 1,
  colors = colorRamp(c("white", "red")),
  xgap = 0.4,
  ygap = 0.4
)%>%
  layout(
    title = "",
    showlegend = T,
    yaxis = y.ax,
    xaxis = x.ax
  )

plots.da <- subplot(empty.plot, p.da.12mos, p.da.24mos,
                    p.da.3mos, p.da.3.12, p.da.3.24, shareX = T, shareY = T, nrows = 2)%>%
  layout(
    showlegend = F
  )

plots.daglu <- subplot(empty.plot, p.daglu.12mos, p.daglu.24mos,
                    p.daglu.3mos, p.daglu.3.12, p.daglu.3.24, shareX = T, shareY = T, nrows = 2)%>%
  layout(
    showlegend = F
  )

plots.th <- subplot(empty.plot, p.th.12mos, p.th.24mos,
                    p.th.3mos, p.th.3.12, p.th.3.24, shareX = T, shareY = T, nrows = 2)%>%
  layout(
    showlegend = F
  )

sp <- subplot(
  empty.plot, p.th.12mos, p.th.24mos,
  p.th.3mos, p.th.3.12, p.th.3.24,
  empty.plot, p.da.12mos, p.da.24mos,
  p.da.3mos, p.da.3.12, p.da.3.24,
  empty.plot, p.daglu.12mos, p.daglu.24mos,
  p.daglu.3mos, p.daglu.3.12, p.daglu.3.24,
  shareX = T, shareY = T,
  nrows = 6
)%>%
  layout(
    showlegend = F
  )
sp

setwd("~/OneDrive - Indiana University/ShaferLab/Others/Jackie/INTRSECT Ant.-Post._Medial-Lateral/Medial-Lateral/HeatmapRaw/")
orca(sp, "allHeatmaps_wStatSig.pdf", width = 768, height = 1024)