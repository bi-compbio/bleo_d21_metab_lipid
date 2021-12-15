#' Summarise an assay
#'
#' This function summarises a data matrix, as stored in a
#' summarizedExperiment object
#'
#' @param mae a \code{MultiAssayExperiment} or a
#' \code{SummarizedExperiment} object
#' @param name.se character, name of the assay to summarise if
#' \code{mae} is a \code{MultiAssayExperiment}
#' @param quote.aes call to \code{aes()}, graphical mapping for the pca plot.
#' It can map graphical parameters from the \code{colData}.
#' @param perc numeric, percentile for the histogram representation
#' @param .transform function to apply to the data matrix
#' (e.g. \code{log10}). Defaults to none.
#' @param comps two (for scatterplot) or more (pairs plot) integers,
#' indicating which components should be plotted.
#' @param return.interim logical, whether to return a list with the interim
#' results intead of just the final ggplot object
#'
#' @examples
#' data(sample.mae)
#' ## MAE
#' summarise_dataset(sample.mae, "GISTIC")
#' summarise_dataset(sample.mae, "GISTIC", comps = 2:1)
#' summarise_dataset(sample.mae, "GISTIC", comps = 1:3,
#' quote.aes = aes(colour = assay))
#' ## SE
#' sample.se <- sample.mae[["GISTIC"]]
#' colData(sample.se)$type <- c("A", "A", "B")
#' summarise_dataset(sample.se, comps = 1:2,
#' quote.aes = aes(colour = type, shape = type))
#' ## list of interim results
#' list.res <- summarise_dataset(sample.se, return.interim = TRUE)
#' head(list.res$df.plot)
#'
#' @import ggplot2
#' @importFrom GGally ggpairs
#' @import magrittr
#' @importFrom graphics hist par
#' @importFrom stats quantile
#' @importFrom plyr join
#' @importFrom dplyr mutate
#' @import MultiAssayExperiment
#' @importFrom SummarizedExperiment assay
#' @importFrom pcaMethods pca scores
#' @export
summarise_dataset <- function(
  mae, name.se, quote.aes = NULL,
  perc = 0.95, .transform = identity,
  comps = c(1, 2), return.interim = FALSE) {
  
  if (inherits(mae, "MultiAssayExperiment")) {
    message("Summarising a MultiAssayExperiment object")
    
    # transpose to get the usual rows = samples, cols = features
    mat <- do.call(.transform,
                   list(t(MultiAssayExperiment::assays(mae)[[name.se]])))
    
    # for PCA plot
    df.join <- as.data.frame(MultiAssayExperiment::sampleMap(mae))
    join.by <- "primary"
    
    df.coldata <- mae %>%
      MultiAssayExperiment::colData() %>%
      as.data.frame %>%
      dplyr::mutate(primary = rownames(.))
    
  } else if (inherits(mae, "SummarizedExperiment")) {
    message("Summarising a SummarizedExperiment object")
    
    mat <- do.call(.transform, list(t(SummarizedExperiment::assay(mae))))
    
    df.join <- data.frame(colname = rownames(mat))
    join.by <- "colname"
    
    df.coldata <- mae %>%
      SummarizedExperiment::colData() %>%
      as.data.frame %>%
      dplyr::mutate(primary = rownames(.),
                    colname = primary)
  } else {
    stop("'mae' must be of type MultiAssayExperiment or SummarizedExperiment")
  }
  
  message("Dimension")
  show(dim(mat))
  
  message("General properties")
  v <- as.vector(mat)
  print(summary(v))
  
  message("Missings")
  message("- Total")
  print(sum(is.na(v)))
  mat.na <- is.na(mat)
  sample.na <- rowSums(mat.na)
  message("- Number of samples with missings")
  print(sum(sample.na > 0))
  message("- Number of missings per sample")
  print(summary(sample.na))
  message("- Number of missings per feature")
  print(summary(colSums(mat.na)))
  
  message("Summary of the features")
  col.means <- apply(mat, 2, mean)
  col.sds <- apply(mat, 2, sd)
  message("- Means")
  print(summary(col.means))
  message("- Standard deviations")
  print(summary(col.sds))
  
  message("Histograms")
  graphics::par(mfrow = c(1, 2))
  graphics::hist(v,
                 xlab = "Numeric value",
                 main = "Whole matrix")
  
  q3 <- stats::quantile(v, perc, na.rm = TRUE)
  graphics::hist(v[v < q3],
                 xlab = "Numeric value",
                 main = paste0("Values under perc. ", perc))
  
  message("PCA")
  nPcs <- max(comps, 2)
  pca <- pcaMethods::pca(mat, method = "nipals",
                         scale = "uv", center = TRUE,
                         nPcs = nPcs)
  show(pca)
  
  name.comp <- paste0("PC", comps)
  name.compvar <- paste0(name.comp, " (", round(pca@R2[comps]*100, 2), "%)")
  
  if (inherits(mae, "MultiAssayExperiment")) {
    df.join <- as.data.frame(MultiAssayExperiment::sampleMap(mae))
  } else {
    
  }
  
  # df.join is a dummy data.frame with one single column
  # if mae is a SE object, otherwise it contains the mapping to the
  # primary sample IDs in the MAE
  df.pca <- pcaMethods::scores(pca) %>% as.data.frame %>%
    dplyr::mutate(colname = rownames(.)) %>%
    plyr::join(df.join)
  
  # check if PCx columns already exist, overwrite if so
  col.existing <- intersect(colnames(df.coldata), name.comp)
  if (length(col.existing) > 0) {
    message("Overwriting already existing PCx columns: ", col.existing)
    df.coldata <- df.coldata[setdiff(colnames(df.coldata), col.existing)]
  }
  
  df.plot <- plyr::join(df.pca, df.coldata, by = join.by, type = "left")
  
  if (length(comps) == 2) {
    # two component: scatterplot
    gg.obj <- ggplot(df.plot, aes_string(x = name.comp[1], y = name.comp[2])) +
      geom_hline(yintercept = 0, lty = 2) +
      geom_vline(xintercept = 0, lty = 2) +
      geom_point(quote.aes) +
      coord_fixed() +
      xlab(name.compvar[1]) +
      ylab(name.compvar[2])
  } else if (length(comps) > 2) {
    # more than two: pairs plot
    gg.obj <- GGally::ggpairs(
      df.plot,
      mapping = quote.aes,
      columns = name.comp,
      columnLabels = name.compvar)
  } else {
    message("ncomps must have length 2 or more to plot the PCA components")
  }
  
  if (return.interim) {
    list(gg.obj = gg.obj, df.plot = df.plot, pca = pca, names.pc = name.compvar)
  } else {
    gg.obj
  }
}




# I am just tired and went for an easy solution for the plot aspect issue
# square PCAs, fixing the scale
# see https://stackoverflow.com/questions/13445753/force-ggplot2-scatter-plot-to-be-square-shaped
# which refers to the preferred solution in 
# https://stackoverflow.com/questions/7705345/how-can-i-extract-plot-axes-ranges-for-a-ggplot2-object/35372274#35372274
# 
# a pity that coord_fixed and theme(aspect.ratio = 1) do not get along..
get_limits_square <- function(plt, prop.margin = 1.1) {
  # y-range
  rg.x <- layer_scales(plt)$y$range$range
  
  # x-range
  rg.y <- layer_scales(plt)$x$range$range
  
  # take most extreme values, add margin
  rg <- c(min(rg.x, rg.y), max(rg.x, rg.y))*prop.margin
  
  rg
}







# spread multiple columns
# https://community.rstudio.com/t/spread-with-multiple-value-columns/5378/2
myspread <- function(df, key, value) {
  # quote key
  keyq <- rlang::enquo(key)
  # break value vector into quotes
  valueq <- rlang::enquo(value)
  s <- rlang::quos(!!valueq)
  df %>% gather(variable, value, !!!s) %>%
    unite(temp, !!keyq, variable) %>%
    spread(temp, value)
}

# shrink area of venn diagram for cut areas
margin_venn <- function(prop = .9) {
  grid::pushViewport(
    grid::viewport(width = unit(prop, "npc"), 
                   height = unit(prop, "npc")))
}

# function from the RAM package, modified to allow viewports
# adds a callback option after drawing the grid
ram_group_venn <- function (
  vectors, quote.callback.grid, cat.cex = 1.5, cex = 1, 
  cat.pos = NULL, cat.dist = NULL, 
  label = TRUE, lab.cex = 1, lab.col = "black", fill = NULL, 
  file = NULL, ext = NULL, width = 8, height = 8) 
{
  save <- !is.null(file)
  if (save) {
    RAM:::.get.dev(file, ext, height = height, width = width)
  }
  if (!requireNamespace("VennDiagram")) {
    stop("package 'VennDiagram' is required for this function")
  }
  if (!requireNamespace("RColorBrewer")) {
    stop("package 'RColorBrewer' is required for this function")
  }
  if (!requireNamespace("grid")) {
    stop("package 'grid' is required to use this function")
  }
  len <- length(vectors)
  if (is.null(fill)) {
    if (len == 2) {
      fill = c("lightpink", "lightblue")
    }
    else {
      fill = RColorBrewer::brewer.pal(len, "Pastel1")
    }
  }
  else {
    if (length(fill) == len) {
      fill = fill
    }
    else if (length(fill) > len) {
      warning(paste("more colors being provided than required, will ignore ", 
                    length(fill) - len, " colors", sep = ""))
      fill = fill[1:len]
    }
    else {
      warning("not enough colors being provided, will use default")
      if (len == 2) {
        fill = c("lightpink", "lightblue")
      }
      else {
        fill = RColorBrewer::brewer.pal(len, "Pastel1")
      }
    }
  }
  if (len > 2 && label) {
    warning("currently only support 2 groups to have actual item labels; will only use numbers")
  }
  else if (len > 5 || len < 2) {
    stop("please provide 2 to 5 vectors")
  }
  alpha = rep(0.5, len)
  if (!is.null(cat.pos) && !is.null(cat.dist)) {
    v <- VennDiagram::venn.diagram(vectors, fill = fill, 
                                   alpha = alpha, cat.dist = cat.dist, cat.pos = cat.pos, 
                                   cat.fontface = "bold", cat.cex = cat.cex, cex = cex, 
                                   filename = NULL)
  }
  else if (!is.null(cat.pos) && is.null(cat.dist)) {
    v <- VennDiagram::venn.diagram(vectors, fill = fill, 
                                   alpha = alpha, cat.pos = cat.pos, cat.fontface = "bold", 
                                   cat.cex = cat.cex, cex = cex, filename = NULL)
  }
  else if (is.null(cat.pos) && !is.null(cat.dist)) {
    v <- VennDiagram::venn.diagram(vectors, fill = fill, 
                                   alpha = alpha, cat.fontface = "bold", cat.dist = cat.dist, 
                                   cat.cex = cat.cex, cex = cex, filename = NULL)
  }
  else {
    v <- VennDiagram::venn.diagram(vectors, fill = fill, 
                                   alpha = alpha, cat.fontface = "bold", cat.cex = cat.cex, 
                                   cex = cex, filename = NULL)
  }
  if (len > 2 && len <= 5) {
    grid::grid.newpage()
    eval(quote.callback.grid)
    grid::grid.draw(v)
  }
  if (len == 2) {
    if (!label) {
      grid::grid.newpage()
      eval(quote.callback.grid)
      grid::grid.draw(v)
    }
    else {
      name <- lapply(v, names)
      lapply(v, function(i) i$label)
      v.labels <- lapply(v, function(i) i$label)
      v.lab <- vector()
      for (i in 1:length(v.labels)) {
        if (length(v.labels[[i]] %in% names(vectors)) != 
            0 && isTRUE(v.labels[[i]] %in% names(vectors))) {
          v.lab <- c(v.lab, v.labels[[i]])
        }
      }
      v1 <- vectors[[v.lab[1]]]
      v2 <- vectors[[v.lab[2]]]
      v[[5]]$label <- paste(c(v[[5]]$label, setdiff(v1, 
                                                    v2)), collapse = "\n")
      v[[5]]$gp$cex <- lab.cex
      v[[5]]$gp$col <- lab.col
      v[[6]]$label <- paste(c(v[[6]]$label, setdiff(v2, 
                                                    v1)), collapse = "\n")
      v[[6]]$gp$cex <- lab.cex
      v[[6]]$gp$col <- lab.col
      v[[7]]$label <- paste(c(v[[7]]$label, intersect(v1, 
                                                      v2)), collapse = "\n")
      v[[7]]$gp$cex <- lab.cex
      v[[7]]$gp$col <- lab.col
      grid::grid.newpage()
      eval(quote.callback.grid)
      grid::grid.draw(v)
    }
  }
  if (save) {
    dev.off()
  }
  if (exists("v")) v
}

get_assay_with_coldata <- function(mae, str.assay, biomarkers) {
  sumexp <- MultiAssayExperiment::assays(mae)[[str.assay]]
  
  coldata <- MultiAssayExperiment::colData(mae) %>% as.data.frame %>% 
    tibble::rownames_to_column("primary")
  samplemap <- MultiAssayExperiment::sampleMap(mae) %>% as.data.frame %>%
    filter(assay == str.assay)
  
  df.cols <- join(coldata, samplemap, type = "left")
  
  sumexp[as.character(biomarkers), , drop = FALSE] %>% 
    t %>%  
    as.data.frame %>%
    tibble::rownames_to_column("colname") %>%
    gather("biomarker", "biomarkervalue", -colname) %>%
    join(df.cols, type = "left")
}


fit_rank_models <- function(mae, formula.mod, fun.mod = lm, ...) {
  # obtain phenotypic data
  df.coldata <- mae %>% colData %>% as.data.frame %>% 
    dplyr::mutate(primary = rownames(.)) 
  
  plyr::ldply(names(mae), function(name.se) {
    # make sure the order of the samples is preserved for the linear models
    # browser()
    se <- assays(mae)[[name.se]]
    
    col.names <- colnames(se)
    df.colnames <- dplyr::filter(
      as.data.frame(sampleMap(mae)), 
      colname %in% col.names & assay == name.se) %>% 
      join(df.coldata, by = "primary", type = "left") %>%
      set_rownames(.$colname)
    
    mat.anova <- plyr::adply(assay(se), 1, function(y) {
      df.colnames$y <- y
      
      mod <- do.call(
        fun.mod, 
        list(formula = formula.mod, data = df.colnames)
      )
      # browser()
      anova(mod) %>%
        as.data.frame %>%
        mutate(coefficient = rownames(.), assay = name.se)
    }, .id = "biomarker", .parallel = TRUE)
  }, .id = "assay", .progress = "text")
}

# Using regular lms
# 
# The data must be log transformed if you want to use coef as lfc
# this function only allows formulae that are linear on x!
# But you can add covariates
compute_pvals_lfcs <- function(se, var.fc, formula.mod = "y ~ x", 
                               fun.mod = lm, db.ids = c("entrez.id", "KEGG"),
                               ...) {
  
  df.coldata <- se %>% colData %>% as.data.frame 
  df.coldata$x <- droplevels(df.coldata[[var.fc]])
  
  df.rowdata <- se %>% rowData
  
  if (nlevels(df.coldata$x) != 2) 
    stop("Number of levels in ", var.fc, " is not exactly 2")
  
  lev.alt <- paste0("x", levels(df.coldata$x)[2])
  which.db <- intersect(db.ids, names(df.rowdata)) %>% head(1)
  
  if (length(which.db) == 0) 
    stop("the columns ", db.ids, " were not found in the rowData")
  
  df.anova <- plyr::adply(assay(se), 1, function(y) {
    df.coldata$y <- y
    
    # models
    mod <- do.call(
      fun.mod, 
      list(formula = formula.mod, data = df.coldata)
    )
    # browser()
    
    coef <- coef(mod)[lev.alt] 
    pval <- anova(mod) %>%
      as.data.frame %>%
      `[`("x", "Pr(>F)")
    
    data.frame(scoef = coef, spval = pval) %>%
      set_names(paste(names(.), var.fc, sep = "."))
  }, .id = "biomarker", .parallel = TRUE)
  df.anova[[paste0("spval.adj.", var.fc)]] <- p.adjust(df.anova$spval, method = "fdr")
  
  dplyr::mutate(df.anova, 
                db.id = as.character(df.rowdata[[which.db]]), 
                db = which.db)
}

# get residuals, test group variances
compute_statistics_lfcs <- function(se, var.fc, formula.mod = "y ~ x", 
                                   fun.mod = lm, 
                                   ...) {
  
  df.coldata <- se %>% colData %>% as.data.frame 
  df.coldata$x <- droplevels(df.coldata[[var.fc]])
  
  df.rowdata <- se %>% rowData
  
  if (nlevels(df.coldata$x) != 2) 
    stop("Number of levels in ", var.fc, " is not exactly 2")
  
  lev.alt <- paste0("x", levels(df.coldata$x)[2])
  
  df.res <- plyr::adply(assay(se), 1, function(y) {
    df.coldata$y <- y
    
    # comparing variances of abundance (y) between groups (x)
    grp.split <- split(y, df.coldata$x)
    vartest <- var.test(grp.split[[1]], grp.split[[2]], ratio = 1, alternative = "two.sided")
    vt <- vartest$p.value
    
    # models
    mod <- do.call(
      fun.mod, 
      list(formula = formula.mod, data = df.coldata)
    )
    # browser()
    
    # data.frame(residual = residuals(mod))
    resid <- residuals(mod)
    
    if (sd(resid) < 1e-10) sw <- NA
    else sw <- shapiro.test(resid)$p.value
    
    # breusch-pagan test
    bp <- lmtest::bptest(mod)$p.value
    
    # levene test
    lp <- car::leveneTest(mod)$`Pr(>F)`
    stopifnot(length(lp) == 2)
    stopifnot(is.na(lp[2]))
    
    data.frame(pval.shapiro = sw, pval.ftest = vt, pval.bptest = bp, pval.ltest = lp[1])
  }, .id = "biomarker", .parallel = TRUE)
  df.res[[paste0("fdr.shapiro.", var.fc)]] <- p.adjust(df.res$pval.shapiro, method = "fdr")
  df.res[[paste0("fdr.ftest.", var.fc)]] <- p.adjust(df.res$pval.ftest, method = "fdr")
  df.res[[paste0("fdr.bptest.", var.fc)]] <- p.adjust(df.res$pval.bptest, method = "fdr")
  df.res[[paste0("fdr.ltest.", var.fc)]] <- p.adjust(df.res$pval.ltest, method = "fdr")
  
  df.res
}

# Using limma
limma_pvals_lfcs <- function(se, var.fc, formula.design, robust = FALSE,
                             method = "ls", 
                             db.ids = c("entrez.id", "KEGG"),
                             format.like.lm = FALSE,
                             ...) {
  
  df.coldata <- se %>% colData %>% as.data.frame 
  design <- model.matrix(formula.design, data = df.coldata)
  
  df.rowdata <- se %>% rowData
  
  if (nlevels(df.coldata[[var.fc]]) != 2) 
    stop("Number of levels in ", var.fc, " is not exactly 2")
  
  # the contrast should be on the non-reference level 
  lev.alt <- paste0(var.fc, levels(df.coldata[[var.fc]])[-1])
  if (!(lev.alt %in% colnames(design)))
    stop("var.fc (", var.fc, ") leads to the contrast of ", lev.alt, 
         " which does not exist after applying the formula ", formula.design)
  
  which.db <- intersect(db.ids, names(df.rowdata)) %>% head(1)
  
  if (length(which.db) == 0) 
    stop("the columns ", db.ids, " were not found in the rowData")
  
  fit <- limma::lmFit(assay(se), design, method = method)
  
  cnt <- limma::makeContrasts(contrasts = lev.alt, levels = design)
  
  fit.cnt <- limma::contrasts.fit(fit, cnt) %>% limma::eBayes(robust = robust)
  
  df.ans <- limma::topTable(fit.cnt, number = Inf, adjust.method = "fdr", sort.by = "none") %>%
    dplyr::mutate(biomarker = rownames(df.rowdata), 
                  db.id = df.rowdata[[which.db]], 
                  db = which.db)
  
  if (format.like.lm) {
    df.ans[[paste0("scoef.", var.fc)]] <- df.ans$logFC
    df.ans[[paste0("spval.", var.fc)]] <- df.ans$P.Value
    df.ans[[paste0("spval.adj.", var.fc)]] <- df.ans$adj.P.Val
    
    df.ans <- dplyr::select(df.ans, -logFC, -AveExpr, -t, 
                            -P.Value, -adj.P.Val, -B)
  }
  
  df.ans
}

# biohelpeR
# import kegg brite file
keggbrite2df <- function(file.keg) 
{
  b <- readLines(file.keg)
  brite <- list()
  i <- 0
  while (i < length(b)) {
    i <- i + 1
    line <- b[i]
    if (substr(line, 1, 1) == "A") {
      name <- stringr::str_match(line, ">.+<")
      name <- as.character(substr(name, 2, nchar(name) - 
                                    1))
      brite[[name]] <- list()
      current <- name
    }
    else if (substr(line, 1, 1) == "B") {
      name <- substr(line, 4, nchar(line))
      brite[[current]][[name]] <- character()
      currentSub <- name
    }
    else if (substr(line, 1, 1) == "C") {
      name <- stringr::str_match(line, "\\d{5}")
      name <- paste0("map", name)
      brite[[current]][[currentSub]] <- c(brite[[current]][[currentSub]], 
                                          name)
    }
  }
  kegg.brite <- reshape2::melt(brite)
  names(kegg.brite)[1] <- "kegg.id"
  kegg.brite
}

# biohelpeR
# list of dfs to tabset in knitr
list_list_dfs_to_knitr <- function (list.list.dfs, header = "###") 
{
  txt <- NULL
  for (grp in names(list.list.dfs)) {
    txt <- paste0(txt, header, " ", grp, " \n\n")
    for (db in names(list.list.dfs[[grp]])) {
      txt <- paste0(txt, db, "\n`r list.list.dfs$'", grp, 
                    "'$'", db, "'`\n")
    }
  }
  paste(knitr::knit(text = txt), collapse = "\n")
}

gg_45 <- function () 
{
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
                                                     vjust = 1, hjust = 1))
}

gg_90 <- function () 
{
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                     vjust = 0.5, hjust = 1))
}


igraph_add_graphical_params <- function (graph = NULL, layout = FALSE, graph.layout = NULL, 
          NamesAsLabels = TRUE, 
          ...) 
{
  if (vcount(graph) == 0) {
    warning("The graph is empty and won't be plotted.")
    return(invisible())
  }
  if ("GO.simil" %in% igraph::list.vertex.attributes(graph)) {
    GO.simil <- V(graph)$GO.simil
    GO.annot <- TRUE
  }
  else {
    GO.annot <- FALSE
  }
  add.vertex.shape("triangle", clip = vertex.shapes("circle")$clip, 
                   plot = FELLA:::mytriangle)
  graph.input <- V(graph)$input
  graph.com <- as.character(V(graph)$com)
  vertex.shape <- rep("circle", vcount(graph))
  vertex.shape[graph.input] <- "square"
  vertex.number <- vcount(graph)
  graph.asp <- 1
  if (is.null(graph.layout)) 
    graph.layout <- layout.auto(graph)
  graph.layout <- layout.norm(graph.layout, xmin = -1, xmax = 1, 
                              ymin = -1, ymax = 1)
  mapSolidColor <- c(`1` = "#CD0000", `2` = "#CD96CD", `3` = "#FFA200", 
                     `4` = "#8DB6CD", `5` = "#548B54")
  vertex.color <- vapply(V(graph), function(y) {
    solidColor <- mapSolidColor[graph.com[y]]
    if (!GO.annot) 
      return(solidColor)
    GO.y <- GO.simil[y]
    if (!is.na(GO.y)) {
      if (GO.y < 0.5) 
        solidColor <- "#FFD500"
      else if (GO.y < 0.7) 
        solidColor <- "#FF5500"
      else if (GO.y < 0.9) 
        solidColor <- "#FF0000"
      else solidColor <- "#B300FF"
    }
    solidColor
  }, FUN.VALUE = character(1))
  vertex.frame.color <- rep("black", vcount(graph))
  if (GO.annot) {
    vertex.frame.color[!is.na(GO.simil)] <- "#CD0000"
    vertex.shape[!is.na(GO.simil)] <- "triangle"
  }
  mapSize <- c(`1` = 7, `2` = 5.5, `3` = 4.25, `4` = 3.5, `5` = 3)
  vertex.size <- mapSize[graph.com]
  vertex.size[graph.input] <- 4
  vertex.size <- vertex.size * (300/vcount(graph))^(1/3)
  vertex.label.dist <- 0.1 * (300/vcount(graph))^(1/3)
  vertex.label.degree <- -pi/2
  if (NamesAsLabels) {
    vertex.label <- V(graph)$label
  }
  else {
    vertex.label <- V(graph)$name
  }
  
  graph$layout = graph.layout
  graph$asp = graph.asp
  
  V(graph)$size = vertex.size
  V(graph)$label = vertex.label
  V(graph)$label.dist = vertex.label.dist
  V(graph)$label.color = vertex.color
  V(graph)$label.degree = vertex.label.degree
  V(graph)$frame.color = vertex.frame.color
  V(graph)$color = vertex.color
  V(graph)$shape = vertex.shape
  
  E(graph)$color = "#000000AA"
  E(graph)$arrow.size = 0.25
  
  graph
}

# helper to define grapical params in FELLA
igraph_to_visnetwork <- function(g, width = NULL, height = NULL) {
  
  id <- V(g)$name
  label <- V(g)$label
  nodes <- data.frame(id, label, stringsAsFactors = FALSE)
  
  # GO labels?
  if ("GO.simil" %in% igraph::list.vertex.attributes(g)) {
    GO.simil <- unlist(V(g)$GO.simil)
    GO.annot <- TRUE
  } else {
    GO.annot <- FALSE
  }
  
  map.com <- c("pathway", "module", "enzyme", "reaction", "compound")
  map.color <- c("#E6A3A3", "#E2D3E2", "#DFC1A3", "#D0E5F2", "#A4D4A4")
  map.labelcolor <- c("#CD0000", "#CD96CD", "#CE6700", 
                      "#8DB6CD", "#548B54")
  map.nodeWidth <- c(40, 30, 25, 22, 22)
  
  # input node?
  nodeShape <- ifelse(
    V(g)$input,
    "box",
    "ellipse"
  )
  nodes$group <- map.com[V(g)$com]
  nodes$color <- map.color[V(g)$com]
  
  # width
  nodes$value <- map.nodeWidth[V(g)$com]
  nodes$shape <- nodeShape
  
  # Change color and label if GO annotations are present
  if (GO.annot) {
    ids <- !is.na(GO.simil)
    GO.semsim <- GO.simil[ids]
    GO.hits <- names(GO.semsim)
    if (!is.null(GO.hits)) {
      newColor <- sapply(
        GO.semsim, 
        function(x) {
          if (x < 0.5) return("#FFD500")
          else if (x < 0.7) return("#FF5500")
          else if (x < 0.9) return("#FF0000")
          return("#B300FF")
        }
      )
      newName <- paste0(nodes$label[ids], "[", GO.hits, "]")
      newShape <- "triangle"
      
      # modify name and color
      nodes$label[ids] <- newName
      nodes$color[ids] <- newColor
      nodes$shape[ids] <- newShape
    }
  }
  
  # tooltip
  nodeLink <- paste0(
    "<a href=\"http://www.genome.jp/dbget-bin/www_bget?",
    V(g)$name, "\"", "\ target=\"_blank", "\">", V(g)$name, "</a>")
  if(vcount(g) == 0) nodeLink <- character(0)
  
  nodes$title <- nodeLink
  
  source <- V(g)[get.edgelist(g)[, 1]]$name
  target <- V(g)[get.edgelist(g)[, 2]]$name
  edges <- data.frame(source, target, stringsAsFactors = FALSE)
  
  names(edges) <- c("from", "to")
  
  net <- list(
    nodes = nodes,
    edges = edges)
  
  visNetwork::visNetwork(width = width, height = height, 
                         nodes = net$nodes, edges = net$edges)  %>%
    visNetwork::visIgraphLayout() %>%
    visNetwork::visEdges(smooth = FALSE) %>% 
    # visNetwork::visExport(type = "pdf") %>% #really bad quality
    visNetwork::visOptions(
      selectedBy = "group", 
      nodesIdSelection = TRUE,
      highlightNearest = TRUE, 
      manipulation = TRUE)
  
}
