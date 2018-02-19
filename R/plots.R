#' @S3method plot statisRes
#' @import ggplot2
#' @import ggrepel

plot.interstructure <- function(x, table.names = NULL, ...) {
  y <- NULL
  .e <- environment()
  #names <-match.arg(names)
  df <- data.frame("y" = rep(0, nrow(x$vectors)),
                   "x" = sqrt(x$values)[1] * x$vectors[,1])
  if (!is.null(table.names)){
    rownames(df) <- table.names
  }else{
    table.names <- as.character(c(1:nrow(df)))
  }
  p <- ggplot(df, aes(x = x, y = y), environment = .e) +
    geom_hline(aes(yintercept = 0), colour = "gray65") +
    geom_point(size = 3, aes(shape = as.factor(c(1:nrow(df))))) +
    expand_limits(x = c(0,1), y = c(-0.1, 0.1)) +
    scale_x_continuous(paste("Comp1 -- Quality =",
                             round(x$values[1] * 100 / sum(x$values), 2), "%",
                             sep = " ")) +
    scale_shape_discrete(name = "Tables",
                         breaks = as.factor(c(1:length(table.names))),
                         labels = table.names) +
    #theme_bw() +
    theme(plot.title = element_text(lineheight = .8, face = "bold"),
          axis.line = element_line(colour = "black"),
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    scale_y_continuous(breaks = NULL)

  print(p)
}



plot.statisRes <- function(x, comp = c(1, 2), view = c("ind", "var"),
                           showing = c("compro", "time-spec"),
                           table.names = NULL, title = NULL, ...) {
  y <- V1 <- V2 <- V3 <- shapes <- xfrom <- yfrom <- xto <- yto <- NULL
  view <- match.arg(view)
  showing <- match.arg(showing)
  if (is.null(title)){
    title <- paste0(showing, " positions of ", view, "\n DIM ",
                 comp[1], " - ", comp[2])
  }else{
    title <- title
  }
  if (is.null(table.names)) {
    table.names <- c(paste0("Table", 1:length(x$observations$tcoord)))
  }
  e. <- environment()
  p <- ggplot(environment = e.) +
    scale_x_continuous(name = paste0("Comp ", comp[1])) +
    scale_y_continuous(name = paste0("Comp ", comp[2])) +
    theme_bw() +
    theme(plot.title = element_text(lineheight = .8, face = "bold"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    geom_hline(aes(yintercept = 0), colour = "gray65") +
    geom_vline(aes(xintercept = 0), colour = "gray65")
  if (view == "ind") {
    if (showing == "compro") {
      if (x$parameters$method == "sir") {
        coord <- data.frame("V1" = x$observations$coord[, comp[1]],
                            "V2" = x$observations$coord[, comp[2]],
                            "V3" = rank(rownames(x$observations$coord))
                            /nrow(x$observations$coord))
        p <- p +
          ggtitle(title) +
          geom_point(data = coord, aes(x = V1, y = V2, colour = V3), size = 5) +
          scale_color_continuous(name = "slices",
                                 breaks = coord$V3,
                                 labels = rownames(coord))
      } else {
        coord <- as.data.frame(x$observations$coord[, c(comp[1], comp[2])])
        p <- p +
          ggtitle(title) +
          geom_point(data = coord, aes(x = V1, y = V2),
                     colour = "grey", size = 3) +
          geom_text(data = coord, aes(x = V1, y = V2),
                    label = rownames(coord), size = 3)
      }

    }
    if (showing == "time-spec") {
      if (x$parameters$method == "sir") {
        for (i in 1:length(x$observations$tcoord)) {
          coord <-
            data.frame("V1" = x$observations$tcoord[[i]][, comp[1]],
                       "V2" = x$observations$tcoord[[i]][, comp[2]],
                       "cols" = rank(rownames(x$observations$tcoord[[i]]))
                       /nrow(x$observations$tcoord[[i]]),
                       "shapes" = as.factor(rep(i,
                                                nrow(x$observations$tcoord[[i]]))))
          posi <- x$observations$tcoord
          p <- p +
            geom_point(data = coord,
                       aes(x = V1, y = V2, colour = cols, shape = shapes),
                       size = 5)
          if (i != 1){
            traj <- data.frame("xfrom" = posi[[i-1]][,comp[1]],
                               "xto" = posi[[i]][,comp[1]],
                               "yfrom" = posi[[i-1]][,comp[2]],
                               "yto" = posi[[i]][,comp[2]])
            p <- p +
              geom_segment(data = traj,
                           aes(x = xfrom, y = yfrom,
                               xend = xto, yend = yto),
                           colour= "grey")
          }
        }
        p <- p +
          ggtitle(title) +
          scale_color_continuous(name = "slices",
                                 breaks = rank(rownames(x$observations$coord))
                                 /nrow(x$observations$coord),
                                 labels = rownames(x$observations$coord)) +
          scale_shape_discrete(name = "Tables",
                               breaks = as.factor(c(1:length(posi))),
                               labels = table.names)
      } else {
        ncols <- length(x$observations$tcoord)
        hues <- seq(15, 375, length = ncols + 1)
        cols <- hcl(h = hues, l = 65, c = 100)[1:ncols]
        for (i in 1:length(x$observations$tcoord)) {
          coord <- as.data.frame(x$observations$tcoord[[i]]
                                 [, c(comp[1], comp[2])])
          posi <- x$observations$tcoord
          coord$cols<-cols[i]
          p <- p +
            geom_point(data = coord, aes(x = V1, y = V2, colour = cols),
                       size = 2, alpha = 0.5) +
            geom_text(data = coord, aes(x = V1, y = V2),
                      label = rownames(coord), size = 3.5)

          if (i != 1){
            traj <- data.frame("xfrom" = posi[[i-1]][, comp[1]],
                               "xto" = posi[[i]][, comp[1]],
                               "yfrom" = posi[[i-1]][, comp[2]],
                               "yto" = posi[[i]][, comp[2]])
            p <- p +
              geom_segment(data = traj, aes(x = xfrom, y = yfrom,
                                            xend = xto, yend = yto),
                           colour= "grey")

          }
        }
        p <- p +
          ggtitle(title) +
          scale_colour_discrete(name = "Tables",
                                breaks = cols,
                                labels = table.names)
      }
    }
  } else {
    tt <- seq(0, 2 * pi, length = 100)
    circle <- data.frame(x = cos(tt), y = sin(tt))
    p <- p +
      geom_path(data = circle, aes(x = x, y = y),
                colour = "gray65") +
      geom_hline(aes(yintercept = 0), colour = "gray65") +
      geom_vline(aes(xintercept = 0), colour = "gray65")
    if (showing == "compro"){
      coord <- data.frame("V1" = x$variables$cos2[, comp[1]],
                          "V2" = x$variables$cos2[, comp[2]])
      p <- p +
        geom_point(data = coord, aes(x = V1, y = V2),
                   colour = "grey", size = 3) +
        geom_text(data = coord, aes(x = V1, y = V2),
                  label = rownames(coord), size = 3)+
        ggtitle(title)
    }
    if (showing == "time-spec"){
      ncols <- length(x$observations$tcoord)
      hues <- seq(15, 375, length = ncols + 1)
      cols <- hcl(h = hues, l = 65, c = 100)[1:ncols]
      p <- p +
        ggtitle(title) +
        scale_colour_discrete(name = "Tables",
                              breaks = cols,
                              labels = table.names)
      for (i in 1:length(x$variables$tcos2)) {
        coord <- data.frame("V1" = x$variables$tcos2[[i]][, comp[1]],
                            "V2" = x$variables$tcos2[[i]][, comp[2]])
        posi <- x$variables$tcos2
        coord$cols <- cols[i]
        if (i != 1){
          traj <- data.frame("xfrom" = posi[[i-1]][, comp[1]],
                             "xto" = posi[[i]][, comp[1]],
                             "yfrom" = posi[[i-1]][, comp[2]],
                             "yto" = posi[[i]][, comp[2]])
          p <- p +
            geom_segment(data = traj, aes(x = xfrom, y = yfrom,
                                          xend = xto, yend = yto),
                         colour = "gray")
        }
      }
#         coord <- data.frame("V1" = x$variables$tcos2[[i]][, comp[1]],
#                             "V2" = x$variables$tcos2[[i]][, comp[2]])
      coord <- as.data.frame(do.call("rbind", lapply(x$variables$tcos2,
                                       function(alist) return(alist[, comp]))))
      colnames(coord) <- c("V1","V2")
        coord$cols <- rep(cols,
                          each = nrow(x$variables$tcos2[[1]]))
        coord$variable <- rep(rownames(x$variables$tcos2[[1]]),
                              length(x$variables$tcos2))
        p <- p +
          geom_point(data = coord, aes(x = V1, y = V2, colour = cols),
                     size = 2)

          p <- p + geom_text_repel(data = coord, aes(x = V1, y = V2),
                    label = rownames(coord), size = 3.5,
                    segment.size=0, force=0.5, segment.color="white")

    }
  }
  print(p)
}
