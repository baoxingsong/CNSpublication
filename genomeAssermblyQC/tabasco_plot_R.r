suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(reshape))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

option.list = list(
  make_option(
    c("-d", "--directory"),
    type = "character",
    default = getwd(),
    help = "directory containing TABASCO *.count files, default current working directory",
    metavar = "path"
  ),
  make_option(
    c("-l", "--labels"),
    type = "character",
    default = NULL,
    help = "optional comma-separated list of axis labels, default is basename of each *.count file",
    metavar = "string"
  ),
  make_option(
    c("-t", "--title"),
    type = "character",
    default = 'TABASCO results',
    help = "optional plot title, default 'TABASCO results'",
    metavar = "string"
  ),
  make_option(
    c("-p", "--prefix"),
    type = "character",
    default = "tabasco_plot",
    help = "prefix of output file, default 'tabasco_plot.pdf'",
    metavar = "string"
  )
  # make_option(
  #   c("--width"),
  #   type = "integer",
  #   default = 250,
  #   help = "horizontal dimension of output plot (mm)",
  #   metavar = "integer"
  # ),
  # make_option(
  #   c("--height"),
  #   type = "integer",
  #   default = 100,
  #   help = "vertical dimension of output plot (mm)",
  #   metavar = "integer"
  # )
)

opt_parser = OptionParser(option_list = option.list, prog = "Rscript ./tabasco_plot.R")
option.list = parse_args(opt_parser)
count.files <-
  list.files(path = option.list$directory,
             pattern = "*.count",
             full.names = T)

if (is.null(option.list$directory)) {
  print_help(opt_parser)
  stop("A directory argument must be provided", call. = F)
} else if (length(count.files) == 0) {
  print_help(opt_parser)
  cat(count.files)
  stop("No files with *.count extension found in provided directory",
       call. = F)
}

if (length(option.list$labels) == 0) {
  count.labels <-
    count.files %>%
    sapply(function(x)
      basename(x) %>% sub(pattern = "*.count", replacement = ""),
      USE.NAMES = F)
} else {
  count.labels <- option.list$labels %>%
    strsplit(",") %>%
    unlist(use.names = F)
  if (length(count.labels) != length(count.files)) {
    print_help(opt_parser)
    cat("Files:\t", count.files, "\n")
    cat("Labels:\t", count.labels, "\n")
    stop(
      "Unequal number of labels and *.count files ",
      paste(
        "(",
        length(count.files),
        ",",
        length(count.labels),
        ")",
        sep = ""
      ),
      call. = F
    )
  }
}

#list of counts by file
count.list <- lapply(count.files, read.delim)

#dataframe for plotting
count.df <- suppressMessages(
  count.list  %>%
    melt() %>%
    `colnames<-`(c("variable", "count", "label")) %>%
    mutate(
      percent = with(., tapply(count, label, function(x)
        x / sum(x) * 100)) %>%
        unlist(use.names = F) %>%
        round(1),
      variable = factor(
        variable,
        levels = levels(variable) %>% rev,
        labels = c("missing (M)", "fragmented (F)", "duplicated (D)", "single-copy (S)")
      ),
      label = factor(label, labels = count.labels)
    )
)

#table of counts by file
count.table <-
  count.list %>% 
  bind_rows() %>% 
  `row.names<-`(count.labels)

# annotation for bars
ann <- count.table %>%
  apply(1, function(x)
    paste(c("S", "D", "F", "M"), x, sep = ":", collapse = "  ")) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  `colnames<-`(c("label", "text"))

# add number of TABASCOs to title if consistent across *counts files
# otherwise, add to annotation for each sample
if (var(rowSums(count.table))) {
  # ann$text <-
  #   ann$text %>%
  #   paste("  n:", rowSums(count.table), sep = "")
  plot.title <-
    option.list$title
} else {
  plot.title <-option.list$title
    # paste(option.list$title, " (n=", sum(count.table[1, ]), ")", sep = "")
}

#outputs to pdf due to formatting issues with png (margins were being cutoff, couldn't troubleshoot successfully)
#replacing "pdf" with "tiff" works also
output.file <- paste(option.list$prefix, "pdf", sep = ".")

# plot.colors <- viridisLite::plasma(10, direction = -1)[c(2,4,9:10)] %>% rev
# plot.colors <- c("dodgerblue4", "dodgerblue2", "firebrick", "orange")
plot.colors <- c("dodgerblue4", "dodgerblue2", "orange", "darkgoldenrod1")
# plot.colors <- viridisLite::plasma(4, end=0.9, direction = -1)

#making plot
p1 <- 
  ggplot(count.df) +
  geom_col(aes(x = label, y = percent, fill = variable), 
           width = 0.67) +
  geom_text(
    aes(x = label, y = 50, label = text), 
    # size=6,
    data = ann,
    color = "black",
    fontface = "bold"
  ) +
  scale_fill_manual(values = plot.colors) +
  # scale_fill_viridis_d(end = 0.9,
  #                      option = "plasma",
  #                      direction = -1) +
  coord_flip(clip = "off", expand = T) +
  guides(fill = guide_legend(reverse = TRUE)) +
  ggtitle(label = plot.title) +
  labs(y = "% TABASCOs") +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = 'bold', size = rel(1.3)),
    legend.title = element_blank(),
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(margin = margin(r = 2, unit = "mm")),
    legend.box.spacing = unit(0, units = "mm"),
    plot.title = element_text(
      hjust = 0.5,
      size = rel(2),
      face = "bold"
    ), 
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(face = 'bold'),
    panel.grid = element_blank()
  )

save_plot(output.file,
          p1,
          base_height = 6,
          base_width = 8)