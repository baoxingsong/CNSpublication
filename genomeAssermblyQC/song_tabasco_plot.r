suppressMessages(library(grid))
suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(reshape))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

option.list = list(
  # make_option(
  #   c("-d", "--directory"),
  #   type = "character",
  #   default = getwd(),
  #   help = "directory containing TABASCO *.count files, default current working directory",
  #   metavar = "path"
  # ),
  # make_option(
  #   c("-l", "--labels"),
  #   type = "character",
  #   default = NULL,
  #   help = "optional comma-separated list of axis labels, default is basename of each *.count file",
  #   metavar = "string"
  # ),
  make_option(
    c("-t", "--title"),
    type = "character",
    default = 'TABASCO Assessment Results',
    help = "optional plot title, default 'TABASCO results'",
    metavar = "string"
  )
  # make_option(
  #   c("-p", "--prefix"),
  #   type = "character",
  #   default = "tabasco_plot",
  #   help = "prefix of output file, default 'tabasco_plot.pdf'",
  #   metavar = "string"
  # )
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
my_title <-  option.list$title


# this code is modified from BUSCO
# !!! CONFIGURE YOUR PLOT HERE !!! 
# Output
my_output <- "song_tabasco_figure"
my_width <- 7.875
my_height <- 5.90625
my_unit <- "in"

# Colors
my_colors <- c(  "darkgoldenrod1", "orange", "dodgerblue2", "dodgerblue4")
# Bar height ratio
my_bar_height <- 0.75

# Legend

# Font
my_family <- "sans"
my_size_ratio <- 1

# !!! SEE YOUR DATA HERE !!! 

files=list.files(path = "./", pattern = "\\.count$")
if( length(files) < 1 ){
  print("there is not result file under this folder")
  quit(status=1)
}

my_species = rep(files, each = 4)
my_species = gsub("\\.count$", "", my_species)
my_species <- factor(my_species)
my_species <- factor(my_species,levels(my_species)[c(length(levels(my_species)):1)]) # reorder your species here just by changing the values in the vector :

data = read.table(files[1], header=TRUE)
v = c(data[1,1], data[1,2], data[1,3], data[1,4])
my_values=v
my_percentage = v/(sum(v))*100
if( length(files) > 1 ){
  for( i in 2:length(files) ){
    data = read.table(files[i], header=TRUE)
    v = c(data[1,1], data[1,2], data[1,3], data[1,4])
    my_values=c(my_values, v)
    my_percentage = c(my_percentage, v/(sum(v))*100)
  }
}

######################################
######################################
######################################
# Code to produce the graph
labsize = 1
if (length(levels(my_species)) > 10){
  labsize = 0.66
}
print("Plotting the figure ...")
category <- c(rep(c("S","D","F","M"),c(1)))
category <-factor(category)
category = factor(category,levels(category)[c(4,1,2,3)])
df = data.frame(my_species,my_percentage,my_values,category)

figure <- ggplot() + 
  
  geom_bar(aes(y = my_percentage, x = my_species, fill = category), data = df, stat="identity", width=my_bar_height) + 
  coord_flip() + 
  theme_gray(base_size = 8) + 
  scale_y_continuous(labels = c("0","20","40","60","80","100"), breaks = c(0,20,40,60,80,100)) + 
  scale_fill_manual(values = my_colors,labels =c(" Complete (C) and single-copy (S)  ",
                                                 " Complete (C) and duplicated (D)",
                                                 " Fragmented (F)  ",
                                                 " Missing (M)")) +   
  ggtitle(my_title) + 
  xlab("") + 
  ylab("\n%TABASCOs") + 
  
  theme(plot.title = element_text(family=my_family, colour = "black", size = rel(2.2)*my_size_ratio, face = "bold")) + 
  theme(legend.position="top",legend.title = element_blank()) + 
  theme(legend.text = element_text(family=my_family, size = rel(1.2)*my_size_ratio)) + 
  theme(panel.background = element_rect(color="#FFFFFF", fill="white")) + 
  theme(panel.grid.minor = element_blank()) + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.y = element_text(family=my_family, colour = "black", size = rel(1.66)*my_size_ratio)) +
  theme(axis.text.x = element_text(family=my_family, colour = "black", size = rel(1.66)*my_size_ratio)) + 
  theme(axis.line = element_blank()) +
  theme(axis.ticks.length = unit(.85, "cm")) + 
  theme(axis.ticks.y = element_line(colour="white", size = 0)) + 
  theme(axis.ticks.x = element_line(colour="#222222")) + 
  theme(axis.ticks.length = unit(0.4, "cm")) + 
  theme(axis.title.x = element_text(family=my_family, size=rel(1.2)*my_size_ratio)) + 
  
  guides(fill = guide_legend(override.aes = list(colour = NULL))) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

for(i in rev(c(1:length(levels(my_species))))){
  detailed_values <- my_values[my_species==my_species[my_species==levels(my_species)[i]]]
  total_tabascos <- sum(detailed_values)
  figure <- figure + 
    annotate("text", label=paste("M:", detailed_values[4], ", F:", detailed_values[3], ", C:", detailed_values[1] + detailed_values[2], " [D:", detailed_values[2], ", S:", detailed_values[1], "], n:",
                                 total_tabascos, sep=""), 
             y=3, x = i, size = labsize*4*my_size_ratio, colour = "black", hjust=0, family=my_family)
}
png(file=paste(my_output, ".png", sep=""), width = my_width, height = my_height, unit = my_unit, res=350)
figure
dev.off()
pdf(file=paste(my_output, ".pdf", sep=""), width = my_width, height = my_height)
figure
dev.off()
print("Done")
