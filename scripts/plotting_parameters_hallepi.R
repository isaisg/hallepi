#Palettes
#Palette 
#Define palette for different phosphate in soils
#palette_pi_soil <- paletteer_d(package = "awtools",palette = "mpalette",direction = -1)[1:4]
palette_pi_soil <- c("#4bb8ff","#87ff1d","#ff4c92","#0d3877")
names(palette_pi_soil) <- c("low","medium","high","low+Pi")

palette_variance <- paletteer_d(package = "dutchmasters",palette = "pearl_earring",11)[c(1:4,7)] %>%
  c(.,"white")
names(palette_variance) <- c("Fraction","Genotype","Phosphate",
                             "Fraction:Phosphate","Fraction:Genotype","Residual")

#Load palette
palette_phenotype <- c("#80A0C7",
                       paletteer_d(package = "dutchmasters",palette = "milkmaid",n = 13)[c(3,5,13)],
                       "white")
names(palette_phenotype) <- c("Treatment","Bacteria","Bacteria:Treatment",
                              "BioRep","Residual")


### Optimized for the oh.ggsave function
size_median <- 6
median_color <- "black"
size_axis_text.x <- 25
size_axis_text.y <- 45
size_axis_title.x <- 0
size_axis_title.y <- 55
size_legend_text <- 45
strip_text_size <- 55
size_title_text <- 55
legend_proportion_size <- 4
size_vline <- 0.8
color_vline <- "#D9D9D9"
size_sina <- 7
alpha_sina <- 0.5
size_anova <- 15
size_point <- 8



theme_hallepi_boxplot <- theme_ohchibi(size_axis_text.x = size_axis_text.x,size_axis_text.y = size_axis_text.y,
              size_axis_title.x = size_axis_title.x,size_axis_title.y = size_axis_title.y,
              size_legend_text = size_legend_text,
              size_title_text = size_title_text, 
              legend_proportion_size = legend_proportion_size,font_family = "Arial"
              )
