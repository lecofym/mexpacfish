# Open libraries and clean console ####
pacman:: p_load(tidyverse, # Data wrangling
                ggplot2, # Graphs
                gridExtra, # Multiple plots
                mFD, elbow, # Functional diversity
                vegan, # Turn to binary data
                geometry, # Functional volume estimation
                svglite,
                ggsci,) # Color palettes

rm(list = ls())
shell('cls')

# Load database of species presence ####
origin <- read.csv('Presence.csv',
                   header = T, stringsAsFactors = T, fileEncoding = 'latin1')

# Visualize content
str(origin)

# Functional diversity analysis ####
# Subset presence data
presence<- origin[, c(4:1238)]
rownames(presence)<- origin[, 3]


# Load database of biological traits ####
traits<- read.csv('Fish_Traits.csv',
                  header = T, stringsAsFactors = T,
                  row.names = 1)
str(traits)
sub.traits<- traits # Not for now

# Set ordinal traits as "ordinal" data in R language
traits$Size<- factor(traits$Size,
                     levels = c('1', '2', '3', '4', '5', '6'),
                     ordered = T)
traits$Mobility<- factor(traits$Mobility,
                         levels = c('1', '2', '3', '4'),
                         ordered = T)
traits$Gregariousness<- factor(traits$Gregariousness,
                               levels = c('1', '2', '3', '4'),
                               ordered = T)
traits$Position<- factor(traits$Position,
                         levels = c('1', '2', '3'),
                         ordered = T)
as_tibble(traits)

# Load database with traits information ####
types<- read.csv('Traits_info.csv',
                 header = T, stringsAsFactors = T)
sp.tr.summary(tr_cat = types,
              sp_tr = traits)


# Gower distance estimation ####
gower<- funct.dist(sp_tr = traits,
                   tr_cat = types,
                   metric = 'gower')

# Functional space quality estimation ####
qual<- quality.fspaces(sp_dist = gower,
                       maxdim_pcoa = 16, # Maximum
                       deviation_weighting = 'absolute')

# Elbow method
var.expl<- qual$details_fspaces$pc_eigenvalues$Eigenvalues/
  sum(qual$details_fspaces$pc_eigenvalues$Eigenvalues)

round(sum(var.expl[1:5]), 2) # Quality using five axes
round(sum(var.expl[1:4]), 2) # Quality using four axes

inflection_points<- data.frame(
  'Dimensions' = 1:length(var.expl),
  'Explained variance' = cumsum(var.expl))
perf_PM<- elbow(
  data = inflection_points[, c('Dimensions',
                               'Explained.variance')])


# Functional indices estimation ####
# Species coordinates
coord<- qual$details_fspaces$sp_pc_coord
funct.space.plot(sp_faxes_coord = coord[, c('PC1', 'PC2',
                                                  'PC3', 'PC4',
                                                  'PC5')],
                 faxes = NULL, faxes_nm = NULL)


# Convex hull ####
# Species coordinates
fd.coord<- coord[, c('PC1', 'PC2', 'PC3', 'PC4', 'PC5')]

# Subset site and state information
transects<- origin[, c(1:2)]
rownames(transects)<- origin[, 3]
provinces<- levels(origin[, 1])

presence.conditions <- lapply(provinces, function(x) {
  
  trans <- rownames(transects[transects$Province == x,])
  
  colSums(presence[rownames(presence) %in% trans, ])
  
})
presence.conditions <- do.call(rbind, presence.conditions)
rownames(presence.conditions) <- provinces
presence.conditions<- t(presence.conditions)


# Regional convex hull ####
# Merging coords and presence data
dim(fd.coord); dim(presence.conditions)

regional.space<- decostand(presence.conditions[, ],
                                   method = 'pa') %>%
  data.frame(fd.coord, .)

# Data per region
provinces
calif.peru <- regional.space %>%
  filter(Californian == '1') # 383 species
cortez <- regional.space %>%
  filter(Cortez == '1') # 1033 species
oceanic <- regional.space %>%
  filter(Oceanic == '1') # 441 species
panamic <- regional.space %>%
  filter(Panamic == '1') # 519 species
backbone<- regional.space %>%
  filter(Californian == '1' & Cortez == '1' &
           Oceanic == '1' & Panamic == '1')


# Functional volume ####
# Global
chg<- convhulln(coord[, 1:5], options = "FA")$vol

# Regional
# California-Peruvian province
chCaliPeru<- convhulln(calif.peru[ , c(1,2, 3, 4, 5)],
                 options = "FA")$vol
chCaliPeru<- chCaliPeru/chg*100 # Percentage

# Cortez province
chCor<- convhulln(cortez[ , c(1,2, 3, 4, 5)],
                            options = "FA")$vol
chCor<- chCor/chg*100 # Percentage

# Ocean Islands province
chOce<- convhulln(oceanic[ , c(1,2, 3, 4, 5)],
                            options = "FA")$vol
chOce<- chOce/chg*100 # Percentage

# Panamic province
chPan<- convhulln(panamic[ , c(1,2, 3, 4, 5)],
                  options = "FA")$vol
chPan<- chPan/chg*100 # Percentage


# Backbone
chBack<- convhulln(backbone[ , c(1,2, 3, 4, 5)],
                  options = "FA")$vol
chBack<- chBack/chg*100 # Percentage


# Convex hull based on the axes 1 and 2 ####
hull12.fn <- function(x)x[chull(x$PC1,x$PC2),]
regional.hull <- hull12.fn(regional.space)
CaliPeru.hull <- hull12.fn(calif.peru)
Cortez.hull <- hull12.fn(cortez)
Ocean.hull <- hull12.fn(oceanic)
Panam.hull <- hull12.fn(panamic)
Back.hull<- hull12.fn(backbone)

# Colors for each province and ultimately the backbone
province.cols<- pal_d3('category10')(5)
province.cols<-c("#D62728FF","#1F77B4FF","#FF7F0EFF","#2CA02CFF","#9467BDFF")
# Graphics ####
# California-Peruvian province
PC12_CalPeru <- ggplot(regional.space, aes (x = PC1, y = PC2))+
  geom_polygon(data = regional.hull, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0,
               col = "black")+
  geom_point(data = regional.space, aes (x = PC1, y = PC2),
             size = 0.5, shape = 3, colour = "grey")+
  geom_polygon(data = CaliPeru.hull, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0.3,
               col = province.cols[1], fill = province.cols[1])+
  labs(x = '', y = 'PCoA 2 (17%)')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(size = 25))

# Cortez province
PC12_Cor <- ggplot(regional.space, aes (x = PC1, y = PC2))+
  geom_polygon(data = regional.hull, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0,
               col = "black")+
  geom_point(data = regional.space, aes (x = PC1, y = PC2),
             size = 0.5, shape = 3, colour = "grey")+
  geom_polygon(data = Cortez.hull, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0.3,
               col = province.cols[2], fill = province.cols[2])+
  labs(x = '', y = '')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(size = 25))

# Ocean Islands province
PC12_Oce <- ggplot(regional.space, aes (x = PC1, y = PC2))+
  geom_polygon(data = regional.hull, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0,
               col = "black")+
  geom_point(data = regional.space, aes (x = PC1, y = PC2),
             size = 0.5, shape = 3, colour = "grey")+
  geom_polygon(data = Ocean.hull, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0.3,
               col = province.cols[3], fill = province.cols[3])+
  labs(x = 'PCoA 1 (29%)', y = '')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(size = 25))

# Panamic province
PC12_Pan <- ggplot(regional.space, aes (x = PC1, y = PC2))+
  geom_polygon(data = regional.hull, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0,
               col = "black")+
  geom_point(data = regional.space, aes (x = PC1, y = PC2),
             size = 0.5, shape = 3, colour = "grey")+
  geom_polygon(data = Panam.hull, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0.3,
               col = province.cols[4], fill = province.cols[4])+
  labs(x = '', y = '')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(size = 25))

# Backbone
PC12_Back <- ggplot(regional.space, aes (x = PC1, y = PC2))+
  geom_polygon(data = regional.hull, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0,
               col = "black")+
  geom_point(data = regional.space, aes (x = PC1, y = PC2),
             size = 0.5, shape = 3, colour = "grey")+
  geom_polygon(data = Back.hull, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0.3,
               col = province.cols[5], fill = province.cols[5])+
  labs(x = '', y = '')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(size = 25))


# Convex hull based on the axes 3 and 4 ####
hull34.fn <- function(x)x[chull(x$PC3,x$PC4),]
regional.hull <- hull34.fn(regional.space)
CaliPeru.hull <- hull34.fn(calif.peru)
Cortez.hull <- hull34.fn(cortez)
Ocean.hull <- hull34.fn(oceanic)
Panam.hull <- hull34.fn(panamic)


# Graphics ####
# California-Peruvian province
PC34_CalPeru <- ggplot(regional.space, aes (x = PC3, y = PC4))+
  geom_polygon(data = regional.hull, aes (x = PC3, y = PC4),
               linetype = 1, linewidth = 1, alpha = 0,
               col = "black")+
  geom_point(data = regional.space, aes (x = PC3, y = PC4),
             size = 0.5, shape = 3, colour = "grey")+
  geom_polygon(data = Cortez.hull, aes (x = PC3, y = PC4),
               linetype = 1, linewidth = 1, alpha = 0.3,
               col = province.cols[1], fill = province.cols[1])+
  labs(x = '', y = 'PCoA 4 (10%)')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(size = 25))

# Cortez province
PC34_Cor <- ggplot(regional.space, aes (x = PC3, y = PC4))+
  geom_polygon(data = regional.hull, aes (x = PC3, y = PC4),
               linetype = 1, linewidth = 1, alpha = 0,
               col = "black")+
  geom_point(data = regional.space, aes (x = PC3, y = PC4),
             size = 0.5, shape = 3, colour = "grey")+
  geom_polygon(data = Cortez.hull, aes (x = PC3, y = PC4),
               linetype = 1, linewidth = 1, alpha = 0.3,
               col = province.cols[2], fill = province.cols[2])+
  labs(x = '', y = '')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(size = 25))

# Ocean Islands province
PC34_Oce <- ggplot(regional.space, aes (x = PC3, y = PC4))+
  geom_polygon(data = regional.hull, aes (x = PC3, y = PC4),
               linetype = 1, linewidth = 1, alpha = 0,
               col = "black")+
  geom_point(data = regional.space, aes (x = PC3, y = PC4),
             size = 0.5, shape = 3, colour = "grey")+
  geom_polygon(data = Ocean.hull, aes (x = PC3, y = PC4),
               linetype = 1, linewidth = 1, alpha = 0.3,
               col = province.cols[3], fill = province.cols[3])+
  labs(x = 'PCoA 3 (13%)', y = '')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(size = 25))

# Panamic province
PC34_Pan <- ggplot(regional.space, aes (x = PC3, y = PC4))+
  geom_polygon(data = regional.hull, aes (x = PC3, y = PC4),
               linetype = 1, linewidth = 1, alpha = 0,
               col = "black")+
  geom_point(data = regional.space, aes (x = PC3, y = PC4),
             size = 0.5, shape = 3, colour = "grey")+
  geom_polygon(data = Panam.hull, aes (x = PC3, y = PC4),
               linetype = 1, linewidth = 1, alpha = 0.3,
               col = province.cols[4], fill = province.cols[4])+
  labs(x = '', y = '')+
  theme_bw()+ 
  theme(panel.grid = element_blank(),
        text = element_text(size = 25))

# Backbone
PC34_Back <- ggplot(regional.space, aes (x = PC3, y = PC4))+
  geom_polygon(data = regional.hull, aes (x = PC3, y = PC4),
               linetype = 1, linewidth = 1, alpha = 0,
               col = "black")+
  geom_point(data = regional.space, aes (x = PC3, y = PC4),
             size = 0.5, shape = 3, colour = "grey")+
  geom_polygon(data = Back.hull, aes (x = PC3, y = PC4),
               linetype = 1, linewidth = 1, alpha = 0.3,
               col = province.cols[5], fill = province.cols[5])+
  labs(x = '', y = '')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(size = 25))


# Trait convex hull ####
# Merging trait and presence data
dim(fd.coord); dim(traits)
trait.space<- data.frame(fd.coord, traits)

# Data per traits
# Size
size.1 <- subset(trait.space, trait.space$Size == "1")
size.2 <- subset(trait.space, trait.space$Size == "2")
size.3 <- subset(trait.space, trait.space$Size == "3")
size.4 <- subset(trait.space, trait.space$Size == "4")
size.5 <- subset(trait.space, trait.space$Size == "5")
size.6 <- subset(trait.space, trait.space$Size == "6")


# Convex hull based on the axes 1 and 2
regional.hull <- hull12.fn(trait.space)
hulls_S1 <- hull12.fn(size.1)
hulls_S2 <- hull12.fn(size.2)
hulls_S3 <- hull12.fn(size.3)
hulls_S4 <- hull12.fn(size.4)
hulls_S5 <- hull12.fn(size.5)
hulls_S6 <- hull12.fn(size.6)


# Graphics ####
# Size 1
PC12_S1 <- ggplot(trait.space, aes (x = PC1, y = PC2))+
  geom_polygon(data = regional.hull, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0,
               col = "black")+
  geom_point(data = trait.space, aes (x = PC1, y = PC2),
             size = 0.5, shape = 3, colour = "grey")+
  geom_polygon(data = hulls_S1, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0.3, col = "green", fill = "green")+
  theme_bw()+ theme(panel.grid = element_blank())

# Size 2
PC12_S2 <- ggplot(trait.space, aes (x = PC1, y = PC2))+
  geom_polygon(data = regional.hull, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0,
               col = "black")+
  geom_point(data = trait.space, aes (x = PC1, y = PC2),
             size = 0.5, shape = 3, colour = "grey")+
  geom_polygon(data = hulls_S2, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0.3,
               col = "blue", fill = "blue")+
  theme_bw()+ theme(panel.grid = element_blank())

# Size 3
PC12_S3 <- ggplot(trait.space, aes (x = PC1, y = PC2))+
  geom_polygon(data = regional.hull, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0,
               col = "black")+
  geom_point(data = trait.space, aes (x = PC1, y = PC2),
             size = 0.5, shape = 3, colour = "grey")+
  geom_polygon(data = hulls_S3, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0.3,
               col = "yellow", fill = "yellow")+
  theme_bw()+ theme(panel.grid = element_blank())

# Size 4
PC12_S4 <- ggplot(trait.space, aes (x = PC1, y = PC2))+
  geom_polygon(data = regional.hull, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0,
               col = "black")+
  geom_point(data = trait.space, aes (x = PC1, y = PC2),
             size = 0.5, shape = 3, colour = "grey")+
  geom_polygon(data = hulls_S4, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0.3,
               col = "purple", fill = "purple")+
  theme_bw()+ theme(panel.grid = element_blank())

# Size 5
PC12_S5 <- ggplot(trait.space, aes (x = PC1, y = PC2))+
  geom_polygon(data = regional.hull, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0,
               col = "black")+
  geom_point(data = trait.space, aes (x = PC1, y = PC2),
             size = 0.5, shape = 3, colour = "grey")+
  geom_polygon(data = hulls_S5, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0.3,
               col = "red", fill = "red")+
  theme_bw()+ theme(panel.grid = element_blank())

# Size 6
PC12_S6 <- ggplot(trait.space, aes (x = PC1, y = PC2))+
  geom_polygon(data = regional.hull, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0,
               col = "black")+
  geom_point(data = trait.space, aes (x = PC1, y = PC2),
             size = 0.5, shape = 3, colour = "grey")+
  geom_polygon(data = hulls_S6, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0.3,
               col = "orange", fill = "orange")+
  theme_bw()+ theme(panel.grid = element_blank())

# All in one graph
grid.arrange(PC12_S1, PC12_S2, PC12_S3, PC12_S4,
             PC12_S5, PC12_S6, nrow = 2)


# We remove all objects that are no longer useful ####
rm(regional.hull, Cortez.hull, CaliPeru.hull,Ocean.hull,
   Panam.hull, Back.hull, size.1, size.2, size.3, size.4, size.5,
   size.6, hulls_S1, hulls_S2, hulls_S3, hulls_S4, hulls_S5,
   hulls_S6, PC12_S1, PC12_S2, PC12_S3, PC12_S4, PC12_S5, PC12_S6,
   trait.space, qual, perf_PM, gower)


# Dominant traits per province ####
traits.cali.peru<- sub.traits[(rownames(sub.traits) %in%
                          rownames(calif.peru)), ] %>%
  add_column(Province = 'Californian')
traits.cort<- sub.traits[(rownames(sub.traits) %in%
                           rownames(cortez)), ] %>%
  add_column(Province = 'Cortez')
traits.oce<- sub.traits[(rownames(sub.traits) %in%
                          rownames(oceanic)), ] %>%
  add_column(Province = 'Oceanic')
traits.panam<- sub.traits[(rownames(sub.traits) %in%
                           rownames(panamic)), ] %>%
  add_column(Province = 'Panamic')
traits.backbone<- sub.traits[(rownames(sub.traits) %in%
                                rownames(backbone)), ] %>%
  add_column(Province = 'Backbone')
traits.province<- rbind(traits.cali.peru, traits.cort, traits.oce,
                        traits.panam, traits.backbone)
traits.province$Province <- factor(traits.province$Province,
                     levels = c('Californian', 'Cortez', 'Oceanic',
                                'Panamic', 'Backbone'),
                     ordered = T)


# Trait frequency distribution ####
# Composition per trait
#traits.province %>% group_by(trait) %>% count()

# Size
Size<- ggplot(data = traits.province,
              aes(x = Size, fill = Province)) +
  geom_histogram(aes(y = ifelse(after_stat(count) > 0,
                                after_stat(count / sum(count))*100,
                                NA)),
    binwidth = 0.5, position = position_dodge(0.7), color = 'black', alpha=0.3) +
  labs(x = 'Size', y = 'Relative frequency (%)')+
  scale_y_continuous(limits = c(0, 12),
                     breaks = seq(0, 12, by = 2))+
  scale_x_continuous(breaks = seq(1, 6, by = 1), labels=c("0-7 cm", "7.1-15 cm", "15.1-30 cm", "30.2-50 cm", "50.1-80 cm",">80 cm"))+
# Size 1
  annotate('text', x = 0.72, y = 2.3, label = '45')+
  annotate('text', x = 0.87, y = 4.4, label = '97')+
  annotate('text', x = 1.0, y = 1.3, label = '19')+
  annotate('text', x = 1.15, y = 2.3, label = '44')+
  annotate('text', x = 1.29, y = 0.6, label = '2')+
# Size 2
  annotate('text', x = 1.72, y = 2.7, label = '55')+
  annotate('text', x = 1.85, y = 7.3, label = '167')+
  annotate('text', x = 2.00, y = 2.8, label = '56')+
  annotate('text', x = 2.15, y = 3.5, label = '74')+
  annotate('text', x = 2.28, y = 1.1, label = '14')+
# Size 3
  annotate('text', x = 2.72, y = 4.5, label = '98')+
  annotate('text', x = 2.85, y = 10.8, label = '255')+
  annotate('text', x = 3.0, y = 5.3, label = '118')+
  annotate('text', x = 3.15, y = 5.5, label = '123')+
  annotate('text', x = 3.28, y = 1.3, label = '20')+
# Size 4
  annotate('text', x = 3.72, y = 3.6, label = '76')+
  annotate('text', x = 3.85, y = 9.5, label = '222')+
  annotate('text', x = 4.0, y = 4.1, label = '89')+
  annotate('text', x = 4.15, y = 5.1, label = '113')+
  annotate('text', x = 4.28, y = 1.2, label = '17')+
# Size 5
  annotate('text', x = 4.72, y = 2.8, label = '57')+
  annotate('text', x = 4.85, y = 5.9, label = '134')+
  annotate('text', x = 5.00, y = 3.1, label = '63')+
  annotate('text', x = 5.15, y = 3.3, label = '68')+
  annotate('text', x = 5.28, y = 1.0, label = '12')+
# Size 6
  annotate('text', x = 5.72, y = 2.6, label = '52')+
  annotate('text', x = 5.85, y = 6.9, label = '158')+
  annotate('text', x = 6.00, y = 4.4, label = '96')+
  annotate('text', x = 6.15, y = 4.4, label = '97')+
  annotate('text', x = 6.28, y = 1.1, label = '14')+
  theme_classic()+
  scale_fill_manual(values= province.cols)+
  theme(panel.grid = element_blank(),
        text = element_text(size = 15),
        legend.position = 'none')+
  annotate('text', x = 6, y = 12, label = 'c)')

# Diet
Diet<- ggplot(data = traits.province,
              aes(x = Diet, fill = Province)) +
  geom_bar(aes(y = after_stat(count / sum(count))*100),
           position = position_dodge2(), color = 'black', alpha=0.3) +
  labs(x = 'Diet', y = '')+
  scale_y_continuous(limits = c(0, 20),
                     breaks = seq(0, 20, by = 5))+
  scale_x_discrete(label= c("Piscivores", "Herbivores-Detritivores", "Invertivores (mobile)", "Invertivores (sessile)","Omnivores", "Planktivores"))+
  # FC diet
  annotate('text', x = 0.64, y = 6.1, label = '127')+
  annotate('text', x = 0.82, y = 15.3, label = '352')+
  annotate('text', x = 1.01, y = 8.1, label = '175')+
  annotate('text', x = 1.19, y = 9.1, label = '200')+
  annotate('text', x = 1.36, y = 2.1, label = '28')+
  # HD diet
  annotate('text', x = 1.64, y = 1.4, label = '11')+
  annotate('text', x = 1.82, y = 3.9, label = '71')+
  annotate('text', x = 2.01, y = 2.5, label = '37')+
  annotate('text', x = 2.19, y = 2.3, label = '33')+
  annotate('text', x = 2.36, y = 1.1, label = '3')+
  # IM diet
  annotate('text', x = 2.64, y = 7.3, label = '156')+
  annotate('text', x = 2.82, y = 17.5, label = '407')+
  annotate('text', x = 3.01, y = 6.7, label = '140')+
  annotate('text', x = 3.19, y = 8.2, label = '177')+
  annotate('text', x = 3.36, y = 2.3, label = '32')+
  # IS diet
  annotate('text', x = 3.64, y = 1.2, label = '6')+
  annotate('text', x = 3.82, y = 2.0, label = '24')+
  annotate('text', x = 4.01, y = 1.4, label = '11')+
  annotate('text', x = 4.19, y = 1.5, label = '13')+
  annotate('text', x = 4.36, y = 1.1, label = '3')+
  # OM diet
  annotate('text', x = 4.64, y = 1.6, label = '16')+
  annotate('text', x = 4.82, y = 2.7, label = '43')+
  annotate('text', x = 5.01, y = 1.8, label = '20')+
  annotate('text', x = 5.19, y = 2.1, label = '27')+
  annotate('text', x = 5.36, y = 1.2, label = '6')+
  # Pk diet
  annotate('text', x = 5.64, y = 3.7, label = '67')+
  annotate('text', x = 5.82, y = 6.5, label = '136')+
  annotate('text', x = 6.01, y = 3.3, label = '58')+
  annotate('text', x = 6.19, y = 3.8, label = '69')+
  annotate('text', x = 6.36, y = 1.3, label = '7')+
  theme_classic()+
  scale_fill_manual(values=province.cols)+
  theme(panel.grid = element_blank(),
        text = element_text(size = 15),
        legend.position = 'none')+
  annotate('text', x = 6, y = 20, label = 'd)')

grid.arrange(Size, Diet, nrow = 2) # 1280x1024

# Mobility
Mob<- ggplot(data = traits.province,
             aes(x = Mobility, fill = Province)) +
  geom_histogram(aes(y = ifelse(after_stat(count) > 0,
                                after_stat(count / sum(count))*100,
                                NA)),
                 binwidth = 0.5, position = position_dodge(0.7), color = 'black', alpha=0.3) +
  labs(x = 'Mobility', y = 'Relative frequency (%)')+
  scale_y_continuous(limits = c(0, 17.5),
                     breaks = seq(0, 17.5, by = 2.5))+
  scale_x_continuous(breaks = seq(1, 4, by = 1), label=c("Sedentary", "Mobile within-reef", "Mobile among reefs","Very mobile"))+
# Mobility 1
  annotate('text', x = 0.73, y = 7.3, label = '156')+
  annotate('text', x = 0.86, y = 17.1, label = '398')+
  annotate('text', x = 1.0, y = 6.7, label = '142')+
  annotate('text', x = 1.14, y = 8.1, label = '176')+
  annotate('text', x = 1.28, y = 2.3, label = '33')+
# Mobility 2
  annotate('text', x = 1.73, y = 5.3, label = '107')+
  annotate('text', x = 1.86, y = 10.0, label = '223')+
  annotate('text', x = 2.0, y = 4.6, label = '89')+
  annotate('text', x = 2.14, y = 4.9, label = '97')+
  annotate('text', x = 2.28, y = 1.6, label = '16')+
# Mobility 3
  annotate('text', x = 2.73, y = 3.7, label = '66')+
  annotate('text', x = 2.86, y = 10.6, label = '238')+
  annotate('text', x = 3.0, y = 4.4, label = '84')+
  annotate('text', x = 3.14, y = 6.4, label = '134')+
  annotate('text', x = 3.28, y = 1.6, label = '15')+
# Mobility 4
  annotate('text', x = 3.73, y = 3.2, label = '54')+
  annotate('text', x = 3.86, y = 8.0, label = '174')+
  annotate('text', x = 4.0, y = 6.1, label = '126')+
  annotate('text', x = 4.14, y = 5.5, label = '112')+
  annotate('text', x = 4.28, y = 1.6, label = '15')+
  theme_classic()+
  scale_fill_manual(values=province.cols)+
  theme(panel.grid = element_blank(),
        text = element_text(size = 15),
        legend.position = 'none')+
  annotate('text', x = 4, y = 17.5, label = 'a)')

# Gregariousness
Greg<- ggplot(data = traits.province,
              aes(x = Gregariousness, fill = Province)) +
  geom_histogram(aes(y = ifelse(after_stat(count) > 0,
                                after_stat(count / sum(count))*100,
                                NA)),
                 binwidth = 0.5, position = position_dodge(0.7), color = 'black', alpha=0.3) +
  labs(x = 'Gregariousness', y = '')+
  scale_y_continuous(limits = c(0, 30),
                     breaks = seq(0, 30, by = 10))+
  scale_x_continuous(breaks = seq(1, 4, by = 1), label=c("Solitary", "Pairing", "Small group (3-50)", "large group (>50)"))+
# Gregariousness 1
  annotate('text', x = 0.73, y = 12.6, label = '275')+
  annotate('text', x = 0.86, y = 29.2, label = '684')+
  annotate('text', x = 1.0, y = 12.1, label = '263')+
  annotate('text', x = 1.14, y = 14.8, label = '329')+
  annotate('text', x = 1.28, y = 3.1, label = '40')+
# Gregariousness 2
  annotate('text', x = 1.73, y = 1.6, label = '3')+
  annotate('text', x = 1.86, y = 2.0, label = '12')+
  annotate('text', x = 2.0, y = 1.9, label = '9')+
  annotate('text', x = 2.14, y = 1.8, label = '7')+
  annotate('text', x = 2.28, y = 1.5, label = '1')+
  # Gregariousness 3
  annotate('text', x = 2.73, y = 4.7, label = '80')+
  annotate('text', x = 2.86, y = 11.7, label = '252')+
  annotate('text', x = 3.0, y = 6.8, label = '131')+
  annotate('text', x = 3.14, y = 6.8, label = '130')+
  annotate('text', x = 3.28, y = 2.8, label = '32')+
# Gregariousness 4
  annotate('text', x = 3.73, y = 2.5, label = '25')+
  annotate('text', x = 3.86, y = 4.9, label = '85')+
  annotate('text', x = 4.0, y = 3.0, label = '38')+
  annotate('text', x = 4.14, y = 3.6, label = '53')+
  annotate('text', x = 4.28, y = 1.7, label = '6')+
  theme_classic()+
  scale_fill_manual(values=province.cols)+
  theme(panel.grid = element_blank(),
        text = element_text(size = 15),
        legend.position = 'none')+
  annotate('text', x = 4, y = 30, label = 'b)')

grid.arrange(Mob, Greg, nrow = 2) # 1280x1024


# Position
Pos<- ggplot(data = traits.province,
             aes(x = Position, fill = Province)) +
  geom_histogram(aes(y = ifelse(after_stat(count) > 0,
                                after_stat(count / sum(count))*100,
                                NA)),
                 binwidth = 0.5, position = position_dodge(0.7), color = 'black', alpha=0.3) +
  labs(x = 'Position', y = 'Relative frequency (%)')+
  scale_y_continuous(limits = c(0, 30),
                     breaks = seq(0, 30, by = 5))+
  scale_x_continuous(breaks = seq(1, 3, by = 1), label=c("Benthic", "Benthopelagic", "Pelagic"))+
# Position 1
  annotate('text', x = 0.72, y = 9.8, label = '204')+
  annotate('text', x = 0.86, y = 29.1, label = '682')+
  annotate('text', x = 1.00, y = 10.5, label = '223')+
  annotate('text', x = 1.14, y = 14.2, label = '314')+
  annotate('text', x = 1.28, y = 3.3, label = '45')+
# Position 2
  annotate('text', x = 1.72, y = 6.3, label = '119')+
  annotate('text', x = 1.86, y = 9.4, label = '194')+
  annotate('text', x = 2.00, y = 5.9, label = '108')+
  annotate('text', x = 2.14, y = 6.3, label = '119')+
  annotate('text', x = 2.28, y = 2.3, label = '19')+
# Position 3
  annotate('text', x = 2.72, y = 3.9, label = '60')+
  annotate('text', x = 2.86, y = 7.9, label = '157')+
  annotate('text', x = 3.00, y = 6.0, label = '110')+
  annotate('text', x = 3.14, y = 5.0, label = '86')+
  annotate('text', x = 3.28, y = 2.1, label = '15')+
  theme_classic()+
  scale_fill_manual(values=province.cols)+
  theme(panel.grid = element_blank(),
        text = element_text(size = 15),
        legend.position = 'bottom')+
  annotate('text', x = 3.28, y = 30, label = 'e)')

# Activity
Act<- ggplot(data = traits.province,
             aes(x = Activity, fill = Province)) +
  geom_bar(aes(y = after_stat(count / sum(count))*100),
           position = position_dodge2(), color = 'black', alpha=0.3) +  labs(x = 'Activity', y = '')+
  scale_y_continuous(limits = c(0, 40),
                     breaks = seq(0, 40, by = 5))+
  scale_x_discrete(label=c("Diurnal", "Nocturnal"))+
# Activity 1
  annotate('text', x = 0.64, y = 14.4, label = '306')+
  annotate('text', x = 0.83, y = 35.4, label = '826')+
  annotate('text', x = 1.00, y = 15.5, label = '333')+
  annotate('text', x = 1.18, y = 18.0, label = '394')+
  annotate('text', x = 1.36, y = 4.2, label = '55')+
# Activity 2
  annotate('text', x = 1.64, y = 5.1, label = '77')+
  annotate('text', x = 1.83, y = 10.4, label = '207')+
  annotate('text', x = 2.00, y = 6.4, label = '108')+
  annotate('text', x = 2.18, y = 7.1, label = '125')+
  annotate('text', x = 2.36, y = 3.0, label = '24')+
  theme_classic()+
  scale_fill_manual(values=province.cols)+
  theme(panel.grid = element_blank(),
        text = element_text(size = 15),
        legend.position = 'bottom')+
  annotate('text', x = 2, y = 40, label = 'b)')

grid.arrange(Pos, Act, nrow = 2) # 1280x1024
grid.arrange(Mob, Greg,Size, Diet,Pos,Act, nrow=3) #2500x3072
ggsave("Figure 2.svg", grid.arrange(Mob, Greg,Size, Diet,Pos,Act, nrow=3), width = 3500, height= 2500, units='px', dpi=300, bg="white")

# We remove all objects that are no longer useful ####
rm(Pos, Mob, Act, Greg, Size, Diet, traits.province)


# Functional entities ####
# California-Peru province ####
spp.fes.caliperu<- data.frame(FE = str_c(traits.cali.peru$Size, '',
                                traits.cali.peru$Mobility, '',
                                traits.cali.peru$Activity, '',
                                traits.cali.peru$Gregariousness, '',
                                traits.cali.peru$Position, '',
                                traits.cali.peru$Diet))
rownames(spp.fes.caliperu)<- rownames(traits.cali.peru)
as_tibble(spp.fes.caliperu)

# Avoiding repeated FE and changing row names from
# fd.coord object
sub.calif.peru<- as.matrix(calif.peru[, -c(6:9)])
rownames(sub.calif.peru) = spp.fes.caliperu[, ]
sub.calif.peru<- sub.calif.peru[!duplicated(rownames(sub.calif.peru)), ]

dim(sub.calif.peru) # 193 unique FEs


# Cortez province ####
spp.fes.cort<- data.frame(FE = str_c(traits.cort$Size, '',
                                    traits.cort$Mobility, '',
                                    traits.cort$Activity, '',
                                    traits.cort$Gregariousness, '',
                                    traits.cort$Position, '',
                                    traits.cort$Diet))
rownames(spp.fes.cort)<- rownames(traits.cort)
as_tibble(spp.fes.cort)

# Avoiding repeated FE and changing row names from
# fd.coord object
sub.cortez<- as.matrix(cortez[, -c(6:9)])
rownames(sub.cortez) = spp.fes.cort[, ]
sub.cortez<- sub.cortez[!duplicated(rownames(sub.cortez)), ]

dim(sub.cortez) # 348 unique FEs

# Ocean Islands province ####
spp.fes.oce<- data.frame(FE = str_c(traits.oce$Size, '',
                                    traits.oce$Mobility, '',
                                    traits.oce$Activity, '',
                                    traits.oce$Gregariousness, '',
                                    traits.oce$Position, '',
                                    traits.oce$Diet))
rownames(spp.fes.oce)<- rownames(traits.oce)
as_tibble(spp.fes.oce)

# Avoiding repeated FE and changing row names from
# fd.coord object
sub.oceanic<- as.matrix(oceanic[, -c(6:9)])
rownames(sub.oceanic) = spp.fes.oce[, ]
sub.oceanic<- sub.oceanic[!duplicated(rownames(sub.oceanic)), ]

dim(sub.oceanic) # 226 unique FEs


# Panamic province ####
spp.fes.pana<- data.frame(FE = str_c(traits.panam$Size, '',
                                     traits.panam$Mobility, '',
                                     traits.panam$Activity, '',
                                     traits.panam$Gregariousness, '',
                                     traits.panam$Position, '',
                                     traits.panam$Diet))
rownames(spp.fes.pana)<- rownames(traits.panam)
as_tibble(spp.fes.pana)

# Avoiding repeated FE and changing row names from
# fd.coord object
sub.panamic<- as.matrix(panamic[, -c(6:9)])
rownames(sub.panamic) = spp.fes.pana[, ]
sub.panamic<- sub.panamic[!duplicated(rownames(sub.panamic)), ]

dim(sub.panamic) # 245 unique FEs

# Backbone ####
spp.fes.back<- data.frame(FE = str_c(traits.backbone$Size, '',
                                     traits.backbone$Mobility, '',
                                     traits.backbone$Activity, '',
                                     traits.backbone$Gregariousness, '',
                                     traits.backbone$Position, '',
                                     traits.backbone$Diet))
rownames(spp.fes.back)<- rownames(traits.backbone)
as_tibble(spp.fes.back)

# Avoiding repeated FE and changing row names from
# fd.coord object
sub.backbone<- as.matrix(backbone[, -c(6:9)])
rownames(sub.backbone) = spp.fes.back[, ]
sub.backbone<- sub.backbone[!duplicated(rownames(sub.backbone)), ]

dim(sub.backbone) # 57 unique FEs

# We remove all objects that are no longer useful ####
rm(spp.fes.caliperu, spp.fes.cort, spp.fes.oce, spp.fes.pana,
   spp.fes.back, traits.cali.peru, traits.cort, traits.oce,
   traits.panam, traits.backbone, types, presence,
   presence.conditions)


# Convex hull based on the axes 1 and 2 ####
regional.hull <- hull12.fn(regional.space)
backbone.hull<- hull12.fn(backbone)

# Graphics
PC12_backbone<- ggplot(regional.space, aes (x = PC1, y = PC2))+
  geom_point(data = regional.space, aes (x = PC1, y = PC2),
             size = 0.5, shape = 3, colour = "grey")+
  geom_polygon(data = regional.hull, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0,
               col = "black")+
  geom_polygon(data = backbone.hull, aes (x = PC1, y = PC2),
               linetype = 1, linewidth = 1, alpha = 0.3,
               col = province.cols[5], fill = province.cols[5])+
  labs(x = 'PCoA 1 (29%)', y = '')+
  theme_bw()+annotate('text', x = 0.25, y = 0.5, label = 'b)')+
  theme(panel.grid = element_blank(),
        text = element_text(size = 25))



# Convex hull based on the axes 3 and 4 ####
regional.hull <- hull34.fn(regional.space)
backbone.hull<- hull34.fn(backbone)

# Graphics ####
PC34_backbone<- ggplot(regional.space, aes (x = PC3, y = PC4))+
  geom_point(data = regional.space, aes (x = PC3, y = PC4),
             size = 0.5, shape = 3, colour = "grey")+
  geom_polygon(data = regional.hull, aes (x = PC3, y = PC4),
               linetype = 1, linewidth = 1, alpha = 0,
               col = "black")+
  geom_polygon(data = backbone.hull, aes (x = PC3, y = PC4),
               linetype = 1, linewidth = 1, alpha = 0.3,
               col = province.cols[5], fill = province.cols[5])+
  labs(x = 'PCoA 3 (13.11%)', y = '')+
  theme_bw()+annotate('text', x = 0.6, y = 0.25, label = 'c)')+
  theme(panel.grid = element_blank(),
        text = element_text(size = 25))


# Regional plots ####
# Col 1= Species richness; Col 2= Relative species richness;
# Col 3= Number of FEs; Col 4= Relaive number of FEs
# Col 5= Functional volume
metrics<- c('Richness', 'Relative richness',
            'FEs', 'Relative FEs', 'Volume')
provinces<- c(provinces, 'Backbone')
regional.volume<- data.frame(Province = c(rep(provinces, each = 5)),
                             Metric = c(metrics),
    Value = c(383, (383/1235)*100, 199, (199/405)*100, chCaliPeru,
              1033,(1033/1235)*100, 360, (360/405)*100, chCor,
              441, (441/1235)*100, 236, (236/405)*100, chOce,
              519, (519/1235)*100, 255, (255/405)*100, chPan,
              79, (79/1235)*100, 59, (59/405)*100, chBack))
regional.volume$Metric<- factor(regional.volume$Metric,
                    levels = c('S', 'Relative richness',
                               'FE', 'Relative FEs', 'Fvol'),
                     ordered = T)


# Histogram plots
# California-Peruvian province
hist.caliperu<- regional.volume %>%
  filter(Province == 'Californian') %>% .[c(2,4:5),] %>%
  ggplot(data = ., aes(x = Metric, y = Value, fill = Province)) +
  geom_col(position = position_dodge(), color = 'black',
           linewidth = 1.05, alpha=0.3, show.legend = F) +
  labs(subtitle = expression(Californian),
       x = NULL, y = 'Relative richness (%)')+
  scale_x_discrete(labels = c('No. spp', 'No. FEs', 'Volume'))+
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, by = 20))+
  annotate('text', x = 1, y = 38, label = '383', size = 7)+
  annotate('text', x = 2, y = 56.1, label = '199', size = 7)+
  annotate('text', x = 3, y = 77, label = '70', size = 7)+
  theme_classic()+
  scale_fill_manual(values = province.cols[1])+
  theme(panel.grid = element_blank(),
        text = element_text(size = 25))

# Cortez province
hist.cortez<- regional.volume %>%
  filter(Province == 'Cortez') %>% .[c(2,4:5),] %>%
  ggplot(data = ., aes(x = Metric, y = Value, fill = Province)) +
  geom_col(position = position_dodge(), color = 'black',
           linewidth = 1.05, alpha=0.3, show.legend = F) +
  labs(subtitle = expression(Cortez),
       x = NULL, y = '')+
  scale_x_discrete(labels = c('No. spp', 'No. FEs', 'Volume'))+
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, by = 20))+
  annotate('text', x = 1, y = 90.6, label = '1033', size = 7)+
  annotate('text', x = 2, y = 95.9, label = '360', size = 7)+
  annotate('text', x = 3, y = 100.2, label = '98', size = 7)+
  theme_classic()+
  scale_fill_manual(values = province.cols[2])+
  theme(panel.grid = element_blank(),
        text = element_text(size = 25))

# Ocean Islands province
hist.oceanic<- regional.volume %>%
  filter(Province == 'Oceanic') %>% .[c(2,4:5),] %>%
  ggplot(data = ., aes(x = Metric, y = Value, fill = Province)) +
  geom_col(position = position_dodge(), color = 'black',
           linewidth = 1.05, alpha=0.3, show.legend = F) +
  labs(subtitle = expression(Ocean~Islands),
       x = NULL, y = '')+
  scale_x_discrete(labels = c('No. spp', 'No. FEs', 'Volume'))+
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, by = 20))+
  annotate('text', x = 1, y = 42.7, label = '441', size = 7)+
  annotate('text', x = 2, y = 65.3, label = '236', size = 7)+
  annotate('text', x = 3, y = 84, label = '77', size = 7)+
  theme_classic()+
  scale_fill_manual(values = province.cols[3])+
  theme(panel.grid = element_blank(),
        text = element_text(size = 25))

# Panamic province
hist.panamic<- regional.volume %>%
  filter(Province == 'Panamic') %>% .[c(2,4:5),] %>%
  ggplot(data = ., aes(x = Metric, y = Value, fill = Province)) +
  geom_col(position = position_dodge(), color = 'black',
           linewidth = 1.05, alpha=0.3, show.legend = F) +
  labs(subtitle = expression(Panamic),
       x = NULL, y = '')+
  scale_x_discrete(labels = c('No. spp', 'No. FEs', 'Volume'))+
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, by = 20))+
  annotate('text', x = 1, y = 49, label = '519', size = 7)+
  annotate('text', x = 2, y = 70, label = '255', size = 7)+
  annotate('text', x = 3, y = 87, label = '80', size = 7)+
  theme_classic()+
  scale_fill_manual(values = province.cols[4])+
  theme(panel.grid = element_blank(),
        text = element_text(size = 25))

# Backbone
hist.backbone<- regional.volume %>%
  filter(Province == 'Backbone') %>% .[c(2,4:5),] %>%
  ggplot(data = ., aes(x = Metric, y = Value, fill = Province)) +
  geom_col(position = position_dodge(), color = 'black',
           linewidth = 1.05, alpha=0.3, show.legend = F) +
  labs(subtitle = expression(Backbone),
       x = NULL, y = '')+
  scale_x_discrete(labels = c('No. spp', 'No. FEs', 'Volume'))+
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, by = 20))+
  annotate('text', x = 1, y = 13.4, label = '79', size = 7)+
  annotate('text', x = 2, y = 21.6, label = '59', size = 7)+
  annotate('text', x = 3, y = 42, label = '35', size = 7)+
  theme_classic()+annotate('text', x = 3, y = 100, label = 'a)')+
  scale_fill_manual(values = province.cols[5])+
  theme(panel.grid = element_blank()+annotate('text', x = 3, y = 100, label = 'a)'),
        text = element_text(size = 25))

# Histograms and volume plots ####
grid.arrange(hist.caliperu, hist.cortez, hist.oceanic,
             hist.panamic, hist.backbone, PC12_CalPeru, PC12_Cor,
             PC12_Oce, PC12_Pan, PC12_backbone, PC34_CalPeru,
             PC34_Cor, PC34_Oce, PC34_Pan, PC34_backbone,
             nrow = 3) # 2560x1440
ggsave("Figure 3.svg", width = 2560, height= 1440, units='px', dpi=300, bg="white")

