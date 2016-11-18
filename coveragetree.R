setwd("~/biclub/")
library(ggtree)
library(ape)
library(grid)
library(gridExtra)

tree <- read.tree("Tree1.nwk")

vital1 <- read.csv("vitalF1_R1_df.primermismatch.csv", header = TRUE, sep = "\t")
vital2 <- read.csv('vitalF2_R2_df.primermismatch.csv', header = TRUE, sep = '\t')
vital3 <- read.csv('vitalF3_R3_df.primermismatch.csv', header = TRUE, sep = '\t')

vtot <- merge(vital1, vital2)
vtot <- merge(vtot, vital3)

vtot$vitalFMM <- with(vtot, pmin(vitalF1_MM, vitalF2_MM, vitalF3_MM))
vtot$vitalRMM <- with(vtot, pmin(vitalR1_MM, vitalR2_MM, vitalR3_MM))

flint <- read.csv('Flint_df.primermismatch.csv', header = TRUE, sep = '\t')

mybut <- read.csv('Myprimers_df.primermismatch.csv', header = TRUE, sep = '\t')

metadata <- merge(flint, mybut)
metadata <- merge(metadata, vtot)

metadata$butTOT <- metadata$but672FWD_MM + metadata$but1031REV_MM
metadata$flintTOT <- metadata$BCoATscrF_MM + metadata$BCoATscrR_MM
metadata$vitalTOT <- metadata$vitalFMM + metadata$vitalRMM
colnames(metadata)[1] <- 'SeqID'
colnames(metadata)[4] <- 'funct'

metadata <- metadata[order(metadata$funct),]
metadata$number <- c(1:length(metadata$SeqID))



### My primers ###

p <- ggtree(tree, layout ="circular", size=.5)
p <- p %<+% metadata
p$data$nonfunct <- grepl('other', p$data$funct)
p$data$confirmed <- grepl('confirmed', p$data$funct)
p$data$greens <- grepl('^1$', p$data$butTOT)

p <- p + 
  geom_hilight(node = 1303, fill = 'grey', alpha = .7) + 
  geom_tippoint(aes(color=butTOT), size=2.7) + 
  scale_color_continuous(name="Total Mismatches", limits=c(0,10),
                         oob=scales::squish, low='darkgreen', high = 'red')


x1 <- p$data$x[p$data$confirmed=='TRUE']
y1 <- p$data$y[p$data$confirmed=='TRUE']
cfirmed<- as.data.frame(cbind(x1,y1))
x11 <- p$data$x[p$data$nonfunct=='TRUE']
y11 <- p$data$y[p$data$nonfunct=='TRUE']
nons<- as.data.frame(cbind(x11,y11))
p1 <- p + geom_point(data=cfirmed, x=x1+0.35, y=y1, fill="blue", shape=24, size=4, stroke=1) +
  geom_point(data=nons, x=x11+0.35, y=y11, fill="yellow", shape=24, size=4, stroke=1)
p1


### Vital Primers ###

p <- ggtree(tree, layout ="circular")
p <- p %<+% metadata
p$data$nonfunct <- grepl('other', p$data$funct)
p$data$confirmed <- grepl('confirmed', p$data$funct)
p$data$greens <- grepl('^1$', p$data$butTOT)

p <- p + 
  geom_hilight(node = 1303, fill = 'grey', alpha = .7) +
  geom_tippoint(aes(color=vitalTOT), size=2.7) + 
  scale_color_continuous(name="Total Mismatches", limits=c(0,10),
                         oob=scales::squish, low='darkgreen', high = 'red')


x1 <- p$data$x[p$data$confirmed=='TRUE']
y1 <- p$data$y[p$data$confirmed=='TRUE']
cfirmed<- as.data.frame(cbind(x1,y1))
x11 <- p$data$x[p$data$nonfunct=='TRUE']
y11 <- p$data$y[p$data$nonfunct=='TRUE']
nons<- as.data.frame(cbind(x11,y11))
p2 <- p + geom_point(data=cfirmed, x=x1+0.35, y=y1, fill="blue", shape=24, size=4, stroke=1) +
  geom_point(data=nons, x=x11+0.35, y=y11, fill="yellow", shape=24, size=4, stroke=1)
p2
### Flint Primers ###

p <- ggtree(tree, layout ="circular", size=.5, alpha=.7)
p <- p %<+% metadata

p$data$nonfunct <- grepl('other', p$data$funct)
p$data$confirmed <- grepl('confirmed', p$data$funct)
p$data$greens <- grepl('^1$', p$data$butTOT)

p <- p  + 
  geom_hilight(node = 1303, fill = 'grey', alpha = .7) + 
  geom_tippoint(aes(color=flintTOT), size=2.7) + 
  scale_color_continuous(name="Total Mismatches", limits=c(0,10),
                         oob=scales::squish, low='darkgreen', high = 'red')


x1 <- p$data$x[p$data$confirmed=='TRUE']
y1 <- p$data$y[p$data$confirmed=='TRUE']
cfirmed<- as.data.frame(cbind(x1,y1))
x11 <- p$data$x[p$data$nonfunct=='TRUE']
y11 <- p$data$y[p$data$nonfunct=='TRUE']
nons<- as.data.frame(cbind(x11,y11))
p3 <- p + geom_point(data=cfirmed, x=x1+0.35, y=y1, fill="blue", shape=24, size=4, stroke=1) +
  geom_point(data=nons, x=x11+0.35, y=y11, fill="yellow", shape=24, size=4, stroke=1)
p3






tiff('coveragetrees.tiff', width = 8.5, height = 5.4, units = 'in', res = 400)
multiplot(p3,p1,p2, ncol = 3)
grid.arrange(p3,p1,p2, ncol=3)


dev.off()

###################### FIND OUT WHICH NODE!!!!! ##################



p <- ggtree(tree, layout ="circular", size=.5, alpha=.7)
p <- p %<+% metadata
p +geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_point(data=cfirmed, x=x1+0.35, y=y1, fill="blue", shape=24, size=4, stroke=1) +
  geom_point(data=nons, x=x11+0.35, y=y11, fill="yellow", shape=24, size=4, stroke=1)


#################################### numbers for paper   #######################################
# this section calculates the numbers for primer validation results section
# counts the total number of tips 
library('caper')

total.tips <- sum(p$data$isTip)
good.tips <- clade.members(1661, tree)       # this clade contains all of the verified legitimate But coding sequence

megood <- sum(p$data$butTOT[p$data$node %in% good.tips] < 3)     
vitgood <- sum(p$data$vitalTOT[p$data$node %in% good.tips] < 3)
flintgood <- sum(p$data$flintTOT[p$data$node %in% good.tips] < 3)

metot <- sum(p$data$butTOT[p$data$isTip] < 3)
vittot <- sum(p$data$vitalTOT[p$data$isTip] < 3)
flinttot <- sum(p$data$flintTOT[p$data$isTip] < 3)


vitgood-megood



metot  # total number of sequences with 2 or fewer total mismatches with my primers
vittot # total number of sequences with 2 or fewer total mismatches with Vital's primers
flinttot # total number of sequences with 2 or fewer total mismatches with Flint's primers


metot - megood # the number of sequences that are likely to be amplified by my primers that are outside of the 'good clade'
vittot - vitgood # the number of sequences that are likely to be amplified by Vital's primers that are outside of the 'good clade'
flinttot-flintgood # the number of sequences that are likely to be amplified by Vital's primers that are outside of 'the good clade'


# now the ratios
megood/metot  *100  # 94% of the sequences likely to amplify with my primers are within the 'good clade'
vitgood/vittot *100 # 44% of the sequences likely to amplify with Vital's primers are within the 'good clade'
flintgood/flinttot *100 # 100%


