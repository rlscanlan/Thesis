setwd("/Users/rebekahscanlon/Desktop/CellularSen/datafiles")


#A files
A_DDIS = read.csv("A_DDIS.csv")
A_DDIS_p53 = read.csv("A_DDIS_p53KD.csv")
A_OIS = read.csv("A_OIS.csv")
A_OIS_p53 = read.csv("A_OIS_p53KD.csv")

#B files
B_DDIS = read.csv("B_DDIS.csv")
B_DDIS_p53 = read.csv("B_DDIS_p53kd.csv")
B_DDIS_nfkb = read.csv("B_DDIS_nfkbkd.csv")
B_OIS = read.csv("B_OIS.csv")
B_OIS_p53 = read.csv("B_OIS_p53kd.csv")
B_OIS_nfkb = read.csv("B_OIS_nfkbkd.csv")

#full model files
DDIS = read.csv("DDIS.csv")
DDIS_p53 = read.csv("DDIS_P53.csv")
DDIS_nfkb = read.csv("DDIS_NFKB.csv")
OIS = read.csv("OIS.csv")
OIS_p53 = read.csv("OIS_P53.csv")
OIS_nfkb = read.csv("OIS_NFKB.csv")


sentype = c("OIS")
B_OIS$sentype = sentype

library(tidyverse)
library(tidyr)
library(dplyr)

####newgraphs####
#a
ggplot()+
  geom_line(data = A_DDIS, aes(x = time, y = p38, linetype = sentype, color = " None"), size =1.5) +
  geom_line(data = A_OIS, aes(x = time, y = p38, linetype = sentype, color = " None"), size =1.5)+
  geom_line(data = A_DDIS_p53, aes(x = time, y = p38, linetype = sentype, color = " p53"), size =1.5) +
  geom_line(data = A_OIS_p53, aes(x = time, y = p38, linetype = sentype, color = " p53"), size =1.5)+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 30), axis.title.x=element_text(size=17),
        legend.text=element_text(size=17), axis.title.y=element_text(size=20),
        legend.title=element_text(size=20), axis.text.x=element_text(size=17),
        axis.text.y=element_text(size=17)) +
  labs(title = "p38", x = "Time (AU)", y = "Expression (AU)", linetype = "Senescence", color = "Gene inhibition")+
  ylim(0,5.5)+
  xlim(30,100)

#b and c
ggplot()+
  geom_line(data = DDIS, aes(x = time, y = IL1a, linetype = sentype, color = "r"), size =1.5) +
  geom_line(data = OIS, aes(x = time, y = IL1a, linetype = sentype, colour = "t"), size =1.5)+
  #geom_line(data = DDIS, aes(x = time, y = NICD, linetype = sentype, color = "NICD"), size =1.5) +
  #geom_line(data = OIS, aes(x = time, y = NICD, linetype = sentype, color = "NICD"), size =1.5)+
  #geom_line(data = DDIS_nfkb, aes(x = time, y = pNFkB, linetype = sentype, color = "RELA"), size =1.5) +
  #geom_line(data = OIS_nfkb, aes(x = time, y = pNFkB, linetype = sentype, color = "RELA"), size =1.5)+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 30), axis.title.x=element_text(size=17),
        legend.text=element_text(size=17), axis.title.y=element_text(size=20),
        legend.title=element_text(size=20), axis.text.x=element_text(size=17),
        axis.text.y=element_text(size=17)) +
  labs(title = "IL1a", x = "Time (AU)", y = "Expression (AU)", linetype = "Senescence")+
  ylim(0,5.5)+
  xlim(30,100)


#geom_line(data = B_DDIS_p53, aes(x = time, y = pp38, linetype = sentype, color = " p53"), size =1.5) +
#geom_line(data = B_OIS_p53, aes(x = time, y = pp38, linetype = sentype, color = " p53"), size =1.5)+
####normal sen sims####


ggplot()+
  geom_line(data = DDIS, aes(x = time, y = NOTCH1, linetype = " NOTCH", color = sentype), size =1.5) +
  geom_line(data = OIS, aes(x = time, y = NOTCH1, linetype = " NOTCH",  color = sentype), size =1.5)+
  geom_line(data = DDIS, aes(x = time, y = NICD, linetype = "NICD",  color = sentype), size =1.5) +
  geom_line(data = OIS, aes(x = time, y = NICD, linetype = "NICD", color = sentype), size =1.5)+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 20), axis.title.x=element_text(size=15),
        legend.text=element_text(size=13), axis.title.y=element_text(size=15),
        legend.title=element_text(size=15), axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12)) +
  labs(title = "Notch", x = "Time (AU)", y = "Expression (AU)", linetype = "NOTCH state", colour = "Senescence")+
  ylim(0,5.5)+
  xlim(30,100)


ggplot()+
  geom_line(data = DDIS, aes(x = time, y = IL6, color = sentype), size =1.5) +
  geom_line(data = OIS, aes(x = time, y = IL6, color = sentype), size = 1.5)+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 20), axis.title.x=element_text(size=15),
        legend.text=element_text(size=13), axis.title.y=element_text(size=15),
        legend.title=element_text(size=15), axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12)) +
  labs(title = "IL6", x = "Time (AU)", y = "Expression (AU)", colour = "Senescence")+
  ylim(0,5.5)+
  xlim(30,100)

#6x7.5

####p53kd#####

ggplot()+
  geom_line(data = DDIS, aes(x = time, y = p53, color = sentype, linetype = "None"), size =1.5) +
  geom_line(data = OIS, aes(x = time, y = p53, color = sentype, linetype = "None"), size =1.5)+
  geom_line(data = DDIS_p53, aes(x = time, y = p53, color = sentype, linetype = "p53"), size =1.5) +
  geom_line(data = OIS_p53, aes(x = time, y = p53, color = sentype, linetype = "p53"), size =1.5)+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 20), axis.title.x=element_text(size=15),
        legend.text=element_text(size=13), axis.title.y=element_text(size=15),
        legend.title=element_text(size=15), axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12)) +
  labs(title = "p53", x = "Time (AU)", y = "Expression (AU)", colour = "Senescence", linetype = "Gene inhibition")+
  ylim(0,6.1)+
  xlim(30,100)

####nfkb KD####

ggplot()+
  geom_line(data = DDIS, aes(x = time, y = pp53, color = sentype, linetype = " None"), size =1.5) +
  geom_line(data = OIS, aes(x = time, y = pp53, color = sentype, linetype = " None"), size =1.5)+
  geom_line(data = DDIS_nfkb, aes(x = time, y = pp53, color = sentype, linetype = "NF-kB"), size =1.5) +
  geom_line(data = OIS_nfkb, aes(x = time, y = pp53, color = sentype, linetype = "NF-kB"), size =1.5)+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 20), axis.title.x=element_text(size=15),
        legend.text=element_text(size=13), axis.title.y=element_text(size=15),
        legend.title=element_text(size=15), axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12)) +
  labs(title = "pp53", x = "Time (AU)", y = "Expression (AU)", colour = "Senescence", linetype = "Gene inhibition")+
  ylim(0,5.5)+
  xlim(30,100)
