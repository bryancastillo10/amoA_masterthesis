library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(tidyr)


rel_abun<-read.csv("raw_data/rel_amoA.csv", fileEncoding="UTF-8-BOM")
species_col <- colnames(rel_abun)[2:length(rel_abun)]

dat.aggregate <-pivot_longer(rel_abun,
                             cols=species_col,
                             names_to="amoA_Taxa",
                             values_to="Relative_Abundance")

dat.aggregate$ID<-factor(dat.aggregate$ID,levels=c("T1","T2","T3","T4",
                                                   "C1","C2","C3","C4","C5",
                                                   "PCB","LCD1","LCD2","R2.O","R1.A",
                                                   "R1.O","SBR0","SBR7"))

abundancePlot <- ggplot(dat.aggregate,aes(x=ID,y=Relative_Abundance))+
  geom_bar(aes(fill=amoA_Taxa),stat="identity",color="black",size=0.3)+
  facet_grid(~WWTP,scale="free",space="free_x")+
  scale_fill_manual(name="amoA Taxa",labels=c(
    "Nitrosomonas europaea ATCC 19718,ATCC 25978,Nm35,Nm50",
    "Nitrosomonas eutropha C91,Nm57,Nm56",
    "Nitrosomonas eutropha Nm38 Nm19 Nm14",
    "Nitrosomonas eutropha Nm24",
    "Nitrosomonas halophila Nm1",
    "Nitrosomonas nitrosa Nm91,Nm146,Nm90",
    "Nitrosomonas oligotropha Nm49",
    "Nitrosomonas sp. Nm34",
    "Nitrosomonas stercoris KYUHINS",
    "Nitrosospira multiformis ATCC 25196,Nl13,ATCC 25196",
    "Nitrosospira multiformis Nl18",
    "Nitrosospira sp. Nsp14"
  ),
    values=c(brewer.pal(3,"Dark2"),hcl.colors(5,"Batlow")
                        ,"#E69F00",brewer.pal(3,"Set2"))
  )+
  theme_classic()+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size=9,color="black",face="italic"),
    panel.spacing.x = unit(3, "mm"),
    axis.text.x = element_text(angle = 45, hjust = 1,size=10,color="black"),
    axis.text.y=element_text(size=15),
    axis.title.y=element_text(size=15),
    legend.title=element_text(size=rel(1.0)),
    legend.text=element_text(size=rel(1))
  )+
  scale_x_discrete(expand=c(0.1,0.5))+
  scale_y_continuous(expand=c(0,0))+
  labs(
    x="",
    y="Relative amoA Abundance (%)")

ggsave("plot_output/abundancePlot.png", height=7, width=14)
