library(ggplot2)
library(dplyr)
library(ggpubr)
qPCR<-read.csv("raw_data/genecopyamoA.csv", fileEncoding="UTF-8-BOM")

qPCR$WWTP<- ordered(qPCR$WWTP, levels=c("TFC","CSC", "PCB CSTR","PCB SBR","TSC"))

newqPCR<-qPCR %>%
  group_by(WWTP) %>%
  mutate_at(2,log10) %>%
  summarize(
    GeneCopy= mean(amoA.copy.mL, na.rm=TRUE),
    Removal= mean(NH3.Removal.Efficiency,na.rm=TRUE),
    SD1 =  sd(amoA.copy.mL, na.rm=TRUE),
    SD2 =  sd(NH3.Removal.Efficiency, na.rm=TRUE)
  )

qPCRPlot <-ggplot(newqPCR) + 
    geom_bar(aes(x = WWTP, y = GeneCopy, fill = "amoA Gene Abundance"), 
            stat = "identity", width=1,
            colour = "black") +
    geom_errorbar(aes(x = WWTP, ymin = GeneCopy - SD1, ymax = GeneCopy + SD1),
                  width = 0.2, color="midnightblue",
                  position = "dodge") +
  
    # geom_line(aes(x = WWTP, y = Removal * max(GeneCopy),group=1), 
    #         color = 'black', size = 2) +
    
    geom_point(aes(x = WWTP, y = Removal * max(GeneCopy), color="NH3 Removal Efficiency"), 
                size = 3)+
    
    geom_errorbar(aes(x = WWTP, ymin = Removal * max(GeneCopy) - SD2, 
                      ymax = Removal * max(GeneCopy) + SD2), 
                  colour = "red",
                  width = .2)+
    scale_fill_manual(name="amoA Gene Abundance",
      values=c("amoA Gene Abundance"="dodgerblue")
    )+
    scale_color_manual(name=expression(NH[3]~"Removal Efficiency"),
      values=c("NH3 Removal Efficiency"="black")
    )+

    theme_classic()+
    theme(
      legend.title=element_blank(),
      legend.text=element_text(size=rel(1),face="bold"),
      legend.position="top",
      legend.key.height = unit(0.1,"cm"),
      axis.text.x = element_text(angle=45,vjust=1,hjust=1,color=
                                "black"),
      axis.text.y=element_text(size=rel(0.8))
      
    ) +
    scale_y_continuous(expand=c(0,0),
      position="left", name ="amoA Gene Abundance (log gene copy/mL sludge)",
      sec.axis = sec_axis(~./10,name=expression(NH[3]~"Removal Efficiency"),
                      breaks=c(0,0.2,0.4,0.6,0.8,1.0)))

# To save ggplot
# ggsave("plot_output/qPCR_amoA.png", height=8, width=10)