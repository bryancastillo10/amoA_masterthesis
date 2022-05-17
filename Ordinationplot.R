#Ordination Plot Using R

library(ggplot2)
library(vegan)
library(ggfortify)
library(dplyr)
library(ggsci)

#Import raw data first
spe<-read.csv("working_directory/raw_data/species_abun.csv",  fileEncoding="UTF-8-BOM")
envwater<-read.csv("working_directory/raw_data/envi.csv",  fileEncoding="UTF-8-BOM")


#Prepare the dependent variables by normalizing community data with log-transform
#As suggested by other research articles, log(x+1) and hellinger transformation can normalize community data
spe.mat<- spe[,c(3:12)]
spe.hell<-decostand(log1p(spe.mat),method="hellinger")

#Set site as categorical variable, you can arrange it in your own way by adding
#arguments in the factor() function
envwater$Site<-factor(envwater$Site)
Site<-envwater[,2]
envi<- envwater[,c(3:11)]

#PCA for species-site interaction: unconstrained analysis
dist.spe<-vegdist(spe.hell, method="euclidean")
mod<-betadisper(dist.spe,Site)
df_site<-data.frame(Site=rownames(mod$centroids),data.frame(mod$centroids))


#RDA for species-envi interaction: constrained analysis
rdamodel<-rda(spe.hell~.,envi)
df_species  <- data.frame(summary(rdamodel)$species[,1:2])
df_environ  <- scores(rdamodel, display = "bp") 

#Gather calculation results
#RDA scores and biplot point coordinates
df_species  <- data.frame(summary(rdamodel)$species[,1:2])
df_environ  <- scores(rdamodel, display = 'bp') 


#Returns the value of percentage of RDA scores
RDA1_varex<-round(summary(rdamodel)$cont$importance[2,1]*100,2)
RDA2_varex<-round(summary(rdamodel)$cont$importance[2,2]*100,2)

#for ggplot2 geom_segment scaling
scaling_factor<-2

#ggplot2
ggplot(df_species,aes(x=RDA1, y=RDA2))+
  geom_hline(yintercept=0,linetype="dashed")+
  geom_vline(xintercept=0,linetype="dashed")+
  coord_fixed()+
  
  #Plot Environmental Parameters
  # In the following parts, you can observe that I manually set up the
  #geom functions so that I can adjust their color,size,and shape,
  #per individual coordinates, so the settings would depend on your preferences
  #Check df_environ in the console to manually know the dataframe index/position of the coordinates
  geom_segment(data=df_environ[1,1:2,drop=FALSE], 
               aes(x=0, xend=RDA1*scaling_factor,y=0, 
                   yend=RDA2*scaling_factor,color="BOD")
               ,size=0.8,arrow=arrow(length=unit(0.02,"npc")))+
  geom_segment(data=df_environ[2,1:2,drop=FALSE], 
               aes(x=0, xend=RDA1*scaling_factor,y=0, 
                   yend=RDA2*scaling_factor,color="COD"),  
               size=0.8,arrow=arrow(length=unit(0.02,"npc")))+
  geom_segment(data=df_environ[3,1:2,drop=FALSE], 
               aes(x=0, xend=RDA1*scaling_factor,y=0, 
                   yend=RDA2*scaling_factor,color="NH3")
               ,size=0.8,arrow=arrow(length=unit(0.02,"npc")))+  
  geom_segment(data=df_environ[4,1:2,drop=FALSE], 
               aes(x=0, xend=RDA1*scaling_factor,y=0, 
                   yend=RDA2*scaling_factor,color="NO3")
               ,size=0.8,arrow=arrow(length=unit(0.02,"npc")))+
  geom_segment(data=df_environ[5,1:2,drop=FALSE], 
               aes(x=0, xend=RDA1*scaling_factor,y=0, 
                   yend=RDA2*scaling_factor,color="PO4")
               ,size=0.8,arrow=arrow(length=unit(0.02,"npc")))+
  geom_segment(data=df_environ[6,1:2,drop=FALSE], 
               aes(x=0, xend=RDA1*scaling_factor,y=0, 
                   yend=RDA2*scaling_factor,color="Cl")
               ,size=0.8,arrow=arrow(length=unit(0.02,"npc")))+
  geom_segment(data=df_environ[7,1:2,drop=FALSE], 
               aes(x=0, xend=RDA1*scaling_factor,y=0, 
                   yend=RDA2*scaling_factor,color="pH")  
               ,size=0.8,arrow=arrow(length=unit(0.02,"npc")))+
  geom_segment(data=df_environ[8,1:2,drop=FALSE], 
               aes(x=0, xend=RDA1*scaling_factor,y=0, 
                   yend=RDA2*scaling_factor,color="Temperature")
               ,size=0.8,arrow=arrow(length=unit(0.02,"npc")))+  
  geom_segment(data=df_environ[9,1:2,drop=FALSE], 
               aes(x=0, xend=RDA1*scaling_factor,y=0, 
                   yend=RDA2*scaling_factor,color="DO")
               ,size=0.8,arrow=arrow(length=unit(0.02,"npc")))+
  #Sites    
  geom_point(data=df_site[1,1:3,drop=FALSE],aes(x=PCoA1,y=PCoA2),color="gray24",size=1.5)+
  geom_text(data=df_site[1,1:3,drop=FALSE],aes(x=PCoA1,y=PCoA2),color="gray25",
            label="Bridge",size=2,face="bold",nudge_y=-0.06)+
  geom_point(data=df_site[2,1:3,drop=FALSE],aes(x=PCoA1,y=PCoA2),color="gray24",size=1.5)+
  geom_text(data=df_site[2,1:3,drop=FALSE],aes(x=PCoA1,y=PCoA2),color="gray25",
            label="Delta",size=2,face="bold",nudge_y=-0.06)+
  geom_point(data=df_site[3,1:3,drop=FALSE],aes(x=PCoA1,y=PCoA2),color="gray24",size=1.5)+
  geom_text(data=df_site[3,1:3,drop=FALSE],aes(x=PCoA1,y=PCoA2),color="gray25",
            label="River.A",size=2,face="bold",nudge_y=-0.07)+
  geom_point(data=df_site[4,1:3,drop=FALSE],aes(x=PCoA1,y=PCoA2),color="gray24",size=1.5)+
  geom_text(data=df_site[4,1:3,drop=FALSE],aes(x=PCoA1,y=PCoA2),color="gray25",
            label="River.B",size=2,face="bold",nudge_y=-0.07)+
  
  #Plot species abundance 
  geom_point(data=df_species[1,1:2,drop=FALSE],
             aes(x=RDA1,y=RDA2, shape="Pseudomonas.sp.J4AJ"),size=1.2)+
  geom_point(data=df_species[2,1:2,drop=FALSE],
             aes(x=RDA1,y=RDA2, shape="Kocuria.sp.97H2c"),size=1.2)+
  geom_point(data=df_species[3,1:2,drop=FALSE],
             aes(x=RDA1,y=RDA2, shape="Stenotrophomonas.sp.Hy3tC5"),size=1.2)+
  geom_point(data=df_species[4,1:2,drop=FALSE],
             aes(x=RDA1,y=RDA2, shape="Pseudomonas.alcaligenes.B2M2O"),size=1.2)+
  geom_point(data=df_species[5,1:2,drop=FALSE],
             aes(x=RDA1,y=RDA2, shape="Acinetobacter.sp.A17"),size=1.2)+
  geom_point(data=df_species[6,1:2,drop=FALSE],
             aes(x=RDA1,y=RDA2, shape="Ponticoccus.gilvus.19DR"),size=1.2)+
  geom_point(data=df_species[7,1:2,drop=FALSE],
             aes(x=RDA1,y=RDA2, shape= "Agrobacterium.tumefaciens"),size=1.2)+
  geom_point(data=df_species[8,1:2,drop=FALSE],
             aes(x=RDA1,y=RDA2, shape= "Delftia.sp.C5"),size=1.2)+
  geom_point(data=df_species[9,1:2,drop=FALSE],
             aes(x=RDA1,y=RDA2, shape="Pseudomonas.aeruginosa.SMVIT1"),size=1.2)+
  geom_point(data=df_species[10,1:2,drop=FALSE],
             aes(x=RDA1,y=RDA2, shape="Vibrio.sp.442"),size=1.2)+
  scale_shape_manual(name="16srRNA_Taxa",
                     values=c(0,1,2,3,4,5,6,7,8,9))+  
  #each numbers in the shape values correspond to a specific shape
  
  #Additional customization to organize legends and spacing
  theme_classic()+
  theme(
    panel.spacing.x = unit(0.5,"cm"),
    panel.background = element_rect(color="black"),
    legend.title=element_text(size=rel(0.8)),
    legend.text=element_text(size=rel(0.5),face="bold"),
    legend.spacing.y=unit(0.1,"cm"),
    legend.key.height = unit(0.2,"cm"),
    legend.background = element_rect(color="black")
  )+
  #this color palette is from ggsci package, you can make your own color combinations
  #by scale_color_manual() or find other palettes in the ggsci package
  scale_color_futurama(
    name="Water Quality Parameters"
  )+
  #Axis Label and Scaling, adjustment depends on the output 
  labs(x=paste0("Component 1 (",RDA1_varex," %)"),
       y=paste0("Component 2 (",RDA2_varex," %)"))+
  
  
  scale_x_continuous(limits=c(-3,3)) +
  scale_y_continuous(limits=c(-2.5,2.5))

