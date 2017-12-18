#some interesting links
#http://www.who.int/neglected_diseases/preventive_chemotherapy/sth/en/
#http://gco.iarc.fr/today/online-analysis-map?mode=cancer&mode_population=continents&population=900&sex=0&cancer=7&type=0&statistic=0&prevalence=0&color_palette=default&projection=natural-earth
#http://profiles.ucsf.edu/joe.wiemels #collaborator
#http://iicc.iarc.fr/results/registries.php
#http://iicc.iarc.fr/results/
#http://iicc.iarc.fr/results/comparative.php childhood cancers



setwd("~/Documents/stanford/hcc/data/epi")

cancers = c("bladder", "brain", "breast", "cervix", "colon", "kidney", "leukemia", "liver", "lung", "ovary", 
            "prostate", "stomach", "corpus", "gallbladder", "lyphoma", "myeloma", "pancreatic", "testis", "thyroid")
ps = NULL
cors = NULL
#uncomment to analyze all cancer data
#for (cancer in cancers){
#liver = read.csv(paste(cancer, ".csv", sep=""))

#confounding factors
poverty = read.csv("poverty.csv")
poverty = poverty[poverty$Series.Code == "NY.GDP.MKTP.CD", c("Country.Name", "X2010..YR2010.")]
poverty$country = tolower(poverty$Country.Name)
poverty$poverty = log(poverty$X2010..YR2010.)

cancer = "all cancer"
cancer_incidence = read.csv("all_cancer.csv")

parasite = read.csv("heminstic_all.csv") #Schistosomiasis_all.csv, can change to different diseases.
parasite$mean_coverage = apply(parasite[, c("Programme.coverage.PreSAC", 
                                      "National.coverage.PreSAC",
                              "Programme.coverage.SAC", "National.coverage.SAC")], 1, function(x) max(x, na.rm =T))


#Mbd Pzq  Dec Alb
Mbd_use = sapply(parasite$Drug.used.PreSAC, function(x){
  if (length(grep("Alb|Mbd", x, ignore.case = F))> 0){
    T
  }else{
    F
  }
})

#could subset by year or drug
#parasite = parasite[Mbd_use, ]
#parasite = parasite[parasite$Year == 2012, ]

#we should manually check the mis-matched countries between two datasets
cancer_incidence$country = tolower(cancer_incidence$Population)
parasite$country = tolower(parasite$Country)
cancer_incidence_parasite = merge(cancer_incidence, parasite, by = "country")
cancer_incidence_parasite = merge(cancer_incidence_parasite, poverty, by = "country")
for (i in 6:ncol(cancer_incidence_parasite)){ #7 
  x = as.numeric(as.character(cancer_incidence_parasite[, i]))
  y = cancer_incidence_parasite$Value
  if (sum(!is.na(x)) < 10) next
  print(colnames(cancer_incidence_parasite)[i])
  print(cor.test(x, y))
}

cancer_incidence_parasite_subset = unique(cancer_incidence_parasite[, c("Value", "Year", "country", "National.coverage.PreSAC", "poverty")])
cancer_incidence_parasite_subset = cancer_incidence_parasite_subset[!is.na(cancer_incidence_parasite_subset$National.coverage.PreSAC), ]
cancer_incidence_parasite_subset = aggregate(cbind(Value, National.coverage.PreSAC, poverty) ~ country, cancer_incidence_parasite_subset, mean)
x = as.numeric(as.character(cancer_incidence_parasite_subset[, "National.coverage.PreSAC"]))
y = cancer_incidence_parasite_subset$Value
plot(x, y)
test = cor.test(x, y)
test
t.test(y[x<0.5], y[x>0.5])
mean(y[x<0.5])/mean(y[x>0.5])
boxplot(y[x<0.5], y[x>0.5])
mean(y[x<0.5])/mean(y[x>0.5])

lm_model = lm(y ~ x)

lm_plot = ggplot(cancer_incidence_parasite_subset, aes( as.numeric(as.character(cancer_incidence_parasite_subset[, "National.coverage.PreSAC"])), Value  )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +                                                                                              
  stat_smooth(method="lm", se=F, color="black")  + geom_point(size=3) + 
  annotate("text", label = paste(cancer), 
           x = 0.25, y = 300, size = 6, colour = "black") +
  annotate("text", label = paste("r=", format(-1 * summary(lm_model)$r.squared ^ 0.5, digit=2), ", ",  "P=", format(anova(lm_model)$`Pr(>F)`[1], digit=2), sep=""), 
           x = 0.25, y = 280, size = 6, colour = "black") +
  scale_size(range = c(2, 5)) +
  xlab("National coverage PreSAC") + guides(shape=FALSE, size=FALSE) +
  ylab("Incidence") + coord_cartesian(xlim = c(0, 1), ylim=c(0, 300))
print(lm_plot)


cor.test(cancer_incidence_parasite_subset$Value, cancer_incidence_parasite_subset$poverty)
#write.csv(cancer_incidence_parasite_subset, "liver_parasite.csv")

#considering confounding factors
library(glmnet)
lm_model = glm(Value ~ 
                 National.coverage.PreSAC + poverty
                 , family = gaussian, data = cancer_incidence_parasite_subset)
#                 Number.of.SAC.targeted + Reported.number.of.SAC.treated  + Programme.coverage.SAC + 
#                 Reported.number.of.PreSAC.treated + Programme.coverage.PreSAC + 
#Number.of.PreSAC.targeted +  
summary(lm_model)

cancer_incidence_parasite_subset$coverage_bin = ifelse(cancer_incidence_parasite_subset$National.coverage.PreSAC > 0.5, 1, 0)
lm_model = glm(coverage_bin ~ 
                 Value+ poverty
               ,data = cancer_incidence_parasite_subset)
summary(lm_model)

#ps = c(ps, test$p.value)
#cors = c(cors, test$estimate)
#}
#liver_parasite_subset[order(liver_parasite_subset$National.coverage.PreSAC), ]
#write.csv(cancer_cor, "cancer_cor.csv")
#cancer_cor = data.frame(cancers[1:length(ps)], ps, cors)
#cancer_cor = cancer_cor[order(cancer_cor$ps), ]
