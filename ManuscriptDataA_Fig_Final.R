library(readr)
library(lme4)
library(car)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(gtable)
library(grid)
library(ggnewscale)
library(lmerTest) 
library(onewaytests)
library(r2glmm)
library(effects)


#Color Scheme
cbPalette <- c("#990951", "#CD90A3", "#76A8C3", "#006096")

path = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)


#Behavioural Data (Primary Study)
datafr <- read_csv("PrimaryStudyBehav.csv")
datafr$tp <- as.factor(datafr$tp)
datafr$gen = NA
datafr$gen[datafr$self > datafr$other & datafr$resp == 1] = 0
datafr$gen[datafr$self > datafr$other & datafr$resp == 0] = 1
datafr$gen[datafr$self < datafr$other & datafr$resp == 0] = 0
datafr$gen[datafr$self < datafr$other & datafr$resp == 1] = 1
x = datafr
x$cond[x$tp == 10] = 2
x$cond[x$tp == 1.5] = 1


ppl = unique(datafr$subj)
m = matrix(0, nrow = length(ppl), ncol = 13)
for (i in 1:length(ppl)) {
  m[i,1] = ppl[i]
  for (j in 1:2) {
    y = x[x$subj==ppl[i],]
    y = y[y$cond==j,]
    m[i, j+11] = mean(is.na(y$resp))
    y = y[!is.na(y$gen),]
    m[i, j+1] = mean(y$gen)
    m[i, j+3] = mean(y$rt)
    m[i, j+5] = mean(y$resp)
    m[i, j+7] = mean(y$given)
    m[i, j+9] = median(y$rt)
  }
}

colnames(m) <- c('subj', 'shortgen', 'longgen','shortrt','longrt','shortresp','longresp', 'shortgiven', 'longgiven', 'shortmedrt','longmedrt', 'shortna', 'longna')
m = data.frame(m)


datafr = datafr[!is.na(datafr$resp),]

datafr1 = matrix(0, nrow = 2*length(ppl), ncol = 5)
for (i in 1:length(ppl)) {
  for (j in 1:2) {
    datafr1[(j-1)*length(ppl) +i,1] = ppl[i]
    y = x[x$subj==ppl[i],]
    y = y[y$cond==j,]
    y = y[!is.na(y$gen),]
    datafr1[(j-1)*length(ppl) +i, 2] = mean(y$gen)
    datafr1[(j-1)*length(ppl) +i, 3] = mean(y$rt)
    datafr1[(j-1)*length(ppl) +i, 4] = mean(y$resp)
    datafr1[(j-1)*length(ppl) +i, 5] = j
  }
}
colnames(datafr1) <- c('subj', 'generosity','rt','resp', 'cond')
datafr1 = data.frame(datafr1)
datafr1$tp[datafr1$cond == 1] = 'High'
datafr1$tp[datafr1$cond == 2] = 'Low'


#Supplementary 1 Behavioural Data (Replication Study 1)
supp1A <- read_csv("Rep1ABehav.csv")
supp1B <- read_csv("Rep1BBehav.csv")
supp1B$subj = supp1B$subj + 40
supp1 = rbind(supp1A,supp1B)
supp1$gen = NA
supp1$gen[supp1$self > supp1$other & supp1$resp == 1] = 0
supp1$gen[supp1$self > supp1$other & supp1$resp == 0] = 1
supp1$gen[supp1$self < supp1$other & supp1$resp == 0] = 0
supp1$gen[supp1$self < supp1$other & supp1$resp == 1] = 1
x1 = supp1
x1$cond[x1$tp >5] = 2
x1$cond[x1$tp < 5] = 1

ppl1 = unique(supp1$subj)
msupp1 = matrix(0, nrow = length(ppl1), ncol = 13)
for (i in 1:length(ppl1)) {
  msupp1[i,1] = ppl1[i]
  for (j in 1:2) {
    y1 = x1[x1$subj==ppl1[i],]
    y1 = y1[y1$cond==j,]
    msupp1[i, j+1] = mean(y1$gen,na.rm = T)
    msupp1[i, j+3] = mean(y1$rt,na.rm = T)
    msupp1[i, j+5] = mean(y1$resp,na.rm = T)
    msupp1[i, j+7] = mean(y1$given,na.rm = T)
    msupp1[i, j+9] = median(y1$rt)
    msupp1[i, j+11] = mean(is.na(y1$resp))
  }
}

colnames(msupp1) <- c('subj', 'shortgen', 'longgen','shortrt','longrt','shortresp','longresp', 'shortgiven', 'longgiven', 'shortmedrt','longmedrt', 'shortna', 'longna')
msupp1 = data.frame(msupp1)



supp1 = supp1[!is.na(supp1$resp),]


supp11 = matrix(0, nrow = 2*length(ppl1), ncol = 5)
for (i in 1:length(ppl1)) {
  for (j in 1:2) {
    supp11[(j-1)*length(ppl1) +i,1] = ppl1[i]
    y1 = x1[x1$subj==ppl1[i],]
    y1 = y1[y1$cond==j,]
    y1 = y1[!is.na(y1$gen),]
    supp11[(j-1)*length(ppl1) +i, 2] = mean(y1$gen)
    supp11[(j-1)*length(ppl1) +i, 3] = mean(y1$rt)
    supp11[(j-1)*length(ppl1) +i, 4] = mean(y1$resp)
    supp11[(j-1)*length(ppl1) +i, 5] = j
  }
}
colnames(supp11) <- c('subj', 'generosity','rt','resp', 'cond')
supp11 = data.frame(supp11)
supp11$tp[supp11$cond == 1] = 'High'
supp11$tp[supp11$cond == 2] = 'Low'


#Supplementary 2 Behavioural Data (Replication Study 2)
supp2 <- read_csv("Rep2Behav.csv")
supp2$gen = NA
supp2$gen[supp2$self > supp2$other & supp2$resp == 1] = 0
supp2$gen[supp2$self > supp2$other & supp2$resp == 0] = 1
supp2$gen[supp2$self < supp2$other & supp2$resp == 0] = 0
supp2$gen[supp2$self < supp2$other & supp2$resp == 1] = 1



x2 = supp2
x2$cond[x2$tp == 10] = 2
x2$cond[x2$tp < 5] = 1


ppl2 = unique(supp2$subj)
msupp2 = matrix(0, nrow = length(ppl2), ncol = 13)
for (i in 1:length(ppl2)) {
  msupp2[i,1] = ppl2[i]
  for (j in 1:2) {
    y2 = x2[x2$subj==ppl2[i],]
    y2 = y2[y2$cond==j,]
    msupp2[i, j+1] = mean(y2$gen,na.rm = T)
    msupp2[i, j+3] = mean(y2$rt,na.rm = T)
    msupp2[i, j+5] = mean(y2$resp,na.rm = T)
    msupp2[i, j+7] = mean(y2$given,na.rm = T)
    msupp2[i, j+9] = median(y2$rt)
    msupp2[i, j+11] = mean(is.na(y2$resp))
  }
}

colnames(msupp2) <- c('subj', 'shortgen', 'longgen','shortrt','longrt','shortresp','longresp', 'shortgiven', 'longgiven','shortmedrt','longmedrt', 'shortna', 'longna')
msupp2 = data.frame(msupp2)



supp2 = supp2[!is.na(supp2$resp),]



supp21 = matrix(0, nrow = 2*length(ppl2), ncol = 5)
for (i in 1:length(ppl2)) {
  for (j in 1:2) {
    supp21[(j-1)*length(ppl2) +i,1] = ppl2[i]
    y1 = x2[x2$subj==ppl2[i],]
    y1 = y1[y1$cond==j,]
    y1 = y1[!is.na(y1$gen),]
    supp21[(j-1)*length(ppl2) +i, 2] = mean(y1$gen)
    supp21[(j-1)*length(ppl2) +i, 3] = mean(y1$rt)
    supp21[(j-1)*length(ppl2) +i, 4] = mean(y1$resp)
    supp21[(j-1)*length(ppl2) +i, 5] = j
  }
}
colnames(supp21) <- c('subj', 'generosity','rt','resp', 'cond')
supp21 = data.frame(supp21)
supp21$tp[supp21$cond == 1] = 'High'
supp21$tp[supp21$cond == 2] = 'Low'


#Primary
#Mean proportion of trials fail to respond per condition
mean(m$shortna)
mean(m$longna)

#Rep 1
#Mean proportion of trials fail to respond per condition
mean(msupp1$shortna)
mean(msupp1$longna)

#Rep 2
#Mean proportion of trials fail to respond per condition
mean(msupp2$shortna)
mean(msupp2$longna)

#Primary
#Manipulation Check: T-TEST between conditions of RTs
mean(m$shortrt)
mean(m$longrt)
t.test(m$shortrt-m$longrt)
sd(m$shortrt - m$longrt)/sqrt(60)

#Rep 1
#Manipulation Check: T-TEST between conditions of RTs
mean(msupp1$shortrt)
mean(msupp1$longrt)
t.test(msupp1$shortrt-msupp1$longrt)
sd(msupp1$shortrt - msupp1$longrt)/sqrt(65)

#Rep 2
#Manipulation Check: T-TEST between conditions of RTs
mean(msupp2$shortrt)
mean(msupp2$longrt)
t.test(msupp2$shortrt-msupp2$longrt)
sd(msupp2$shortrt - msupp2$longrt)/sqrt(49)

#Primary Study
#Group Effects of Time Pressure: T-TEST between conditions of Percentage Generous Choice
mean(m$shortgen)
mean(m$longgen)
shapiro.test(m$shortgen-m$longgen)
t.test(m$shortgen-m$longgen)
sd(m$shortgen - m$longgen)/sqrt(60)

#Rep 1
#Group Effects of Time Pressure: T-TEST between conditions of Percentage Generous Choice
mean(msupp1$shortgen)
mean(msupp1$longgen)
shapiro.test(msupp1$shortgen-msupp1$longgen)
t.test(msupp1$shortgen-msupp1$longgen)
sd(msupp1$shortgen - msupp1$longgen)/sqrt(65)

#Rep 2
#Group Effects of Time Pressure: T-TEST between conditions of Percentage Generous Choice
mean(msupp2$shortgen)
mean(msupp2$longgen)
shapiro.test(msupp2$shortgen-msupp2$longgen)
t.test(msupp2$shortgen-msupp2$longgen)
sd(msupp2$shortgen - msupp2$longgen)/sqrt(49)

#Primary Study
#Individual Differences in Effects of Time Pressure:Correlation between Generosity under Time Pressure & Change in Generosity
m$changegen = m$longgen - m$shortgen
cor.test(m$changegen, m$shortgen)
cor.test(-m$changegen, m$longgen)

#Rep 1
#Individual Differences in Effects of Time Pressure:Correlation between Generosity under Time Pressure & Change in Generosity
msupp1$changegen = msupp1$longgen - msupp1$shortgen
cor.test(msupp1$changegen, msupp1$shortgen)
cor.test(-msupp1$changegen, msupp1$longgen)

#Rep 2
#Individual Differences in Effects of Time Pressure:Correlation between Generosity under Time Pressure & Change in Generosity
msupp2$changegen = msupp2$longgen - msupp2$shortgen
cor.test(msupp2$changegen, msupp2$shortgen)
cor.test(-msupp2$changegen, msupp2$longgen)

##Fig.2a
df <- datafr1 %>%
  group_by(tp) %>%
  summarise(mean_generosity = mean(generosity),
            sd_generosity = sd(generosity),
            n_generosity = n(),
            SE_generosity = sd_generosity/sqrt(n()))
summarise_each(funs(mean, sd, se=sd(.)/sqrt(n())), generosity)

genplot1 <- ggplot(datafr1, aes(x=tp, y=generosity,fill = tp)) + 
  geom_violin(alpha = .8,width = .5, aes(fill = tp),color = 'grey40') + 
  scale_fill_manual(values=c(cbPalette[1], cbPalette[4])) +
  ggnewscale::new_scale_fill() +
  geom_crossbar(data = df, aes(x = tp, fill = tp, y = mean_generosity, 
                               ymin = mean_generosity-SE_generosity, 
                               ymax = mean_generosity + SE_generosity),
                fatten = 2, width = .6, alpha = 1)+
  geom_crossbar(data = df, aes(x = tp, fill = tp, y = mean_generosity, 
                               ymin = mean_generosity, 
                               ymax = mean_generosity),
                fatten = 2, width = .8, alpha = 1)+
  scale_fill_manual(values=c(cbPalette[2], cbPalette[3])) +
  annotate('text', x = 1.5, y = 0.9, label = '***', size = .34*8) + 
  theme(text = element_text(size = 8),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major.y = element_line(size = 0.5, linetype = 'solid',
                                          colour = "grey" ), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey"), legend.position="none",
        axis.line = element_blank(),
        axis.text=element_text(size=8))  + ylim(0,1) + 
  labs(title='',x='Time pressure', y='Proportion generosity')  + 
  geom_segment(x = 1, xend = 2, y = 0.85, yend = 0.85)

##Fig.2b
idiffplotmod <- ggplot(m, aes(x=shortgen, y=changegen)) +geom_point(color = cbPalette[1], size = .5) + geom_smooth(method='lm', color='grey40') + 
  labs(title='', x="Genererosity under high time pressure", y="Change in proportion generosity\nfrom high to low time pressure") + 
  xlim(0,1) + ylim(-.25, .3) +
  annotate('text', x = .75, y = 0.2, label = 'r = -0.313*', size = .34*8)+
  theme(text = element_text(size = 8),
        axis.line.x = element_line(color="grey", size = .5),
        axis.line.y = element_line(color="grey", size = .5),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text=element_text(size=8),
        legend.position="none") 

##Fig.2c
idiffplotext <- ggplot(m, aes(x=longgen, y=(-1*changegen))) +geom_point(color = cbPalette[4], size = .5) + geom_smooth(method='lm', color='grey40') + 
  labs(title='', x="Genererosity under low time pressure", y="Change in proportion generosity\nfrom low to high time pressure") + 
  xlim(0,1) + ylim(-.25, .3) +
  annotate('text', x = .75, y = 0.2, label = 'r = -0.101', size = .34*8)+
  theme(text = element_text(size = 8),
        axis.line.x = element_line(color="grey", size = .5),
        axis.line.y = element_line(color="grey", size = .5),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text=element_text(size=8),
        legend.position="none") 



##Fig. 2d
dfsupp1 <- supp11 %>%
  group_by(tp) %>%
  summarise(mean_generosity = mean(generosity),
            sd_generosity = sd(generosity),
            n_generosity = n(),
            SE_generosity = sd_generosity/sqrt(n()))
summarise_each(funs(mean, sd, se=sd(.)/sqrt(n())), generosity)

genplotsupp11 <- ggplot(supp11, aes(x=tp, y=generosity,fill = tp)) + 
  geom_violin(alpha = .8,width = .5, aes(fill = tp),color = 'grey40') + 
  scale_fill_manual(values=c(cbPalette[1], cbPalette[4])) +
  ggnewscale::new_scale_fill() +
  geom_crossbar(data = dfsupp1, aes(x = tp, fill = tp, y = mean_generosity, 
                                    ymin = mean_generosity-SE_generosity, 
                                    ymax = mean_generosity + SE_generosity),
                fatten = 2, width = .6, alpha = 1)+
  geom_crossbar(data = dfsupp1, aes(x = tp, fill = tp, y = mean_generosity, 
                                    ymin = mean_generosity, 
                                    ymax = mean_generosity),
                fatten = 2, width = .8, alpha = 1)+
  scale_fill_manual(values=c(cbPalette[2], cbPalette[3])) +
  annotate('text', x = 1.5, y = 0.95, label = 'n.s.', size = .34*8) + 
  theme(text = element_text(size = 8),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major.y = element_line(size = 0.5, linetype = 'solid',
                                          colour = "grey" ), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey"), legend.position="none",
        axis.line = element_blank(),
        axis.text=element_text(size=8))  + ylim(0,1) + 
  labs(title='',x='Time pressure', y='Proportion generosity')  + 
  geom_segment(x = 1, xend = 2, y = 0.85, yend = 0.85)


##Fig.2e
idiffplotmodsupp1 <- ggplot(msupp1, aes(x=shortgen, y=changegen)) +geom_point(color = cbPalette[1], size = .5) + geom_smooth(method='lm', color='grey40') + 
  labs(title='', x="Genererosity under high time pressure", y="Change in proportion generosity\nfrom high to low time pressure") + 
  xlim(0,1) + ylim(-.25, .3) +
  annotate('text', x = .75, y = 0.2, label = 'r = -0.286*', size = .34*8)+
  theme(text = element_text(size = 8),
        axis.line.x = element_line(color="grey", size = .5),
        axis.line.y = element_line(color="grey", size = .5),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text=element_text(size=8),
        legend.position="none") 

##Fig.2f
idiffplotextsupp1 <- ggplot(msupp1, aes(x=longgen, y=(-1*changegen))) +geom_point(color = cbPalette[4], size = .5) + geom_smooth(method='lm', color='grey40') + 
  labs(title='', x="Genererosity under low time pressure", y="Change in proportion generosity\nfrom low to high time pressure") + 
  xlim(0,1) + ylim(-.25, .3) +
  annotate('text', x = .75, y = 0.2, label = 'r = -0.004', size = .34*8)+
  theme(text = element_text(size = 8),
        axis.line.x = element_line(color="grey", size = .5),
        axis.line.y = element_line(color="grey", size = .5),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text=element_text(size=8),
        legend.position="none") 



##Fig.2g
dfsupp2 <- supp21 %>%
  group_by(tp) %>%
  summarise(mean_generosity = mean(generosity),
            sd_generosity = sd(generosity),
            n_generosity = n(),
            SE_generosity = sd_generosity/sqrt(n()))
summarise_each(funs(mean, sd, se=sd(.)/sqrt(n())), generosity)


genplotsupp21 <- ggplot(supp21, aes(x=tp, y=generosity,fill = tp)) + 
  geom_violin(alpha = .8,width = .5, aes(fill = tp),color = 'grey40') + 
  scale_fill_manual(values=c(cbPalette[1], cbPalette[4])) +
  ggnewscale::new_scale_fill() +
  geom_crossbar(data = dfsupp2, aes(x = tp, fill = tp, y = mean_generosity, 
                                    ymin = mean_generosity-SE_generosity, 
                                    ymax = mean_generosity + SE_generosity),
                fatten = 2, width = .6, alpha = 1)+
  geom_crossbar(data = dfsupp2, aes(x = tp, fill = tp, y = mean_generosity, 
                                    ymin = mean_generosity, 
                                    ymax = mean_generosity),
                fatten = 2, width = .8, alpha = 1)+
  scale_fill_manual(values=c(cbPalette[2], cbPalette[3])) +
  annotate('text', x = 1.5, y = 1.15, label = 'n.s.', size = .34*8) + 
  theme(text = element_text(size = 8),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major.y = element_line(size = 0.5, linetype = 'solid',
                                          colour = "grey" ), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey"), legend.position="none",
        axis.line = element_blank(),
        axis.text=element_text(size=8))  + 
  labs(title='',x='Time pressure', y='Proportion generosity')  + 
  geom_segment(x = 1, xend = 2, y = 1.05, yend = 1.05) + scale_y_continuous(limits = c(0, 1.2),breaks = c(0,.25,.5,.75,1), minor_breaks = c(.125,.375,.625,.875))



##Fig.2h
idiffplotmodsupp2 <- ggplot(msupp2, aes(x=shortgen, y=changegen)) +geom_point(color = cbPalette[1], size = .5) + geom_smooth(method='lm', color='grey40') + 
  labs(title='', x="Genererosity under high time pressure", y="Change in proportion generosity\nfrom high to low time pressure") + 
  xlim(0,1) + ylim(-.25, .3) +
  annotate('text', x = .75, y = 0.2, label = 'r = -0.300*', size = .34*8)+
  theme(text = element_text(size = 8),
        axis.line.x = element_line(color="grey", size = .5),
        axis.line.y = element_line(color="grey", size = .5),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text=element_text(size=8),
        legend.position="none") 

##Fig.2i
idiffplotextsupp2 <- ggplot(msupp2, aes(x=longgen, y=(-1*changegen))) +geom_point(color = cbPalette[4], size = .5) + geom_smooth(method='lm', color='grey40') + 
  labs(title='', x="Genererosity under low time pressure", y="Change in proportion generosity\nfrom low to high time pressure") + 
  xlim(0,1) + ylim(-.25, .3) +
  annotate('text', x = .75, y = 0.2, label = 'r = -0.101', size = .34*8)+
  theme(text = element_text(size = 8),
        axis.line.x = element_line(color="grey", size = .5),
        axis.line.y = element_line(color="grey", size = .5),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text=element_text(size=8),
        legend.position="none") 

#Fig.2
theme_set(theme_cowplot(font_size=8)) 
finalfig2 = plot_grid(genplot1,idiffplotmod,idiffplotext, 
                      genplotsupp11,idiffplotmodsupp1,idiffplotextsupp1, 
                      genplotsupp21,idiffplotmodsupp2,idiffplotextsupp2,
                      rel_widths = c(.6,1,1),
                      labels = 'auto', align = 'h', ncol = 3, nrow = 3,label_size = 8)
ggsave(paste(path,'/ATPEgrpindieffectviol.pdf',sep=''),plot=finalfig2, device = 'pdf',scale = 1,  width = 6.5, height = 6, units = "in",
       dpi = 600)

######################  EYE TRACKING  ###########################
datafr2 <- read_csv("PrimaryStudyEyetracking.csv")
earlyperiod = 286

datafr2$cond = 0
datafr2$cond[datafr2$tp > 5] = -1
datafr2$cond[datafr2$tp < 5] = 1
datafr2$ffix[datafr2$ffix == 2] = 0
datafr2[datafr2 == 'NaN'] = NA
datafr2$gen = NA
datafr2$gen[datafr2$self > datafr2$other & datafr2$resp == 1] = 0
datafr2$gen[datafr2$self > datafr2$other & datafr2$resp == 0] = 1
datafr2$gen[datafr2$self < datafr2$other & datafr2$resp == 0] = 0
datafr2$gen[datafr2$self < datafr2$other & datafr2$resp == 1] = 1



poorcalib = c(17,34,37)

for (i in poorcalib){
  datafr2 = datafr2[datafr2$subj != i,]
}


datafr2$rt[is.na(datafr2$rt)] = datafr2$tp[is.na(datafr2$rt)]
datafr2$selfgaze = NA
datafr2$selfgaze[datafr2$rt < (earlyperiod*.001)] = (datafr2$gazeself[datafr2$rt < (earlyperiod*.001)]-datafr2$gazeother[datafr2$rt < (earlyperiod*.001)])/(ceiling(datafr2$rt[datafr2$rt < (earlyperiod*.001)]*1000))
datafr2$selfgaze[datafr2$rt >= (earlyperiod*.001)] = (datafr2$gazeself[datafr2$rt >= (earlyperiod*.001)]-datafr2$gazeother[datafr2$rt >= (earlyperiod*.001)])/(earlyperiod)

restperiod = ceiling(datafr2$rt*1000)-earlyperiod
datafr2$restselfgaze = 0
datafr2$restselfgaze[restperiod > 0] = (datafr2$gazelateself[restperiod > 0] -datafr2$gazelateother[restperiod > 0])/restperiod[restperiod > 0]

summary(datafr2$ffix)
summary(datafr2$selfgaze)
summary(datafr2$restselfgaze)

datafr2$lr = (((datafr2$subj -1) %% 4)+ 1 < 3)*2 - 1

interactions::interact_plot(lm2.1, pred = 'cond', modx = 'lr')
datafr2$gazeself[datafr2$rt < (earlyperiod*.001)] = datafr2$gazeself[datafr2$rt < (earlyperiod*.001)]/(ceiling(datafr2$rt[datafr2$rt < (earlyperiod*.001)]*1000))
datafr2$gazeself[datafr2$rt >= (earlyperiod*.001)] = datafr2$gazeself[datafr2$rt >= (earlyperiod*.001)]/(earlyperiod)
datafr2$gazeother[datafr2$rt < (earlyperiod*.001)]= datafr2$gazeother[datafr2$rt < (earlyperiod*.001)]/(ceiling(datafr2$rt[datafr2$rt < .(earlyperiod*.001)]*1000))
datafr2$gazeother[datafr2$rt >= (earlyperiod*.001)]= datafr2$gazeother[datafr2$rt >= (earlyperiod*.001)]/(earlyperiod)

ppl = unique(datafr2$subj)

eyemat2 = matrix(0, nrow = length(ppl), ncol =  2)
colnames(eyemat2) <- c('subj', 'shortgen')
eyemat2 = data.frame(eyemat2)
eyemat2$subj = ppl
eyemat2$longgen = NA
for (i in 1:length(ppl)){
  eyemat2$shortgen[i] = mean(datafr2$gen[datafr2$subj == ppl[i] & datafr2$tp < 5],na.rm = T)
  eyemat2$longgen[i] = mean(datafr2$gen[datafr2$subj == ppl[i] & datafr2$tp > 5],na.rm = T)
  
  temp = datafr2[!is.na(datafr2$resp),]
  eyemat2$shortgenc[i] = sum(temp$gen[temp$subj == ppl[i] & temp$tp < 5] ==1)
  eyemat2$longgenc[i] = sum(temp$gen[temp$subj == ppl[i] & temp$tp > 5]==1)
  
  eyemat2$shortselfc[i] = sum(temp$gen[temp$subj == ppl[i] & temp$tp < 5]==0)
  eyemat2$longselfc[i] = sum(temp$gen[temp$subj == ppl[i] & temp$tp > 5]==0)
  
  eyemat2$shortgiven[i] = mean(datafr2$given[datafr2$subj == ppl[i] & datafr2$tp < 5],na.rm = T)
  eyemat2$longgiven[i] = mean(datafr2$given[datafr2$subj == ppl[i] & datafr2$tp > 5],na.rm = T)
  eyemat2$shortgaze[i] = mean(datafr2$selfgaze[datafr2$subj == ppl[i]& datafr2$tp < 5])
  eyemat2$longgaze[i] = mean(datafr2$selfgaze[datafr2$subj == ppl[i] & datafr2$tp > 5])
  
  eyemat2$restshortgaze[i] = mean(datafr2$restselfgaze[datafr2$subj == ppl[i]& datafr2$tp < 5])
  eyemat2$restlonggaze[i] = mean(datafr2$restselfgaze[datafr2$subj == ppl[i] & datafr2$tp > 5])
 }

eyemat2$lr = (((eyemat2$subj -1) %% 4)+ 1 < 3)*2 - 1
eyemat2$changegaze = eyemat2$shortgaze-eyemat2$longgaze
eyemat2$changerestgaze = eyemat2$restshortgaze-eyemat2$restlonggaze
eyemat2$changegen = eyemat2$shortgen-eyemat2$longgen
eyemat2$changegiven = eyemat2$shortgiven-eyemat2$longgiven

cor.test(eyemat2$shortgaze, eyemat2$shortgen)
cor.test(eyemat2$shortgaze, eyemat2$shortgiven)
cor.test(eyemat2$restshortgaze, eyemat2$shortgen)
cor.test(eyemat2$restshortgaze, eyemat2$shortgiven)

cor.test(eyemat2$longgaze, eyemat2$longgen)
cor.test(eyemat2$longgaze, eyemat2$longgiven)
cor.test(eyemat2$restlonggaze, eyemat2$longgen)
cor.test(eyemat2$restlonggaze, eyemat2$longgiven)


cor.test(eyemat2$lr, eyemat2$longgaze)
cor.test(eyemat2$lr, eyemat2$restlonggaze)
cor.test(eyemat2$lr, eyemat2$shortgaze)
cor.test(eyemat2$lr, eyemat2$restshortgaze)

cor.test(eyemat2$changegen , eyemat2$changegaze)
cor.test(eyemat2$changegen , eyemat2$changerestgaze)


cor.test(eyemat2$changegiven , eyemat2$changegaze)
cor.test(eyemat2$changegiven , eyemat2$changerestgaze)

cor.test(-1*eyemat2$changegaze,eyemat2$shortgaze)
cor.test(-1*eyemat2$changerestgaze,eyemat2$restshortgaze)
plot(eyemat2$shortgaze, eyemat2$changegaze)
cor.test(eyemat2$changegaze,eyemat2$longgaze)
cor.test(-1*eyemat2$changerestgaze,eyemat2$restlonggaze)
plot(eyemat2$longgaze, eyemat2$changegaze)


eyemat2$avegaze =(eyemat2$shortgaze+eyemat2$longgaze)/2
eyemat2$averestgaze =(eyemat2$restshortgaze+eyemat2$restlonggaze)/2

t.test(eyemat2$shortgaze-eyemat2$longgaze)
wilcox.test(eyemat2$changegaze,conf.int = TRUE)
t.test(eyemat2$restshortgaze-eyemat2$restlonggaze)
wilcox.test(eyemat2$changerestgaze,conf.int = TRUE)

cor.test(eyemat2$changegaze,eyemat2$shortgaze)
cor.test(eyemat2$changegaze,eyemat2$longgaze)


m3 <- matrix(0, nrow = 2*length(eyemat2$subj), ncol = 3)
colnames(m3) <- c('subj','fixbias', 'tp')
m3 <- data.frame(m3)
m3$subj = rep(eyemat2$subj,2)
m3$gazebias  = c(eyemat2$shortgaze,eyemat2$longgaze)
m3$lategazebias  = c(eyemat2$restshortgaze,eyemat2$restlonggaze)
m3$given = c(eyemat2$shortgiven, eyemat2$longgiven)
m3$genc = c(eyemat2$shortgenc, eyemat2$longgenc)
m3$selfc = c(eyemat2$shortselfc, eyemat2$longselfc)
m3$tp <- c(rep('High', length(eyemat2$subj)), rep('Low', length(eyemat2$subj)))
m3$generosity <- c(eyemat2$shortgen, eyemat2$longgen)


m3$subj <- as.factor(m3$subj)
m3$lr = ((as.numeric(levels(m3$subj))[m3$subj] -1) %% 4 + 1 < 3)*2 -1
m3$cond = c(rep(1, length(eyemat2$subj)),rep(-1,length(eyemat2$subj)))

m3$tp <- as.factor(m3$tp)
m3$tp <- relevel(as.factor(m3$tp), 'High')
summary(bf.test(gazebias ~ tp, data = m3))

summary(glm(cbind(longgenc,longselfc) ~ longgaze,eyemat2, family = 'binomial'))
summary(glm(cbind(shortgenc,shortselfc) ~ shortgaze,eyemat2, family = 'binomial'))

lm5.1 <- glmer(cbind(genc, selfc) ~ (1|subj) + gazebias*cond + lategazebias*cond, m3, family = 'binomial',
               glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
shapiro.test(residuals(lm5.1))
summary(lm5.1)
Anova(lm5.1,type =3)
x = r2beta(lm5.1, partial = T, data = m3) 
r2beta.earlygaze = x$Rsq[2]
r2beta.lategaze = x$Rsq[3]
r2beta.early.interaction = x$Rsq[4]
r2beta.cond = x$Rsq[5]


#Fig 4
ef.1 = Effect(c("gazebias", "cond"),lm5.1, SE=true,xlevels = list(gazebias = seq(-1,1, .1), cond = c(-1,1)))
df.ef1 = data.frame(ef.1)
colnames(df.ef1) <- c('gazebias', 'cond', 'generosity','se','lower','upper')
df.ef1$tp = NA
df.ef1$tp[df.ef1$cond== 1] = 'High'  
df.ef1$tp[df.ef1$cond== -1] = 'Low'
gazebycondint = ggplot(m3, aes(y = generosity, x = gazebias, color = tp)) + geom_point(size = 1) + 
  #geom_line(aes(y = generosity, x = gazebias, group = subj),color = 'black') +
  scale_color_manual(values=c(cbPalette[1], cbPalette[4])) +
  scale_fill_manual(values=c(cbPalette[2], cbPalette[3])) +
  geom_ribbon(data = df.ef1, 
              aes(x = gazebias, ymin = generosity - 1.96*se, ymax = generosity + 1.96*se,
                  fill = tp),alpha = .2,linetype = 0)+
  geom_line(data = df.ef1, aes(x =gazebias, y = generosity, color = tp), size = 1) +
  labs(title='', x="Gaze biases towards self-information", y="Proportion generosity",color= 'Time pressure',fill = 'Time pressure') + 
  theme(text = element_text(size = 8),
        axis.line.x = element_line(color="grey", size = .5),
        axis.line.y = element_line(color="grey", size = .5),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) #+ xlim(-1,1)

ggsave(paste(path,'/ATPEgazebycondint.pdf',sep=''),plot=gazebycondint,device = 'pdf',scale = 1,  width = 4, height = 3, units = "in",
       dpi = 600)

####################COMPUTATIONAL MODELS#########################

datafr4 <- read_csv("migrateMCMCEADDMposterior.csv")
excl <- c(17, 34, 37)
for (i in excl){
  datafr4 = datafr4[datafr4$subj != i,]
}
datafr4$r.rt = log(datafr4$r.rt)
datafr4$avert = log(datafr4$avert)
datafr4$comp.rt = log(datafr4$comp.rt)
subjects = unique(datafr4$subj)

datafr3 <- matrix(0, nrow = length(subjects), ncol =19)
colnames(datafr3) <- c('subj','yes', 'simyes','gen','simgen','rt','simrt','lyes','lsimyes','syes','ssimyes','lgen','lsimgen','sgen','ssimgen','lrt','lsimrt','srt','ssimrt')
datafr3 <- data.frame(datafr3)

datafr3$subj = subjects

for (s in 1:length(subjects)){
  datafr3$yes[s] = mean(datafr4$r.accept[datafr4$subj == subjects[s]],na.rm = T)
  datafr3$simyes[s] = mean(datafr4$accept[datafr4$subj == subjects[s]],na.rm = T)
  datafr3$gen[s] = mean(datafr4$r.gen[datafr4$subj == subjects[s]],na.rm = T)
  datafr3$simgen[s] = mean(datafr4$gen[datafr4$subj == subjects[s]],na.rm = T)
  datafr3$rt[s] = mean(datafr4$r.rt[datafr4$subj == subjects[s]],na.rm = T)
  datafr3$simrt[s] = mean(datafr4$avert[datafr4$subj == subjects[s]],na.rm = T)
  
  datafr3$lyes[s] = mean(datafr4$r.accept[datafr4$subj == subjects[s] & datafr4$time > 5],na.rm = T)
  datafr3$lsimyes[s] = mean(datafr4$accept[datafr4$subj == subjects[s]& datafr4$time > 5],na.rm = T)
  datafr3$lgen[s] = mean(datafr4$r.gen[datafr4$subj == subjects[s]& datafr4$time > 5],na.rm = T)
  datafr3$lsimgen[s] = mean(datafr4$gen[datafr4$subj == subjects[s]& datafr4$time > 5],na.rm = T)
  datafr3$lrt[s] = mean(datafr4$r.rt[datafr4$subj == subjects[s]& datafr4$time > 5],na.rm = T)
  datafr3$lsimrt[s] = mean(datafr4$avert[datafr4$subj == subjects[s]& datafr4$time > 5],na.rm = T)
  
  datafr3$syes[s] = mean(datafr4$r.accept[datafr4$subj == subjects[s] & datafr4$time < 5],na.rm = T)
  datafr3$ssimyes[s] = mean(datafr4$accept[datafr4$subj == subjects[s]& datafr4$time < 5],na.rm = T)
  datafr3$sgen[s] = mean(datafr4$r.gen[datafr4$subj == subjects[s]& datafr4$time < 5],na.rm = T)
  datafr3$ssimgen[s] = mean(datafr4$gen[datafr4$subj == subjects[s]& datafr4$time < 5],na.rm = T)
  datafr3$srt[s] = mean(datafr4$r.rt[datafr4$subj == subjects[s]& datafr4$time < 5],na.rm = T)
  datafr3$ssimrt[s] = mean(datafr4$avert[datafr4$subj == subjects[s]& datafr4$time < 5],na.rm = T)
}

genmat = matrix(0, nrow = length(datafr3$subj)*2, ncol = 4)
colnames(genmat) = c('subj', 'cond', 'simgen', 'gen')
genmat = data.frame(genmat)
genmat$subj = rep(datafr3$subj,2)
genmat$cond[1:length(datafr3$subj)] = 'High'
genmat$cond[(length(datafr3$subj)+1):(length(datafr3$subj)*2)] = 'Low'
genmat$simgen = c(datafr3$ssimgen, datafr3$lsimgen)
genmat$gen = c(datafr3$sgen, datafr3$lgen)
genmat$simyes = c(datafr3$ssimyes, datafr3$lsimyes)
genmat$yes = c(datafr3$syes, datafr3$lyes)
genmat$simrt = c(datafr3$ssimrt, datafr3$lsimrt)
genmat$rt = c(datafr3$srt, datafr3$lrt)

###Pearson's r between observed change in generosity and model-predicted change in  generosity
cor.test(datafr3$ssimgen-datafr3$lsimgen, datafr3$sgen-datafr3$lgen)
###Paired T between observed change in generosity and model-predicted change in generosity
t.test(datafr3$ssimgen-datafr3$lsimgen-(datafr3$sgen-datafr3$lgen))

####Figure 5a: plot of observed change in generosity and model-predicted change in generosity
changegenindi = ggplot(datafr3, aes(x = ssimgen-lsimgen, y = sgen-lgen)) + geom_point(size = 1)  + 
  xlim(-.3, .3) + ylim(-.3, .3)+ geom_abline(slope = 1, intercept = 0,linetype=2, color = 'grey40') + 
  labs(x="Simulated change in generosity", y="Observed change in generosity", title = '') + 
  theme(axis.line.x = element_line(color="grey", size = .5),
        axis.line.y = element_line(color="grey", size = .5),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 8),
        axis.text = element_text(size= 8)) 



###Pearson's r between observed change in rt and model-predicted change in t
cor.test(datafr3$ssimrt-datafr3$lsimrt, datafr3$srt-datafr3$lrt)
###Paired T between observed change in rt and model-predicted change in rt
t.test(datafr3$ssimrt-datafr3$lsimrt-(datafr3$srt-datafr3$lrt))

####Figure 5b:plot of observed change in rt and model-predicted change in rt
changeRTindi = ggplot(datafr3, aes(x = ssimrt-lsimrt, y = srt-lrt)) + 
  geom_point(size = 1) + geom_abline(slope = 1, intercept = 0,linetype=2, color = 'grey40') +
  ylim(-2.5, -.5) + xlim(-2.5,-.5)+ 
  labs(x="Simulated change in logRT (s)", y="Observed change in logRT (s)", title = '') + 
  theme(axis.line.x = element_line(color="grey", size = .5),
        axis.line.y = element_line(color="grey", size = .5),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 8),
        axis.text = element_text(size= 8)) 
cor.test(datafr3$ssimrt-datafr3$lsimrt,(datafr3$srt-datafr3$lrt))
t.test(datafr3$ssimrt-datafr3$lsimrt-(datafr3$srt-datafr3$lrt))


#Test of whether Model captures Intraindividual Differences
mat = matrix(0, nrow = 40*length(unique(datafr4$subj)),6)
mat = data.frame(mat)
colnames(mat) <- c('subj', 'time','simreal', 'prob', 'rt', 'quant') 
##Quantising Oberved Data into 10 Quantiles according to subject-specific Model Predicted Quantiles for each subject by Condition

for (i in 1:length(subjects)) {
  temp = datafr4[datafr4$subj == subjects[i],]
  indexstart = ((i-1)*40)
  mat$subj[(indexstart+1):(i*40)] = subjects[i]
  temp1 = temp[temp$time < 5,]
  temp2 = temp[temp$time > 5,]
  quantilesHTPacc = quantile(temp1$accept, seq(0,1,.1),na.rm = T)
  quantilesLTPacc = quantile(temp2$accept, seq(0,1,.1),na.rm = T)
  quantilesHTPrt = quantile(temp1$avert, seq(0,1,.1),na.rm = T)
  quantilesLTPrt = quantile(temp2$avert, seq(0,1,.1),na.rm = T)
  for (j in 1:9) {
    for (k in 1:2) {
      mat$time[indexstart+(k-1)*10+j] = 1.5
      mat$time[indexstart+20+(k-1)*10+j] = 10
      mat$quant[indexstart+(k-1)*10+j] = j
      mat$quant[indexstart+20+(k-1)*10+j] = j
      if (k == 1) {
        mat$simreal[indexstart+(k-1)*10+j] = 1
        mat$simreal[indexstart+20+(k-1)*10+j] = 1
        mat$prob[indexstart+(k-1)*10+j] = mean(temp1$accept[temp1$accept >= quantilesHTPacc[j] & temp1$accept < quantilesHTPacc[j+1]],na.rm = T)
        mat$prob[indexstart+20+(k-1)*10+j] = mean(temp2$accept[temp2$accept >= quantilesLTPacc[j] & temp2$accept < quantilesLTPacc[j+1]],na.rm = T)
        mat$rt[indexstart+(k-1)*10+j] = mean(temp1$avert[temp1$avert >= quantilesHTPrt[j] & temp1$avert < quantilesHTPrt[j+1]],na.rm = T)
        mat$rt[indexstart+20+(k-1)*10+j] = mean(temp2$avert[temp2$avert >= quantilesLTPrt[j] & temp2$avert < quantilesLTPrt[j+1]],na.rm = T)}
      else {
        mat$simreal[indexstart+(k-1)*10+j] = 0
        mat$simreal[indexstart+20+(k-1)*10+j] = 0
        mat$prob[indexstart+(k-1)*10+j] = mean(temp1$r.accept[temp1$accept >= quantilesHTPacc[j] & temp1$accept < quantilesHTPacc[j+1]],na.rm = T)
        mat$prob[indexstart+20+(k-1)*10+j] = mean(temp2$r.accept[temp2$accept >= quantilesLTPacc[j] & temp2$accept < quantilesLTPacc[j+1]],na.rm = T)
        mat$rt[indexstart+(k-1)*10+j] = mean(temp1$r.rt[temp1$avert >= quantilesHTPrt[j] & temp1$avert < quantilesHTPrt[j+1]],na.rm = T)
        mat$rt[indexstart+20+(k-1)*10+j] = mean(temp2$r.rt[temp2$avert >= quantilesLTPrt[j] & temp2$avert < quantilesLTPrt[j+1]],na.rm = T)}
    }
  }
  j = 10
  for (k in 1:2) {
    mat$time[indexstart+(k-1)*10+j] = 1.5
    mat$time[indexstart+20+(k-1)*10+j] = 10
    mat$quant[indexstart+(k-1)*10+j] = j
    mat$quant[indexstart+20+(k-1)*10+j] = j
    if (k == 1) {
      mat$simreal[indexstart+(k-1)*10+j] = 1
      mat$simreal[indexstart+20+(k-1)*10+j] = 1
      mat$prob[indexstart+(k-1)*10+j] = mean(temp1$accept[temp1$accept >= quantilesHTPacc[j]],na.rm = T)
      mat$prob[indexstart+20+(k-1)*10+j] = mean(temp2$accept[temp2$accept >= quantilesLTPacc[j]],na.rm = T)
      mat$rt[indexstart+(k-1)*10+j] = mean(temp1$avert[temp1$avert >= quantilesHTPrt[j]],na.rm = T)
      mat$rt[indexstart+20+(k-1)*10+j] = mean(temp2$avert[temp2$avert >= quantilesLTPrt[j]],na.rm = T)}
    else {
      mat$simreal[indexstart+(k-1)*10+j] = 0
      mat$simreal[indexstart+20+(k-1)*10+j] = 0
      mat$prob[indexstart+(k-1)*10+j] = mean(temp1$r.accept[temp1$accept >= quantilesHTPacc[j]],na.rm = T)
      mat$prob[indexstart+20+(k-1)*10+j] = mean(temp2$r.accept[temp2$accept >= quantilesLTPacc[j] ],na.rm = T)
      mat$rt[indexstart+(k-1)*10+j] = mean(temp1$r.rt[temp1$avert >= quantilesHTPrt[j] ],na.rm = T)
      mat$rt[indexstart+20+(k-1)*10+j] = mean(temp2$r.rt[temp2$avert >= quantilesLTPrt[j] ],na.rm = T)}
    
  }
}

subjects = unique(mat$subj)
subcors = matrix(NA, nrow = length(subjects), ncol = 5)
colnames(subcors) = c('subj', 'htpcoracc', 'ltpcoracc', 'htpcorrt', 'ltpcorrt')
subcors = data.frame(subcors)


##Obtaining Correlation coefficients between model predicted quantile means and observed quantile means for rt and acceptance for each condition
subcors$subj = subjects
htpvarrt = ltpvarrt = vector('numeric', length = length(subjects))
for (i in 1:length(subjects)) {
  tempmat = mat[mat$subj == subjects[i], ]
  subcors$htpcoracc[i] = cor.test(tempmat$prob[tempmat$simreal==0 & tempmat$time < 5],tempmat$prob[tempmat$simreal==1 & tempmat$time < 5])$estimate
  subcors$ltpcoracc[i] = cor.test(tempmat$prob[tempmat$simreal==0 & tempmat$time > 5],tempmat$prob[tempmat$simreal==1 & tempmat$time > 5])$estimate
  subcors$htpcorrt[i] = cor.test(tempmat$rt[tempmat$simreal==0 & tempmat$time < 5],tempmat$rt[tempmat$simreal==1 & tempmat$time < 5])$estimate
  subcors$ltpcorrt[i] = cor.test(tempmat$rt[tempmat$simreal==0 & tempmat$time > 5],tempmat$rt[tempmat$simreal==1 & tempmat$time > 5])$estimate
  htpvarrt[i] = var(tempmat$rt[tempmat$simreal == 0 & tempmat$time < 5])
  ltpvarrt[i] = var(tempmat$rt[tempmat$simreal == 0 & tempmat$time > 5])
}


#One sample T-test of distribution of correlation coefficients for acceptance under high time pressure
t.test(subcors$htpcoracc)
sd(subcors$htpcoracc)/sqrt(length(subcors$subj))
#One sample T-test of distribution of correlation coefficients for acceptance under low time pressure
t.test(subcors$ltpcoracc)
sd(subcors$ltpcoracc)/sqrt(length(subcors$subj))
#One sample T-test of distribution of correlation coefficients for reaction times under high time pressure
t.test(subcors$htpcorrt)
sd(subcors$htpcorrt)/sqrt(length(subcors$subj))
#One sample T-test of distribution of correlation coefficients for reaction times under low time pressure
t.test(subcors$ltpcorrt)
sd(subcors$ltpcorrt)/sqrt(length(subcors$subj))


df31 <- mat %>%
  group_by(quant,simreal,time) %>%
  summarise(mean_acc = mean(prob,na.rm=T),
            sd_acc = sd(prob,na.rm = T),
            n_acc = n(),
            SE_acc = sd_acc/sqrt(n()))
summarise_each(funs(mean, sd, se=sd(.)/sqrt(n())), prob)



df32 <- mat %>%
  group_by(quant,simreal,time) %>%
  summarise(mean_rt = mean(rt,na.rm=T),
            sd_rt = sd(rt,na.rm = T),
            n_rt = n(),
            SE_rt = sd_rt/sqrt(n()))
summarise_each(funs(mean, sd, se=sd(.)/sqrt(n())), rt)


df31$time[df31$time > 5] = 'Low'

df31$time[df31$time < 5] = 'High'

df32$time[df32$time > 5] = 'Low'

df32$time[df32$time < 5] = 'High'



####Figure 5c: plot of intradividual differences, model predicted quantile means for acceptance against observed quantile means for acceptance for each condition
accplot_intraindi = ggplot(df31[df31$simreal == 0,], aes(x = quant, y = mean_acc, fill = time)) + geom_col(width = .8, position = position_dodge(.8))+ 
  geom_errorbar(data = df31[df31$simreal==0,],aes(ymin=mean_acc - SE_acc, ymax=mean_acc + SE_acc), position = position_dodge(.8), width=0.2,color = 'grey40') + 
  geom_line(data = df31[df31$simreal==1,], aes(x = quant, y = mean_acc, color = time),  size = .75,linetype = 2, position = position_dodge(.8)) + geom_point(data = df31[df31$simreal==1,], aes(x = quant, y = mean_acc,color = time), position = position_dodge(.8),size = 1) + 
  theme(text = element_text(size = 8),panel.background = element_rect(fill = "white",
                                                                      colour = "white",
                                                                      size = 0.5, linetype = "solid"),
        panel.grid.major.y = element_line(size = 0.25, linetype = 'solid',
                                          colour = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_line(size = 0.25, linetype = 'solid',
                                          colour = "grey"),legend.position = 'none',
        axis.line = element_blank(),
        axis.ticks = element_blank()) + labs(x = 'Simulated quantile', y = 'Acceptance rate',title = '',color = 'Time\npressure',fill = 'Time\npressure') + scale_x_continuous(breaks = seq(1,10,1)) +  scale_fill_manual(values = c(cbPalette[2],cbPalette[3])) + scale_color_manual(values = c(cbPalette[1],cbPalette[4]))


####Figure 5d: plot of intradividual differences, model predicted quantile means for reaction times against observed quantile means for reaction times for each condition
rtplot_intraindi = ggplot(df32[df32$simreal == 0,], aes(x = quant, y = mean_rt,fill = time)) + geom_col(width = .8, position = position_dodge(.8)) + 
  geom_errorbar(data = df32[df32$simreal==0,],aes(ymin=mean_rt - SE_rt, ymax=mean_rt + SE_rt), width=0.2,color = 'grey40', position = position_dodge(.8)) + 
  geom_line(data = df32[df32$simreal==1,], aes(x = quant, y = mean_rt,color = time), position = position_dodge(.8),linetype = 2, size = .75) + geom_point(data = df32[df32$simreal==1,], aes(x = quant, y = mean_rt,color = time), position = position_dodge(.8),size = 1) + 
  theme(text = element_text(size = 8),panel.background = element_rect(fill = "white",
                                                                      colour = "white",
                                                                      size = 0.5, linetype = "solid"),
        panel.grid.major.y = element_line(size = 0.25, linetype = 'solid',
                                          colour = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_line(size = 0.25, linetype = 'solid',
                                          colour = "grey"),
        axis.line = element_blank(),
        axis.ticks = element_blank()) + labs(x = 'Simulated quantile', y = 'logRT (s)',title = '',color = 'Time\npressure',fill = 'Time\npressure') + scale_x_continuous(breaks = seq(1,10,1)) +  scale_fill_manual(values = c(cbPalette[2],cbPalette[3])) + scale_color_manual(values = c(cbPalette[1],cbPalette[4]))



####Figure 5
interindifit = plot_grid(changegenindi,changeRTindi, labels = c('a','b'), align = 'h',label_size = 8, rel_widths = c(5,5))
interindifit = grid.arrange(interindifit, ncol = 2,
                            heights = c(1), 
                            widths = c(1,.1))
intraindifit = plot_grid(accplot_intraindi,rtplot_intraindi, labels = c('c','d'), align = 'h',label_size = 8, rel_widths = c(1,1.2))
mcmceaddmmodelchecks = plot_grid(interindifit,intraindifit, labels = c('','',''),label_size = 8, rel_widths = c(1,1), nrow = 2)

ggsave(paste(path,'/ATPEMCMCEADDMmodelchecks.pdf',sep=''),plot=mcmceaddmmodelchecks,device = 'pdf',scale = 1,  width = 6.5, height = 6, units = "in",
       dpi = 600)

#Comparison of gaze-informed DDM to simple DDM
datafr5 <- read_csv("migratesimpleDDMposterior.csv")
excl <- c(17, 34, 37)
for (i in excl){
  datafr5 = datafr5[datafr5$subj != i,]
}
datafr5$r.rt = log(datafr5$r.rt)
datafr5$avert = log(datafr5$avert)

ppl = unique(datafr5$subj)

modelcomp = matrix(NA, nrow = length(ppl)*4, ncol = 2)
colnames(modelcomp) <- c('subj', 'cond')
modelcomp = data.frame(modelcomp)
timelim = c(1.5,10)


modelcomp$subj = rep(ppl,4)
modelcomp$cond = c(rep(1, length(ppl)*2), rep(-1, length(ppl)*2))
modelcomp$gen = 0
modelcomp$genSE = 0
modelcomp$rt = 0
modelcomp$rtSE = 0
modelcomp$model = 0
for (i in 1:length(ppl)){
  temp = datafr4[datafr4$subj == ppl[i],]
  temp2 = datafr5[datafr5$subj == ppl[i],]
  for (j in 1:2){
    temp.2 = temp[temp$time == timelim[j],]
    temp2.2 = temp2[temp2$time == timelim[j],]
    modelcomp$model[(j-1)*(length(ppl)*2)+((i-1)*2)+1] = 1     #gazeinformed
    modelcomp$model[(j-1)*(length(ppl)*2)+((i-1)*2)+2] = -1     #simple
    modelcomp$gen[(j-1)*(length(ppl)*2)+((i-1)*2)+1] = mean(abs(temp.2$gen - temp.2$r.gen),na.rm = T)
    modelcomp$genSE[(j-1)*(length(ppl)*2)+((i-1)*2)+1] = sd(abs(temp.2$gen - temp.2$r.gen),na.rm = T)/sqrt(length(abs(temp.2$gen - temp.2$r.gen)))
    modelcomp$gen[(j-1)*(length(ppl)*2)+((i-1)*2)+2] = mean(abs(temp2.2$gen - temp2.2$r.gen),na.rm = T)
    modelcomp$genSE[(j-1)*(length(ppl)*2)+((i-1)*2)+2] = sd(abs(temp2.2$gen - temp2.2$r.gen),na.rm = T)/sqrt(length(abs(temp2.2$gen - temp2.2$r.gen)))
    modelcomp$rt[(j-1)*(length(ppl)*2)+((i-1)*2)+1] = mean(abs(temp.2$avert - temp.2$r.rt),na.rm = T)
    modelcomp$rtSE[(j-1)*(length(ppl)*2)+((i-1)*2)+1] = sd(abs(temp.2$avert - temp.2$r.rt),na.rm = T)/sqrt(length(abs(temp.2$avert - temp.2$r.rt)))
    modelcomp$rt[(j-1)*(length(ppl)*2)+((i-1)*2)+2] = mean(abs(temp2.2$avert - temp2.2$r.rt),na.rm = T)
    modelcomp$rtSE[(j-1)*(length(ppl)*2)+((i-1)*2)+2] = sd(abs(temp2.2$avert - temp2.2$r.rt),na.rm = T)/sqrt(length(abs(temp2.2$avert - temp2.2$r.rt)))
    
  }
}

modelcomp$cond1 = modelcomp$cond
genAE = lmer(gen ~ (1|subj) + model*cond, modelcomp)
summary(genAE)
anova(genAE, Type = '3', test.statistic = 'F')

interactions::interact_plot(genAE, pred = 'cond', modx = 'model', plot.points = TRUE)

rtAE = lmer(rt ~ (1|subj) + model*cond, modelcomp)
summary(rtAE)
Anova(rtAE, Type = 3) 


interactions::interact_plot(rtAE, pred = 'cond', modx = 'model', plot.points = TRUE)


AEcomp = matrix(0, nrow = length(ppl)*2, ncol = 2)
colnames(AEcomp) <- c('subj', 'cond')
AEcomp <- data.frame(AEcomp)

AEcomp$subj = rep(ppl,2)
AEcomp$cond = c(rep('High', length(ppl)),rep('Low', length(ppl)))
AEcomp$gengaze = c(modelcomp$gen[modelcomp$model==1])
AEcomp$genSEgaze = c(modelcomp$genSE[modelcomp$model==1])
AEcomp$gensimple = c(modelcomp$gen[modelcomp$model==-1])
AEcomp$genSEsimple = c(modelcomp$genSE[modelcomp$model==-1])
AEcomp$rtgaze = c(modelcomp$rt[modelcomp$model==1])
AEcomp$rtSEgaze = c(modelcomp$rtSE[modelcomp$model==1])
AEcomp$rtsimple = c(modelcomp$rt[modelcomp$model==-1])
AEcomp$rtSEsimple = c(modelcomp$rtSE[modelcomp$model==-1])


#Figure 6
##Figure 6A
meanAEgen = ggplot(AEcomp, aes(x = gengaze, y = gensimple,color = cond)) + 
  geom_point(alpha = .6) + geom_abline(intercept = 0, slope = 1) + 
  geom_errorbarh(aes(xmin=gengaze - genSEgaze, xmax=gengaze + genSEgaze),alpha = .6) + 
  geom_errorbar(aes(ymin=gensimple - genSEsimple, ymax=gensimple + genSEsimple),alpha = .6) + 
  scale_color_manual(values = c(cbPalette[1],cbPalette[4])) + 
  theme(text = element_text(size = 8),panel.background = element_rect(fill = "white",
                                                                      colour = "white",
                                                                      size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.position = 'none') + labs(x = 'Mean absolute error of\ngaze-informed ADDM generosity simulations', 
                                                                            y = 'Mean absolute error of\nsimple DDM generosity simulations',title = '',
                                                                            color = 'Time\npressure') + xlim(0,.67) +ylim(0,.67)


##Figure 6B
meanAErt = ggplot(AEcomp, aes(x = rtgaze, y = rtsimple,color = cond))  + 
  geom_point(alpha = .6) + geom_abline(intercept = 0, slope = 1) + 
  geom_errorbarh(aes(xmin=rtgaze - rtSEgaze, xmax=rtgaze + rtSEgaze),alpha = .6) + 
  geom_errorbar(aes(ymin=rtsimple - rtSEsimple, ymax=rtsimple + rtSEsimple),alpha = .6) + 
  scale_color_manual(values = c(cbPalette[1],cbPalette[4])) + 
  theme(text = element_text(size = 8),panel.background = element_rect(fill = "white",
                                                                      colour = "white",
                                                                      size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + labs(x = 'Mean absolute error of\ngaze-informed ADDM logRT simulations', 
                                                   y = 'Mean absolute error of\nsimple DDM logRT simulations',title = '',
                                                   color = 'Time\npressure') + xlim(0, 1.1) + ylim(0, 1.1)

meanAEplots = plot_grid(meanAEgen, meanAErt, labels = "auto",label_size = 8, rel_widths = c(1,1.2), nrow = 1)

ggsave(paste(path,'/ATPEgazevsimpleDDM.pdf',sep=''),plot=meanAEplots,device = 'pdf',scale = 1,  width = 6.5, height = 3, units = "in",
       dpi = 600)

#Analysis of Model Parameters
mcmcparams <-read_csv("migrateMCMCEADDMparams.csv")
excl = c(17,34,37)
for (i in excl){
  mcmcparams = mcmcparams[mcmcparams$subj != i, ]
}

mcmcparams$avegen = mcmcparams$changegen =  NA
for (j in mcmcparams$subj){
  mcmcparams$avegen[mcmcparams$subj==j] = (m$shortgen[m$subj==j]+m$longgen[m$subj==j])/2
  mcmcparams$changegen[mcmcparams$subj==j] = (m$shortgen[m$subj==j]-m$longgen[m$subj==j])
}

#Paired T-Test of wself 
t.test(mcmcparams$dwself)
sd(mcmcparams$dwself)/sqrt(length(mcmcparams$dwself))
t.test(mcmcparams$dwself)$p.value*8
#Paired T-Test of wother
t.test(mcmcparams$dwother)
sd(mcmcparams$dwother)/sqrt(length(mcmcparams$dwother))
t.test(mcmcparams$dwother)$p.value*8
#Paired T-Test of wfair
t.test(mcmcparams$dwfair)
sd(mcmcparams$dwfair)/sqrt(length(mcmcparams$dwfair))
t.test(mcmcparams$dwfair)$p.value*8
#Paired T-Test of bound
t.test(mcmcparams$dbound)
sd(mcmcparams$dbound)/sqrt(length(mcmcparams$dbound))
t.test(mcmcparams$dbound)$p.value*8
#Paired T-Test of collapse
t.test(mcmcparams$dcollapse)
sd(mcmcparams$dcollapse)/sqrt(length(mcmcparams$dcollapse))
t.test(mcmcparams$dcollapse)$p.value*8
#Paired T-Test of stbias
t.test(mcmcparams$dstbias)
sd(mcmcparams$dstbias)/sqrt(length(mcmcparams$dstbias))
t.test(mcmcparams$dstbias)$p.value*8
#Paired T-Test of genbias
t.test(mcmcparams$dgenbias)
sd(mcmcparams$dgenbias)/sqrt(length(mcmcparams$dgenbias))
t.test(mcmcparams$dgenbias)$p.value*8

t.test(mcmcparams$genbias)
sd(mcmcparams$genbias)/sqrt(length(mcmcparams$genbias))
#Paired T-Test of theta
t.test(mcmcparams$dtheta)
sd(mcmcparams$dtheta)/sqrt(length(mcmcparams$dtheta))
t.test(mcmcparams$dtheta)$p.value*8


mean(mcmcparams$wself - mcmcparams$dwself/2)
sd(mcmcparams$wself - mcmcparams$dwself/2)
mean(mcmcparams$wself + mcmcparams$dwself/2)
sd(mcmcparams$wself + mcmcparams$dwself/2)

mean(mcmcparams$wother - mcmcparams$dwother/2)
sd(mcmcparams$wother- mcmcparams$dwother/2)
mean(mcmcparams$wother + mcmcparams$dwother/2)
sd(mcmcparams$wother + mcmcparams$dwother/2)

mean(mcmcparams$wfair - mcmcparams$dwfair/2)
sd(mcmcparams$wfair- mcmcparams$dwfair/2)
mean(mcmcparams$wfair + mcmcparams$dwfair/2)
sd(mcmcparams$wfair + mcmcparams$dwfair/2)

mean(mcmcparams$bound - mcmcparams$dbound/2)
sd(mcmcparams$bound- mcmcparams$dbound/2)
mean(mcmcparams$bound + mcmcparams$dbound/2)
sd(mcmcparams$bound + mcmcparams$dbound/2)

mean(mcmcparams$collapse - mcmcparams$dcollapse/2)
sd(mcmcparams$collapse- mcmcparams$dcollapse/2)
mean(mcmcparams$collapse + mcmcparams$dcollapse/2)
sd(mcmcparams$collapse + mcmcparams$dcollapse/2)

mean(mcmcparams$stbias - mcmcparams$dstbias/2)
sd(mcmcparams$stbias- mcmcparams$dstbias/2)
mean(mcmcparams$stbias + mcmcparams$dstbias/2)
sd(mcmcparams$stbias + mcmcparams$dstbias/2)

mean(mcmcparams$genbias - mcmcparams$dgenbias/2)
sd(mcmcparams$genbias- mcmcparams$dgenbias/2)
mean(mcmcparams$genbias + mcmcparams$dgenbias/2)
sd(mcmcparams$genbias + mcmcparams$dgenbias/2)


mean(mcmcparams$theta - mcmcparams$dtheta/2)
sd(mcmcparams$theta- mcmcparams$dtheta/2)
mean(mcmcparams$theta + mcmcparams$dtheta/2)
sd(mcmcparams$theta + mcmcparams$dtheta/2)


#One-Sample T-test on genbias
t.test(mcmcparams$genbias)
sd(mcmcparams$genbias)/sqrt(length(mcmcparams$genbias))

##Combining Model Parameters with Fixation Analyses (After Exclusion, n = 50)
m2 = eyemat2
m2$mcmcwself =  m2$mcmcwother =  m2$mcmcwfair =  
  m2$mcmcgenbias =  m2$mcmcdwself =  m2$mcmcdwother = 
  m2$mcmcdwfair =  m2$mcmcdgenbias =
  m2$mcmctheta = m2$mcmcdtheta = 
  m2$mcmcbound= m2$mcmcdbound = 
  m2$mcmccollapse= m2$mcmcdcollapse = NA
for (i in 1:length(m2$subj)) {
  p = mcmcparams$subj ==  m2$subj[i]
  if(any(p)){
    m2$mcmcwself[i] = mcmcparams$wself[mcmcparams$subj ==  m2$subj[i]]
    m2$mcmcwother[i] = mcmcparams$wother[mcmcparams$subj ==  m2$subj[i]]
    m2$mcmcwfair[i] = mcmcparams$wfair[mcmcparams$subj ==  m2$subj[i]]
    m2$mcmcgenbias[i] = mcmcparams$genbias[mcmcparams$subj ==  m2$subj[i]]
    m2$mcmcdwself[i] = mcmcparams$dwself[mcmcparams$subj ==  m2$subj[i]]
    m2$mcmcdwother[i] = mcmcparams$dwother[mcmcparams$subj ==  m2$subj[i]]
    m2$mcmcdwfair[i] = mcmcparams$dwfair[mcmcparams$subj ==  m2$subj[i]]
    m2$mcmcdgenbias[i] = mcmcparams$dgenbias[mcmcparams$subj ==  m2$subj[i]]
    m2$mcmctheta[i] = mcmcparams$theta[mcmcparams$subj ==  m2$subj[i]]
    m2$mcmcdtheta[i] = mcmcparams$dtheta[mcmcparams$subj ==  m2$subj[i]]
    m2$mcmcbound[i] = mcmcparams$bound[mcmcparams$subj ==  m2$subj[i]]
    m2$mcmccollapse[i] = mcmcparams$collapse[mcmcparams$subj ==  m2$subj[i]]
    m2$mcmcdbound[i] = mcmcparams$dbound[mcmcparams$subj ==  m2$subj[i]]
    m2$mcmcdcollapse[i] = mcmcparams$dcollapse[mcmcparams$subj ==  m2$subj[i]]
  }
  else {
    m2$subj[i] = NA
  }
}
m2 =  m2[!is.na( m2$subj),]


m4 <- matrix(0, nrow = length(m2$subj)*2, ncol = 2)
colnames(m4) <- c('subj', 'tp')
m4 <- data.frame(m4)
m4$subj  = rep(m2$subj,2)
m4$tp = c(rep('High',length(m2$subj)), rep('Low',length(m2$subj)))
m4$cond = c(rep(1,length(m2$subj)), rep(-1,length(m2$subj)))
m4$gaze= c(m2$shortgaze, m2$longgaze)
m4$wself = rep(m2$mcmcwself,2)
m4$wother =  rep(m2$mcmcwother,2)


#Predicting early gaze biases from individual differences in social preference (average weights)
lm6 <- lmer(gaze ~ (1|subj)+ wself*cond + wother*cond,m4)
summary(lm6)
anova(lm6, type = '3')
r2beta(lm6)

#Predicting changes in gaze biases from individual differences in social preference
lm6.2 <- lm(changegaze ~ mcmcwself+ mcmcwother, m2)
summary(lm6.2)
r2beta(lm6.2)



#Predicting changes in generosity as a function of individual differences in social preference and attention
full = lm(changegen ~ 
            mcmcwself*avegaze*changegaze+
            mcmcwother*avegaze*changegaze+
            mcmcwself*averestgaze +
            mcmcwother*averestgaze +
            mcmcwfair+
            mcmcgenbias+
            mcmcdwfair + 
            mcmcdwself +
            mcmcdwother,m2)


step.model <- step(full, direction = "both", 
                   trace = 1000,
                   steps = 100000)


Anova(full, Type=3) 

summary(step.model)
Anova(step.model, type = '3')
r2beta(step.model)
extractAIC(step.model)
extractAIC(full)
r2beta(step.model)
m5 = m2

#2-way interaction avegaze*changegaze @ wother +1SD
m5$mcmcwother1 = m5$mcmcwother -sd(m5$mcmcwother)
m5$avegaze1 = m5$avegaze + sd(m5$avegaze)
m5$changegaze1 = m5$changegaze
step.model.1 = lm(changegen ~ mcmcwself +  mcmcwother1*avegaze1*changegaze1 + mcmcwother1*averestgaze + mcmcdwother, m5)
summary1 = summary(step.model.1)
summary1$coefficients[10,]

#Simple effect of changegaze @ -1SD avegaze & wother +1SD
summary1$coefficients[5,]

#Simple effect of changegaze @ mean avegaze & wother +1SD
m5$avegaze1 = m5$avegaze 
step.model.1 = lm(changegen ~ mcmcwself +  mcmcwother1*avegaze1*changegaze1 + mcmcwother1*averestgaze + mcmcdwother, m5)
summary1 = summary(step.model.1)
summary1$coefficients[5,]

#Simple effect of changegaze @ + 1SD avegaze & wother +1SD
m5$avegaze1 = m5$avegaze - sd(m5$avegaze)
step.model.1 = lm(changegen ~ mcmcwself +  mcmcwother1*avegaze1*changegaze1 + mcmcwother1*averestgaze + mcmcdwother, m5)
summary1 = summary(step.model.1)
summary1$coefficients[5,]


#2-way interaction avegaze*changegaze @ wother mean
m5$mcmcwother1 = m5$mcmcwother 
m5$avegaze1 = m5$avegaze 
m5$changegaze1 = m5$changegaze
step.model.1 = lm(changegen ~ mcmcwself +  mcmcwother1*avegaze1*changegaze1 + mcmcwother1*averestgaze + mcmcdwother, m5)
summary1 = summary(step.model.1)
summary1$coefficients[10,]


#2-way interaction avegaze*changegaze @ wother -1SD
m5$mcmcwother1 = m5$mcmcwother + sd(m5$mcmcwother)
m5$avegaze1 = m5$avegaze 
m5$changegaze1 = m5$changegaze
step.model.1 = lm(changegen ~ mcmcwself +  mcmcwother1*avegaze1*changegaze1 + mcmcwother1*averestgaze + mcmcdwother, m5)
summary1 = summary(step.model.1)
summary1$coefficients[10,]
