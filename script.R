#Scenarios
require(tidyr)
require(dplyr)
require(ggplot2)
require(cowplot)
require(viridis)

source('model.R')

sc1 <- list(beta=0.15,
            delta=4,
            Ea=0.52, 
            gamma=40, 
            h=20, 
            k=1.2,
            lO=1, 
            m=100, 
            M=2 ,
            n=0.88, 
            O2crit=2, 
            omega=0.03, 
            p=0.8 ,
            P50=4.2, 
            phi=0.25, 
            q=0.75 ,
            shape=0.1,
            temp=seq(8,26,l=100),
            Topt=20, 
            Tref=15 ,
            v=1)

sc2 <- list(beta=0.15,
            delta=4,
            Ea=0.52, 
            gamma=40, 
            h=20, 
            k=1.2,
            lO=1, 
            m=100, 
            M=2 ,
            n=0.88, 
            O2crit=2, 
            omega=0.03, 
            p=0.8 ,
            P50=4.2, 
            phi=0.25, 
            q=0.75 ,
            shape=1,
            temp=seq(8,26,l=100),
            Topt=20, 
            Tref=15 ,
            v=1)

sc3 <- list(beta=0.15,
            delta=4,
            Ea=0.52, 
            gamma=40, 
            h=20, 
            k=1.2,
            lO=1, 
            m=100,  
            M=2 ,
            n=0.88, 
            O2crit=2, 
            omega=0.03, 
            p=0.8 ,
            P50=4.2, 
            phi=0.25, 
            q=0.75 ,
            shape=5,
            temp=seq(8,26,l=100),
            Topt=20, 
            Tref=15 ,
            v=1)

scenarios <- bind_cols(data_frame(Scenario = as.factor(rep(1:3,each=100))),
                       bind_rows(plot_data_temp(sc1),
                       plot_data_temp(sc2),
                       plot_data_temp(sc3)))

ggp <-lst()

sc_act <- tidyr::gather(scenarios,Activity,value,Optimum, Limit,Realised) %>%
  mutate(Activity = relevel(as.factor(Activity), "Realised"))

ggp[['Activity']] <- ggplot(sc_act) + 
  geom_line(aes(x=Temperature, y=value, linetype=Activity, col=Scenario),alpha=0.7) +
  geom_point(data=sc_act[seq.int(1,nrow(sc_act),by = 5),],aes(x=Temperature, y=value, group=Activity, col=Scenario,shape=Scenario),alpha=0.5,size=2) +
  theme_cowplot()+
  #scale_size_discrete(guide='none',range = c(1.2,0.6))+
  theme(legend.justification=c(0,1), 
        legend.position=c(0.55,1),
        legend.box = "horizontal",
        legend.text  = element_text(size=8),
        legend.title = element_text(size=10),
        legend.spacing = unit(0.005,'npc'),
        legend.background = element_rect(fill='white')) + 
  labs(x = expression("Mean temperature " ( degree*C)),
       y=expression('Activity ' (tau)))

sco <- tidyr::gather(scenarios,Rate,value,MMR, `Active M.`, `Std M.`)

ggp[['Oxygen']] <- sco %>%
  ggplot() + 
  geom_line(aes(x=Temperature, y=value, linetype=Rate, col=Scenario),alpha=0.7) +
  geom_point(data=sco[seq.int(1,nrow(sco),by = 5),],aes(x=Temperature, y=value, group=Rate, col=Scenario,shape=Scenario),alpha=0.5,size=2) +
  theme_cowplot()+
  theme(legend.justification=c(0,1), 
        legend.position=c(0.05,0.85),
        legend.text  = element_text(size=8),
        legend.title = element_text(size=10),
        legend.spacing = unit(0.005,'npc')#legend.background = element_rect(fill='white')
        )+
  scale_color_discrete(guide='none')+
  scale_shape_discrete(guide='none')+
  labs(x = expression("Mean temperature " ( degree*C)),
       y=expression(Oxygen~(mg.yr^-1))) 

sc <- function(x) x/max(x)
sc_feed <- scenarios %>% 
  mutate(Consumption=sc(Consumption),
         Efficiency=sc(Efficiency),
         `Available Energy` = sc(`C for growth`),
         `Available Energy` = ifelse(`Available Energy`<0,0,`Available Energy`)) %>% 
  tidyr::gather(Rate,value,`Feeding level`,`Available Energy`,M) %>%
  mutate(Rate = relevel(as.factor(Rate),'Feeding level'))

ggp[['Feeding']] <- sc_feed  %>% 
  ggplot() + 
  geom_line(aes(x=Temperature, y=value, linetype=Rate, col=Scenario),alpha=0.7) +
  geom_point(data=sc_feed[seq.int(1,nrow(sc_feed),by = 5),],aes(x=Temperature, y=value, group=Rate, col=Scenario,shape=Scenario),alpha=0.5,size=2) +
  theme_cowplot()+
  lims(y=c(0,1))+
  theme(legend.justification=c(1,1), 
        legend.position=c(0.8,0.8),
        legend.text  = element_text(size=8),
        legend.title = element_text(size=10),
        legend.spacing = unit(0.005,'npc'))+
  scale_color_discrete(guide='none')+
  scale_shape_discrete(guide='none')+
  labs(x = expression("Mean temperature " ( degree*C)),
       y='Relative rate') 

lm=10^seq(0,6,l=1000)
n_int = 50
growth_scenarios <- bind_cols(data_frame(Scenario = as.factor(rep(1:3,each=n_int))),
                       bind_rows(plot_data_growth_tO2(sc1, lm = lm, n_int = n_int),
                                 plot_data_growth_tO2(sc2, lm = lm, n_int = n_int),
                                 plot_data_growth_tO2(sc3, lm = lm, n_int = n_int)))


lw <- function(w) (w/(10^-2))^(1/3)
growth_scenarios <- growth_scenarios %>%
  mutate(length=lw(Temp))

ggp[['winf']] <- ggplot(growth_scenarios) + 
  geom_line(aes(x=Temperature, y=Temp, col=Scenario),alpha=0.7) +
  geom_point(data=growth_scenarios[seq.int(1,nrow(growth_scenarios),by = 5),],aes(x=Temperature, y=Temp, col=Scenario,shape=Scenario),alpha=0.5,size=2) +
  theme_cowplot()+
  theme(legend.justification=c(0,1), 
        legend.position=c(0.4,0.8),
        legend.text  = element_text(size=8),
        legend.title = element_text(size=10),
        legend.spacing = unit(0.005,'npc'),
        legend.background = element_rect(fill='white'))+
  scale_color_discrete(guide='none')+
  scale_shape_discrete(guide='none')+
  scale_y_log10(limits=c(1,1000))+
  labs(x = expression("Mean temperature " ( degree*C)),
       y='Adult weight (g)')

cowplot::plot_grid(plotlist = ggp,labels = 'auto')
ggsave('activity.pdf',width = 8.5,height = 8.5)

#################################################
############# scale     #######################
#################################################

winf_sc1 <- sc1
winf_sc2 <- sc2
winf_sc3 <- sc3

vsc1 <- 0.5
vsc2 <- 2

winf_sc1.1 <- winf_sc1
winf_sc1.1$v <- vsc1
winf_sc1.2 <- winf_sc1
winf_sc1.2$v <- vsc2
winf_sc2.1 <- winf_sc2
winf_sc2.1$v <- vsc1
winf_sc2.2 <- winf_sc2
winf_sc2.2$v <- vsc2
winf_sc3.1 <- winf_sc3
winf_sc3.1$v <- vsc1
winf_sc3.2 <- winf_sc3
winf_sc3.2$v <- vsc2

lm=10^seq(0,6,l=1000)
n_int = 200

winf_scenarios_shape <- bind_cols(data_frame(Scenario = as.factor(rep(1:3,each=n_int*3)),
                                  Shape = as.factor(rep(rep(c(1,1/2,2),each=n_int),3))),
                              bind_rows(plot_data_growth_tO2(winf_sc1,lm,n_int = n_int),
                                        plot_data_growth_tO2(winf_sc1.1,lm,n_int = n_int),
                                        plot_data_growth_tO2(winf_sc1.2,lm,n_int = n_int),
                                        plot_data_growth_tO2(winf_sc2,lm,n_int = n_int),
                                        plot_data_growth_tO2(winf_sc2.1,lm,n_int = n_int),
                                        plot_data_growth_tO2(winf_sc2.2,lm,n_int = n_int),
                                        plot_data_growth_tO2(winf_sc3,lm,n_int = n_int),
                                        plot_data_growth_tO2(winf_sc3.1,lm,n_int = n_int),
                                        plot_data_growth_tO2(winf_sc3.2,lm,n_int = n_int)))

ggw <- lst()
ggw[['winf_shape']] <- winf_scenarios_shape %>%
  mutate(Shape = relevel(as.factor(Shape),'1')) %>% 
  ggplot() + 
  geom_line(aes(x=Temperature, 
                y=Temp, 
                col=Scenario,
                linetype=Shape), 
            alpha=1,
            position=position_dodge(width=1)) +
  theme_cowplot()+
  scale_y_log10(limits=c(1,1000))+
  theme(legend.justification=c(0,1), 
        legend.position=c(0.65,0.95),
        legend.text  = element_text(size=8),
        legend.title = element_text(size=10),
        legend.spacing = unit(0.005,'npc'),
        legend.background = element_rect(fill='white'))+
  scale_color_discrete(guide='none')+
  labs(x = expression("Mean temperature " ( degree*C)),
       y='Adult weight (g)')

ggw[['winf_shape']] 

#################################################
############# prey avail. #######################
#################################################

gammasc1 <- 30
gammasc2 <- 50

winf_sc1.1 <- winf_sc1
winf_sc1.1$gamma <- gammasc1
winf_sc1.2 <- winf_sc1
winf_sc1.2$gamma <- gammasc2
winf_sc2.1 <- winf_sc2
winf_sc2.1$gamma <- gammasc1
winf_sc2.2 <- winf_sc2
winf_sc2.2$gamma <- gammasc2
winf_sc3.1 <- winf_sc3
winf_sc3.1$gamma <- gammasc1
winf_sc3.2 <- winf_sc3
winf_sc3.2$gamma <- gammasc2

winf_scenarios_gamma <- bind_cols(data_frame(Scenario = as.factor(rep(1:3,each=n_int*3)),
                                       `Consumption rate`= as.factor(rep(rep(c(40,30,50),each=n_int),3))),
                                  bind_rows(plot_data_growth_tO2(winf_sc1,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc1.1,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc1.2,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc2,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc2.1,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc2.2,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc3,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc3.1,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc3.2,lm,n_int = n_int)))

ggw[['winf_prey']] <- winf_scenarios_gamma %>%
  mutate(`Consumption rate` = relevel(as.factor(`Consumption rate`),'40')) %>% 
  ggplot() + 
  geom_line(aes(x=Temperature, 
                y=Temp, 
                col=Scenario,
                linetype=`Consumption rate`), 
            alpha=1,
            position=position_dodge(width=1)) +
  theme_cowplot()+
  theme(legend.justification=c(0,1), 
        legend.position=c(0.65,0.95),
        legend.text  = element_text(size=8),
        legend.title = element_text(size=10),
        legend.spacing = unit(0.005,'npc'),
        legend.background = element_rect(fill='white'))+
  scale_color_discrete(guide='none')+
  scale_y_log10(limits=c(1,1000))+
  labs(x = expression("Mean temperature " ( degree*C)),
       y='Adult weight (g)')


ggw[['winf_prey']]

#################################################
############# metabolics  #######################
#################################################

deltasc1 <- 3
deltasc2 <- 5

winf_sc1.1 <- winf_sc1
winf_sc1.1$delta <- deltasc1
winf_sc1.2 <- winf_sc1
winf_sc1.2$delta <- deltasc2
winf_sc2.1 <- winf_sc2
winf_sc2.1$delta <- deltasc1
winf_sc2.2 <- winf_sc2
winf_sc2.2$delta <- deltasc2
winf_sc3.1 <- winf_sc3
winf_sc3.1$delta <- deltasc1
winf_sc3.2 <- winf_sc3
winf_sc3.2$delta <- deltasc2

winf_scenarios_delta <- bind_cols(data_frame(Scenario = as.factor(rep(1:3,each=3*n_int)),
                                             `Activity scaling` = as.factor(rep(rep(c(4,3,5),each=n_int),3))),
                                  bind_rows(plot_data_growth_tO2(winf_sc1,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc1.1,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc1.2,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc2,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc2.1,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc2.2,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc3,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc3.1,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc3.2,lm,n_int = n_int)))

ggw[['winf_delta']] <- winf_scenarios_delta %>%
  mutate(FAS = relevel(as.factor(`Activity scaling`),'4')) %>% 
  ggplot() + 
  geom_line(aes(x=Temperature, 
                y=Temp, 
                col=Scenario,
                linetype=FAS), 
            alpha=1,
            position=position_dodge(width=1)) +
  theme_cowplot()+
  scale_y_log10(limits=c(1,1000))+
  theme(legend.justification=c(0,1), 
        legend.position=c(0.65,0.95),
        legend.text  = element_text(size=8),
        legend.title = element_text(size=10),
        legend.spacing = unit(0.005,'npc'),
        legend.background = element_rect(fill='white'))+
  scale_color_discrete(guide='none')+
  labs(x = expression("Mean temperature " ( degree*C)),
       y='Adult weight (g)')


#################################################
############# exponents n  #######################
#################################################
lm=10^seq(0,6,l=1000)
nsc1 <- 0.85
nsc2 <- 0.9

winf_sc1.1 <- winf_sc1
winf_sc1.1$n <- nsc1
winf_sc1.2 <- winf_sc1
winf_sc1.2$n <- nsc2 
winf_sc2.1 <- winf_sc2
winf_sc2.1$n <- nsc1
winf_sc2.2 <- winf_sc2
winf_sc2.2$n <- nsc2 
winf_sc3.1 <- winf_sc3
winf_sc3.1$n <- nsc1
winf_sc3.2 <- winf_sc3
winf_sc3.2$n <- nsc2 

winf_scenarios_n <- bind_cols(data_frame(Scenario = as.factor(rep(1:3,each=3*n_int)),
                                             `Exponent (n)` = as.factor(rep(rep(c(0.88,0.85,0.9),each=n_int),3))),
                                  bind_rows(plot_data_growth_tO2(winf_sc1,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc1.1,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc1.2,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc2,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc2.1,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc2.2,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc3,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc3.1,lm,n_int = n_int),
                                            plot_data_growth_tO2(winf_sc3.2,lm,n_int = n_int)))

ggw[['winf_n']] <- winf_scenarios_n %>%
  mutate(`Exponent (n)` = relevel(as.factor(`Exponent (n)`),'0.88')) %>% 
  ggplot() + 
  geom_line(aes(x=Temperature, 
                y=Temp, 
                col=Scenario,
                linetype=`Exponent (n)`), 
            alpha=1,
            position=position_dodge(width=1)) +
  scale_y_log10(limits=c(1,1000))+
  theme_cowplot()+
  theme(legend.justification=c(0,1), 
        legend.position=c(0.65,0.95),
        legend.text  = element_text(size=8),
        legend.title = element_text(size=10),
        legend.spacing = unit(0.005,'npc'),
        legend.background = element_rect(fill='white'))+
  scale_color_discrete(guide='none')+
  labs(x = expression("Mean temperature " ( degree*C)),
       y='Adult weight (g)')

#################################################
############# exponents q  #######################
#################################################
lm=10^seq(0,6,l=1000)
qsc1 <- 0.66
qsc2 <- 0.8

winf_sc1.1 <- winf_sc1
winf_sc1.1$q <- qsc1
winf_sc1.2 <- winf_sc1
winf_sc1.2$q <- qsc2 
winf_sc2.1 <- winf_sc2
winf_sc2.1$q <- qsc1
winf_sc2.2 <- winf_sc2
winf_sc2.2$q <- qsc2 
winf_sc3.1 <- winf_sc3
winf_sc3.1$q <- qsc1
winf_sc3.2 <- winf_sc3
winf_sc3.2$q <- qsc2 

winf_scenarios_q <- bind_cols(data_frame(Scenario = as.factor(rep(1:3,each=3*n_int)),
                                         `Exponent (q)` = as.factor(rep(rep(c(0.75,0.66,0.8),each=n_int),3))),
                              bind_rows(plot_data_growth_tO2(winf_sc1,lm,n_int = n_int),
                                        plot_data_growth_tO2(winf_sc1.1,lm,n_int = n_int),
                                        plot_data_growth_tO2(winf_sc1.2,lm,n_int = n_int),
                                        plot_data_growth_tO2(winf_sc2,lm,n_int = n_int),
                                        plot_data_growth_tO2(winf_sc2.1,lm,n_int = n_int),
                                        plot_data_growth_tO2(winf_sc2.2,lm,n_int = n_int),
                                        plot_data_growth_tO2(winf_sc3,lm,n_int = n_int),
                                        plot_data_growth_tO2(winf_sc3.1,lm,n_int = n_int),
                                        plot_data_growth_tO2(winf_sc3.2,lm,n_int = n_int)))

ggw[['winf_q']] <- winf_scenarios_q %>%
  mutate(`Exponent (q)` = relevel(as.factor(`Exponent (q)`),'0.75')) %>% 
  ggplot() + 
  geom_line(aes(x=Temperature, 
                y=Temp, 
                col=Scenario,
                linetype=`Exponent (q)`), 
            alpha=1,
            position=position_dodge(width=1)) +
  scale_y_log10(limits=c(1,1000))+
  theme_cowplot()+
  theme(legend.justification=c(0,1), 
        legend.position=c(0.65,0.95),
        legend.text  = element_text(size=8),
        legend.title = element_text(size=10),
        legend.spacing = unit(0.005,'npc'),
        legend.background = element_rect(fill='white'))+
  scale_color_discrete(guide='none')+
  labs(x = expression("Mean temperature " ( degree*C)),
       y='Adult weight (g)')

cowplot::plot_grid(plotlist = ggw[1:4],labels = 'auto',ncol = 2,hjust = -10)
ggsave('winf.pdf',width = 8,height = 8)

save.image(file='source_data.RData')


