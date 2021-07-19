####Final Simulation Step 3: make plots
#March 10, 2021
#Guanghao Zhang
#ghzhang@umich.edu

rm(list=ls())

library(plotrix)
library(ggplot2)
library(reshape)
library(ggpubr)
library(gg.gap)
library(patchwork)
library(gridExtra)

#########################RE: np=200, ne close to ne'
load("sim_bignp_close_final_total.RData")

#data for var()/var(2)
RE_bignp_close_alt=data.frame(t(RE_bignp))
colnames(RE_bignp_close_alt) <- c("0.1","0.5","0.9")
RE_bignp_close_alt$Method=c("1","2","3a","3b","3c","3d","4")#c("2","3a","3b","3c","3d","4")
RE_bignp_close_alt <- melt(RE_bignp_close_alt, id.vars = "Method")

RE_plot_bignp_close_alt <- ggplot(data = RE_bignp_close_alt, aes(x=Method, y=value, group=variable)) +
  geom_line(aes(linetype=variable)) + geom_point(aes(shape=Method)) + theme_classic() + xlab("Approach") + ylab("Relative Efficiency") + 
  scale_linetype_discrete(name = "PVE") + 
  scale_shape_discrete(name = "Approach") +
  scale_shape_manual(values=seq(1,7)) + 
  labs(subtitle = expression(paste("A: n'"[e], " < n"[e]))) +
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))

RE_plot_bignp_close_alt

RE_plot_bignp_close_alt_1 <- gg.gap(plot=RE_plot_bignp_close_alt,
       segments=c(1.2,1.7),
       tick_width = c(0.5,1),
       ylim=c(0,3.5))
RE_plot_bignp_close_alt_1

#########################RE: np=200, ne far from ne'
load("sim_bignp_far_final_total.RData")

#data for var()/var(2)

RE_bignp_far_alt=data.frame(t(RE_bignp))
colnames(RE_bignp_far_alt) <- c("0.1","0.5","0.9")
RE_bignp_far_alt$Method=c("1","2","3a","3b","3c","3d","4")
RE_bignp_far_alt <- melt(RE_bignp_far_alt, id.vars = "Method")

RE_plot_bignp_far_alt <- ggplot(data = RE_bignp_far_alt, aes(x=Method, y=value, group=variable)) +
  geom_line(aes(linetype=variable)) + geom_point(aes(shape=Method)) + theme_classic() + xlab("Approach") + ylab("Relative Efficiency") + 
  scale_linetype_discrete(name = "PVE") + 
  scale_shape_discrete(name = "Approach") +
  scale_shape_manual(values=seq(1,7)) + 
  labs(subtitle = expression(paste("B: n'"[e], " << n"[e]))) +
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))

RE_plot_bignp_far_alt

RE_plot_bignp_far_alt_1 <- gg.gap(plot=RE_plot_bignp_far_alt,
                          segments=c(1.1,1.6),
                          tick_width = c(0.5,3),
                          ylim=c(0,6))
RE_plot_bignp_far_alt_1

########################################################################################################################
#########################RE: np=50, ne close to ne'
load("sim_smallnp_close_final_total.RData")

RE_smallnp_close_alt=data.frame(t(RE_smallnp))
colnames(RE_smallnp_close_alt) <- c("0.1","0.5","0.9")
RE_smallnp_close_alt$Method=c("1","2","3a","3b","3c","3d","4")
RE_smallnp_close_alt <- melt(RE_smallnp_close_alt, id.vars = "Method")

RE_plot_smallnp_close_alt <- ggplot(data = RE_smallnp_close_alt, aes(x=Method, y=value, group=variable)) +
  geom_line(aes(linetype=variable)) + geom_point(aes(shape=Method)) + theme_classic() + xlab("Approach") + ylab("Relative Efficiency") + 
  scale_linetype_discrete(name = "PVE") + 
  scale_shape_discrete(name = "Approach") +
  scale_shape_manual(values=seq(1,7)) + 
  labs(subtitle = expression(paste("A: n'"[e], " < n"[e]))) +
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))

RE_plot_smallnp_close_alt

#########################RE: np=50, ne far from ne'
load("sim_smallnp_far_final_total.RData")

RE_smallnp_far_alt=data.frame(t(RE_smallnp))
colnames(RE_smallnp_far_alt) <- c("0.1","0.5","0.9")
RE_smallnp_far_alt$Method=c("1", "2","3a","3b","3c","3d","4")
RE_smallnp_far_alt <- melt(RE_smallnp_far_alt, id.vars = "Method")

RE_plot_smallnp_far_alt <- ggplot(data = RE_smallnp_far_alt, aes(x=Method, y=value, group=variable)) +
  geom_line(aes(linetype=variable)) + geom_point(aes(shape=Method)) + theme_classic() + xlab("Approach") + ylab("Relative Efficiency") + 
  scale_linetype_discrete(name = "PVE") + 
  scale_shape_discrete(name = "Approach") +
  scale_shape_manual(values=seq(1,7)) + 
  labs(subtitle = expression(paste("B: n'"[e], " << n"[e]))) +
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))

RE_plot_smallnp_far_alt

###########################combine subfigs
sum_plot1 <- ggarrange(RE_plot_bignp_close_alt_1, RE_plot_bignp_far_alt_1,#beta_plot_bignp_close,beta_plot_bignp_far,
                           ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")
sum_plot1
ggsave("sim_fig1.pdf", width = 9, height = 4.5, sum_plot1)

sum_plot2 <- ggarrange(RE_plot_smallnp_close_alt, RE_plot_smallnp_far_alt,#beta_plot_smallnp_close,beta_plot_smallnp_far,
                       ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")
sum_plot2
ggsave("sim_fig2.pdf", width = 9, height = 4.5, sum_plot2)
