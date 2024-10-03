#### Statistical Modelling for the SWI project ####

set.seed(1809081)

#### Requirements ####

source("./_Scripts/3_Descriptives.R")

ifelse("tidyverse" %in% (.packages()),
       "Tidyverse loaded from 3_Descriptives.R.", library(tidyverse))

library(brms)
library(emmeans)
library(tidybayes)
library(parameters)
library(modelr)

#### Plotting Variables ####

if(file.exists("Vis")){
  print("Vizualization directory previously created.")
}else{
  dir.create("Vis")
  print("Vizualization directory created.")
}

if(any(!exists('pd') & !exists('makinStuffPretty') & !exists('plot_colors3'))){
  pd <- position_dodge(0.3)
  
  makinStuffPretty <-   theme_classic() +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", size=1),
      #legend.key = element_blank(),
      strip.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      text = element_text(size = 22),
      axis.text = element_text(colour = "black", size = 22),
      legend.title = element_text(size = 26),
      axis.ticks = element_line(colour = "black", size = 1),
      legend.text = element_text(size = 22)
    )
  
  plot_colors3 <- c("navyblue", "orangered1")
} else {"Plotting variables previously created in 3_Descriptives.R."}

#### STATS ####

#define contrasts and prep data for model

halfsum_contrasts = function (...) contr.sum(...) * 0.5

modeldat <- goodtrialmax %>%
  filter(modality != "Training") %>%
  mutate(
    scaled_maxFT3 = scale(maxFT3)[,1],
    scaled_maxFT3roc = scale(maxFT3roc)[,1],
    trial = as.double(even_trial),
    trial = case_when(trial > 5 ~ trial - 20,
                      TRUE ~ trial),
    cen_trial = trial - 3,
    trialroot = trial ^ (1/2),
    cen_trialroot = trialroot - min(trialroot) - 
      ((max(trialroot) - min(trialroot)) / 2),
    scaled_trial = scale(cen_trial)[,1],
    scaled_trialroot = scale(cen_trialroot)[,1],
    group = fct_relevel(as.factor(group), c("MI-2", "MI-10", "OE")),
    modality = fct_relevel(as.factor(modality), c("Pre", "Post")),
    brick = as.factor(brick))

contrasts(modeldat$brick) <- halfsum_contrasts

contrasts(modeldat$group) <- matrix(c(0, 0.5, -0.5, 2/3, -1/3, -1/3),ncol = 2)

premodeldat <- filter(modeldat, modality == "Pre" & group != 'MI-2') %>%
  mutate(group = fct_drop(group))

contrasts(premodeldat$group) <- halfsum_contrasts

premodeldatMI2 <- filter(modeldat,modality == "Pre" & 
                           trial == 1 & group == 'MI-2') %>%
  mutate(group = fct_drop(group))


#define available CPU cores
n_cores <- floor(parallel::detectCores())

#pre-training
model <- scaled_maxFT3 ~ group * scaled_trial * brick + 
  (1 + scaled_trial * brick | db_id)
model1 <- scaled_maxFT3 ~ group * scaled_trialroot * brick + 
  (1 + scaled_trialroot * brick | db_id)
model2 <- scaled_maxFT3 ~ group * scaled_trialroot * brick * order + 
  (1 + scaled_trialroot * brick | db_id)

altmod <- brm(model,
              data = premodeldat,
              prior = c(
                prior(normal(0, 2), class = b),
                prior(exponential(1), class = sd),
                prior(exponential(1), class = sigma)
              ),
              iter = 10000,
              chains = n_cores, cores = n_cores, seed = 64,
              control = list(adapt_delta = .99, max_treedepth = 20),
              save_pars = save_pars(all = TRUE)
)
altmod <- add_criterion(altmod, "loo", moment_match = TRUE, reloo = TRUE)

altmod1 <- brm(model1,
               data = premodeldat,
               prior = c(
                 prior(normal(0, 2), class = b),
                 prior(exponential(1), class = sd),
                 prior(exponential(1), class = sigma)
               ), 
               iter = 10000,
               chains = n_cores, cores = n_cores, seed = 64,
               control = list(adapt_delta = .99, max_treedepth = 20),
               save_pars = save_pars(all = TRUE)
)
altmod1 <- add_criterion(altmod1, "loo", moment_match = TRUE, reloo = TRUE)

altmod2 <- brm(model2,
               data = premodeldat,
               prior = c(
                 prior(normal(0, 2), class = b),
                 prior(exponential(1), class = sd),
                 prior(exponential(1), class = sigma)
               ),
               iter = 10000,
               chains = n_cores, cores = n_cores, seed = 64,
               control = list(adapt_delta = .99, max_treedepth = 20),
               save_pars = save_pars(all = TRUE)
)
altmod2 <- add_criterion(altmod2, "loo", moment_match = TRUE, reloo = TRUE)

loo_compare(altmod, altmod1, altmod2) #altmod 1 is best

res <- tibble(model_parameters(altmod1, ci = 0.95, ci_method = "hdi",
                               test = c("pd", "rope"),rope_ci = 0.95))

altmod1_fitted <- ungroup(premodeldat) %>%
  data_grid(
    scaled_trialroot = unique(scaled_trialroot),
    group = unique(group),
    brick = unique(brick)
  )%>%
  add_linpred_draws(altmod1, re_formula = NA, ndraws = 2000)%>%
  mutate(
    maxFT3 = .linpred * sd(modeldat$maxFT3) + mean(modeldat$maxFT3),
    cen_trialroot = scaled_trialroot * sd(modeldat$cen_trialroot) + 
      mean(modeldat$cen_trialroot),
    trialroot = cen_trialroot + min(modeldat$trialroot) + 
      ((max(modeldat$trialroot) - min(modeldat$trialroot)) / 2),
    trial = trialroot^2
  )

ME_MI10pre <- ggplot(altmod1_fitted, 
                     aes(y = maxFT3, x = trial, linetype = brick, 
                         fill = brick, group = brick)) +
  stat_lineribbon(.width = .90, point_interval = mean_hdi, 
                  colour = NA,alpha = 0.5) +
  stat_summary(fun = function(x) mean(x), geom = "line", size = 1.5)  +
  scale_linetype_manual(name = "Brick", 
                        values = c(Large = 'solid',Small = 'dotted')) +
  scale_fill_manual(name = "Brick", values = alpha(c("grey20", "grey70"), 0.5)) +
  facet_wrap(~ group, ncol = 2) +
  xlab("Trial") +
  ylab("Maximum Load Force (N)") +
  labs(linetype = "Brick") +
  guides(fill = "none")+
  makinStuffPretty


show(ME_MI10pre)

ggsave(filename = "ME_MI10pre.pdf", #.eps broken with transparency
       path = "./Vis/",
       plot = ME_MI10pre,
       dpi = 1200,
       width = 20,
       height = 20,
       units = "cm")

#in MI-2 group

modelMI2 <- scaled_maxFT3 ~ brick + (1 | db_id)

altmodMI2 <- brm(modelMI2,
                 data = premodeldatMI2,
                 prior = c(
                   prior(normal(0, 2), class = b),
                   prior(exponential(1), class = sd),
                   prior(exponential(1), class = sigma)
                 ),
                 iter = 10000,
                 chains = n_cores, cores = n_cores, seed = 64,
                 control = list(adapt_delta = .99, max_treedepth = 20),
                 save_pars = save_pars(all = TRUE)
)

resMI2 <- tibble(model_parameters(altmodMI2, ci = 0.95, ci_method = "hdi",
                                  test = c("pd", "rope"),rope_ci = 0.95))

altmodMI2_emm <- emmeans(altmodMI2, ~  brick)
altmodMI2_draws <- gather_emmeans_draws(altmodMI2_emm)

group_mean_draws <- altmodMI2_draws %>%
  mutate(
    maxFT3 = .value * sd(modeldat$maxFT3) + mean(modeldat$maxFT3))

MI2pre <- ggplot(group_mean_draws, aes(x = brick, y = maxFT3)) +
  stat_pointinterval(.width = c(0.6, 0.9), point_interval = median_hdi)+
  xlab("Brick Size") +
  ylab("Maximum Load Force (N)") +
  makinStuffPretty

show(MI2pre)

ggsave(filename = "MI2pre.eps",
       path = "./Vis/",
       plot = MI2pre,
       dpi = 1200,
       width = 15,
       height = 20,
       units = "cm")


#what happens post training

altmod1post <- brm(model1,
                   data = filter(modeldat, modality == "Post"),
                   prior = c(
                     prior(normal(0, 2), class = b),
                     prior(exponential(1), class = sd),
                     prior(exponential(1), class = sigma)
                   ),
                   iter = 10000,
                   chains = n_cores, cores = n_cores, seed = 64,
                   control = list(adapt_delta = .99, max_treedepth = 20),
                   save_pars = save_pars(all = TRUE)
)

respost <- tibble(parameters(altmod1post, ci = 0.95, ci_method = "hdi",
                             test = c("pd", "rope"),rope_ci = 0.95))
# group brick interaction
altmod1post_emm <- emmeans(altmod1post, ~  group * brick)
altmod1post_draws <- gather_emmeans_draws(altmod1post_emm)

group_mean_draws <- altmod1post_draws %>%
  mutate(
    maxFT3 = .value * sd(modeldat$maxFT3) + mean(modeldat$maxFT3),
    group = factor(group, levels = c("MI-10", "MI-2", "OE")),
    brick = ifelse(brick == "Large", "Large", "Small"))

all_respost <- ggplot(group_mean_draws, aes(x = brick, y = maxFT3)) +
  stat_pointinterval(.width = c(0.6, 0.9), point_interval = median_hdi)+
  facet_wrap(~ group, ncol = 3) +
  xlab("Brick Size") +
  ylab("Mean Maximum Force (N)") +
  makinStuffPretty


show(all_respost)

ggsave(filename = "all_respost.eps",
       path = "./Vis/",
       plot = all_respost,
       dpi = 1200,
       width = 20,
       height = 20,
       units = "cm")

#repeat for FT3maxROC

#pre data
model1roc <- scaled_maxFT3roc ~ group * scaled_trialroot * brick + 
  (1 + scaled_trialroot | db_id / brick)

altmod1roc <- brm(model1roc,
                  data = premodeldat,
                  prior = c(
                    prior(normal(0, 2), class = b),
                    prior(exponential(1), class = sd),
                    prior(exponential(1), class = sigma)
                  ),
                  iter = 10000,
                  chains = n_cores, cores = n_cores, seed = 64,
                  control = list(adapt_delta = .99, max_treedepth = 20),
                  save_pars = save_pars(all = TRUE)
)

resroc <- tibble(parameters(altmod1roc, ci = 0.95, ci_method = "hdi", 
                            test = c("pd", "rope"), rope_ci = 0.95))

altmod1roc_fitted <- ungroup(premodeldat) %>%
  data_grid(
    scaled_trialroot = unique(scaled_trialroot),
    group = unique(group),
    brick = unique(brick)
  )%>%
  add_linpred_draws(altmod1roc, re_formula = NA, ndraws = 2000)%>%
  mutate(
    maxFT3roc = .linpred * sd(modeldat$maxFT3roc) + mean(modeldat$maxFT3roc),
    cen_trialroot = scaled_trialroot * sd(modeldat$cen_trialroot) + 
      mean(modeldat$cen_trialroot),
    trialroot = cen_trialroot + min(modeldat$trialroot) + 
      ((max(modeldat$trialroot) - min(modeldat$trialroot)) / 2),
    trial = trialroot^2
  )

ME_MI10rocpre <- ggplot(altmod1roc_fitted, 
                        aes(y = maxFT3roc, x = trial, linetype = brick, 
                            fill = brick, group = brick)) +
  stat_lineribbon(.width = .90, point_interval = mean_hdi, 
                  colour = NA, alpha = 0.5) +
  stat_summary(fun = function(x) mean(x), geom = "line", size = 1.5) +
  scale_linetype_manual(name = "Brick", 
                        values = c(Large = 'solid', Small = 'dotted')) +
  scale_fill_manual(name = "Brick", values = alpha(c("grey20", "grey70"), 0.5)) +
  facet_wrap(~ group, ncol = 2) +
  xlab("Trial") +
  ylab("Maximum Load Force Rate (N/s)") +
  labs(linetype = "Brick") +
  guides(fill = "none")+
  makinStuffPretty

show(ME_MI10rocpre)

ggsave(filename = "ME_MI10rocpre.pdf", #.eps broken with transparency
       path = "./Vis/",
       plot = ME_MI10rocpre,
       dpi = 1200,
       width = 20,
       height = 20,
       units = "cm")

#MI-2 pre

modelMI2roc <- scaled_maxFT3roc ~ brick + (1 | db_id)

altmodMI2roc <- brm(modelMI2roc,
                    data = premodeldatMI2,
                    prior = c(
                      prior(normal(0, 2), class = b),
                      prior(exponential(1), class = sd),
                      prior(exponential(1), class = sigma)
                    ),
                    iter = 10000,
                    chains = n_cores, cores = n_cores, seed = 64,
                    control = list(adapt_delta = .99, max_treedepth = 20),
                    save_pars = save_pars(all = TRUE)
)

resMI2roc <- tibble(parameters(altmodMI2roc, ci = 0.95, ci_method = "hdi", 
                               test = c("pd", "rope"), rope_ci = 0.95))

altmodMI2roc_emm <- emmeans(altmodMI2roc, ~ brick)
altmodMI2roc_draws <- gather_emmeans_draws(altmodMI2roc_emm)

roc_group_mean_draws <- altmodMI2roc_draws %>%
  mutate(
    maxFT3 = .value * sd(modeldat$maxFT3roc) + mean(modeldat$maxFT3roc))

MI2preroc <- ggplot(roc_group_mean_draws, aes(x = brick, y = maxFT3)) +
  stat_pointinterval(.width = c(0.6, 0.9), point_interval = median_hdi)+
  xlab("Brick Size") +
  ylab("Maximum Load Force Rate (N/s)") +
  makinStuffPretty

show(MI2preroc)

ggsave(filename = "MI2preROC.eps",
       path = "./Vis/",
       plot = MI2preroc,
       dpi = 1200,
       width = 15,
       height = 20,
       units = "cm")

#what happens post training

altmod1rocpost <- brm(model1roc,
                      data = filter(modeldat,modality=="Post"),
                      prior = c(
                        prior(normal(0, 2), class = b),
                        prior(exponential(1), class = sd),
                        prior(exponential(1), class = sigma)
                      ),
                      iter = 10000,
                      chains = n_cores, cores = n_cores, seed = 64,
                      control = list(adapt_delta = .99, max_treedepth = 20),
                      save_pars = save_pars(all = TRUE)
)

resrocpost <- tibble(parameters(altmod1rocpost, ci = 0.95, ci_method = "hdi", 
                                test = c("pd", "rope"), rope_ci = 0.95))

#no credible effects

#rating
modelRatedat <- goodtrialmax %>%
  mutate(
    scale_rating = scale(rating)[,1],
    trial = as.double(even_trial),
    trialroot = trial^(1 / 2),
    cen_trialroot = trialroot - min(trialroot) - 
      ((max(trialroot) - min(trialroot)) / 2),
    scaled_trialroot = scale(cen_trialroot)[,1],
    group = fct_relevel(as.factor(group), c("OE", "MI-10", "MI-2")),
    brick = as.factor(brick),
    type = as.factor(case_when(trial < 11 & group == "MI-10" | 
                                 trial > 40 & group != "MI-10" ~ "OE",
                               trial < 3 & group == "MI-2" | 
                                 trial > 40 & group != "MI-2" ~ "OE",
                               group == "OE" ~ "OE",
                               TRUE ~ "MI"))) %>%
  select(db_id, group, type, brick, trial, trialroot, 
         cen_trialroot, scaled_trialroot, rating, scale_rating)

contrasts(modelRatedat$brick) <- halfsum_contrasts
contrasts(modelRatedat$type) <- halfsum_contrasts


modelRateMIdat <- filter(modelRatedat, group != "OE") %>%
  mutate(group = fct_drop(group))
contrasts(modelRateMIdat$group) <- halfsum_contrasts

modelRateMEdat <- filter(modelRatedat,group == "OE")

#for MI

modelRate <- scale_rating ~ group * brick * type  + (1 | db_id / brick / type)


altmodRateMI <- brm(modelRate,
                    data = modelRateMIdat,
                    prior = c(
                      prior(normal(0, 2), class = b),
                      prior(exponential(1), class = sd),
                      prior(exponential(1), class = sigma)
                    ),
                    iter = 10000,
                    chains = n_cores, cores = n_cores, seed = 64,
                    control = list(adapt_delta = .99, max_treedepth = 20),
                    save_pars = save_pars(all = TRUE)
)

resrate <- tibble(model_parameters(altmodRateMI, ci = 0.95, ci_method = "hdi", 
                                   test = c("pd", "rope"), rope_ci = 0.95))

altmodRate_emm <- emmeans(altmodRateMI, ~  brick * group * type)
altmodRate_draws <- gather_emmeans_draws(altmodRate_emm)

group_mean_rate_draws <- altmodRate_draws %>%
  mutate(
    rate = .value * sd(modelRatedat$rating) + mean(modelRatedat$rating)
  )

Brick_MIrate <- ggplot(group_mean_rate_draws, aes(x = brick, y = rate)) +
  stat_pointinterval(.width = c(0.6, 0.9), point_interval = median_hdi, 
                     position = pd, fatten_point = 2, aes(shape = group)) +
  xlab("Brick Size") +
  scale_y_continuous(name = "Heaviness Rating (1-10)", 
                     breaks = c(2, 4, 6, 8, 10), limits = c(0, 10)) +
  scale_colour_discrete(type = c("grey20", "grey70"), name = "Brick", 
                        labels = c("Large", "Small")) +
  scale_shape_discrete(name = "Group") +
  facet_wrap(~ type, ncol = 2) +
  makinStuffPretty

show(Brick_MIrate)

ggsave(filename = "Brick_MIrate.eps",
       path = "./Vis/",
       plot = Brick_MIrate,
       dpi = 1200,
       width = 20,
       height = 20,
       units = "cm")

#for ME

modelRateME <- scale_rating ~ brick * scaled_trialroot + 
  (1 + scaled_trialroot | db_id / brick)

altmodRateME <- brm(modelRateME,
                    data = modelRateMEdat,
                    prior = c(
                      prior(normal(0, 2), class = b),
                      prior(exponential(1), class = sd),
                      prior(exponential(1), class = sigma)
                    ),
                    iter = 10000,
                    chains = n_cores, cores = n_cores, seed = 64,
                    control = list(adapt_delta = .99, max_treedepth = 20),
                    save_pars = save_pars(all = TRUE)
)

resMErate<- tibble(model_parameters(altmodRateME, ci = 0.95, ci_method = "hdi", 
                                    test=c("pd", "rope"), rope_ci = 0.95))

#main effect of brick
altmodMERate_emm <- emmeans(altmodRateME, ~  brick)
altmodMERate_draws <- gather_emmeans_draws(altmodMERate_emm)

group_mean_rate_draws <- altmodMERate_draws %>%
  mutate(
    rate = .value * sd(modelRatedat$rating) + mean(modelRatedat$rating))

Brick_MErate <- ggplot(group_mean_rate_draws, aes(x = brick, y = rate)) +
  stat_pointinterval(.width = c(0.6, 0.9), point_interval = median_hdi,
                     position = pd) +
  scale_y_continuous(name = "Heaviness Rating (1-10)", 
                     breaks = c(2, 4, 6, 8, 10), limits = c(0, 10)) +
  xlab("Brick Size") +
  scale_color_discrete(type = c("grey20", "grey70"), name = "Brick", 
                       labels = c("Large", "Small")) +
  makinStuffPretty

show(Brick_MErate)

ggsave(filename = "Brick_MErate.eps",
       path = "./Vis/",
       plot = Brick_MErate,
       dpi = 1200,
       width = 20,
       height = 20,
       units = "cm")
