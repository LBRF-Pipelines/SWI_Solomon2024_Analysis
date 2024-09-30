#### Summarize the descriptive results and demographics for the SWI project ####

#### Requirements ####

source("./_Scripts/1_Analysis.R")

ifelse("tidyverse" %in% (.packages()),
       "Tidyverse loaded from 1_Analysis.R.", library(tidyverse))

#### Plotting Variables ####

if(file.exists("Vis")){
  print("Vizualization directory previously created.")
}else{
  dir.create("Vis")
  print("Vizualization directory created.")
}

vlines_data <- tibble(
  group=c("MI-2", "MI-2", "MI-10", "MI-10", "OE", "OE"),
  lines=c(1.5, 20.5, 5.5, 20.5, 5.5, 20.5)
)

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

#### DEMOGRAPHICS ####

finalparticipants <- demographics %>%
  filter(`Participant ID` %in% unique(goodtrialmax$db_id))

finaln <- finalparticipants %>%
  group_by(Group) %>%
  summarize(n = n(),
            female = sum(Sex == "F"),
            male = sum(Sex == "M"),
            meanAge = mean(Age),
            sdAge = sd(Age),
            meanKVIQ_V = mean(`KVIQ-V`),
            sdKVIQ_V = sd(`KVIQ-V`),
            meanKVIQ_K = mean(`KVIQ-K`),
            sdKVIQ_K = sd(`KVIQ-K`)
  )

#dropped:

#P2 only 29 trials
#P5 issue in data collection
#P6 lost trial 1
#P12 only 40 trials
#P15 illusion did not work
#P44 illusion did not work
#P52 illusion did not work

#### RESULT VISUALIZATIONS ####

#max of FT3 on T1

lift_1 <- ggplot(data = filter(goodtrialmax, even_trial==1),
                 aes(x = brick, y = maxFT3)) +
  stat_boxplot(outlier.alpha = 0)+
  geom_point(position = "jitter") +
  facet_wrap(~ group, ncol = 1) +
  xlab("Brick") +
  ylab("Maximum Load Force (N)") +
  makinStuffPretty

show(lift_1)

ggsave(filename = "lift_1.eps",
       path = "./Vis/",
       plot = lift_1,
       dpi = 1200,
       width = 15,
       height = 20,
       units = "cm")

#max of FT3
FT3 <- ggplot(data = goodtrialmax, 
              aes(x = even_trial,y = maxFT3, color = brick, group = brick)) +
  stat_summary(fun = mean,
               fun.min = function(Score) mean(Score) - sd(Score), 
               fun.max = function(Score) mean(Score) + sd(Score), 
               geom = "pointrange", position = pd) +
  geom_vline(data = vlines_data, aes(xintercept = lines))+
  facet_wrap(~ group, ncol = 1) +
  xlab("Exposure") +
  ylab("Maximum Load Force (N)") +
  scale_color_discrete(name = "Brick",type = c("grey20", "grey70"), 
                       labels = c("Large", "Small")) +
  makinStuffPretty

show(FT3)

ggsave(filename = "FT3_max.eps",
       path = "./Vis/",
       plot = FT3,
       dpi = 1200,
       width = 20,
       height = 20,
       units = "cm")

#max of FT3roc

FT3roc <- ggplot(data = goodtrialmax,
                 aes(x = even_trial,y = maxFT3roc, 
                     color = brick, group = brick)) +
  stat_summary(fun = mean,
               fun.min = function(Score) mean(Score) - sd(Score), 
               fun.max = function(Score) mean(Score) + sd(Score), 
               geom = "pointrange", position = pd) +
  geom_vline(data = vlines_data, aes(xintercept = lines))+
  facet_wrap(~ group, ncol = 1) +
  xlab("Exposure") +
  ylab("Maximum Load Force Rate (N/s)") +
  scale_color_discrete(name = "Brick", type=c("grey20", "grey70"), 
                       labels = c("Large", "Small")) +
  makinStuffPretty

show(FT3roc)

ggsave(filename = "FT3roc_max.eps",
       path = "./Vis/",
       plot = FT3roc,
       dpi = 1200,
       width = 20,
       height = 20,
       units = "cm")

#Plot ratings

ratings <- ggplot(data = goodtrialmax,
                  aes(x = even_trial, y = rating, 
                      color = brick, group = brick)) +
  stat_summary(fun = mean,
               fun.min = function(Score) mean(Score) - sd(Score), 
               fun.max = function(Score) mean(Score) + sd(Score), 
               geom = "pointrange",
               position = pd)+
  geom_vline(data = vlines_data, aes(xintercept = lines))+
  facet_wrap(~ group, ncol = 1)+
  ylim(c(0, 10))+
  xlab("Exposure") +
  ylab("Heaviness Rating (1-10)") +
  labs(group = "Brick", color = "Brick")+
  scale_color_discrete(type = c("grey20", "grey70"), 
                       labels = c("Large", "Small")) +
  makinStuffPretty

show(ratings)

ggsave(filename = "ratings.eps",
       path = "./Vis/",
       plot = ratings,
       dpi = 1200,
       width = 40,
       height = 20,
       units = "cm")

# Chronometry Plots

timing <- ggplot(data = filter(MC_timing, movement_type=="MI"),
                 aes(x = trials, y = reaction_time
                     )) + 
  stat_summary(fun = mean,
               fun.min = function(Score) mean(Score) - sd(Score), 
               fun.max = function(Score) mean(Score) + sd(Score), 
               geom = "pointrange",
               position = pd) +
  xlab("Trial") +
  ylab("Movement Duration (ms)") +
  makinStuffPretty

show(timing)

ggsave(filename = "MI_timing.eps",
       path = "./Vis/",
       plot = timing,
       dpi = 1200,
       width = 40,
       height = 20,
       units = "cm")
