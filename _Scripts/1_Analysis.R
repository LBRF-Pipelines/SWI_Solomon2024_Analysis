#### Analyse and clean data for the SWI project ####

#file I/O

path <- "~/Documents/Dalhousie/PhD-Neuroscience/Comps/Comp #1 - David Westwood/SWI_Analysis"
setwd(path)

#### Requirements ####

source("./_Scripts/_Functions/_Analysis.R")
source("./_Scripts/0_Import.R")

ifelse("tidyverse" %in% (.packages()),
       "Tidyverse loaded from 0_Import.R.", library(tidyverse))
library(slider)

#### ANALYSIS ####

fs = 1000
trial_len = 3
trial_samp = trial_len * fs
trial_and_buffer = 5 * fs

trialdat <- taskdat %>%
  group_by(db_id) %>%
  mutate(trial_index = c(0, (goggles - lag(goggles))[-1])) %>%
  # This following line ensures the goggles are open in & close outside of
  # the trial_len window.
  mutate(trial_index = case_when( 
    db_id < 3 & 1:length(trial_index) %in% 
      which(trial_index %in% 1)[
        which(which(trial_index %in% 1) - which(trial_index %in% - 1) >
                (trial_samp - 1201)
        )
      ] ~ 0,
    db_id < 3 & 1:length(trial_index) %in%
      which(trial_index %in% - 1)[
        which(which(trial_index %in% 1) - which(trial_index %in% - 1) > 
                (trial_samp-1201)
        )
      ] ~ 0,
    db_id >= 3 & 1:length(trial_index) %in%
      which(trial_index %in% - 1)[
        which(which(trial_index %in% - 1) - which(trial_index %in% 1) >
                -(trial_samp-1)
        )
      ] ~ 0,
    db_id >= 3 & 1:length(trial_index) %in% 
      which(trial_index %in% 1)[
        which(which(trial_index %in% - 1) - which(trial_index %in% 1) >
                -(trial_samp-1)
        )
      ] ~ 0,
    TRUE~as.double(trial_index))) %>%
  mutate(trial_index = case_when(db_id < 3 & trial_index == -1 ~ 0,
                                 db_id >= 3 & trial_index == 1 ~ 0,
                                 TRUE ~ trial_index)) %>%
  group_by(db_id) %>%
  mutate(trial_index = abs(trial_index),
         trialcount = cumsum(trial_index)) %>%
  mutate(trial_samples = ifelse(slide_dbl(trial_index, 
                                          ~ mean(.x), 
                                          .before = ifelse(unique(db_id) < 3,
                                              trial_and_buffer * (600 / fs) -1,
                                              trial_and_buffer - 1)) > 0,
                                1, 0)
  )

rm(taskdat)

# filter force transducer
cutoff <- c(8)
trialdat <- trialdat %>%
  group_by(db_id) %>%
  mutate(FT1=filtfiltnew(FT1, type = 'low', hz = cutoff, 
                         samples_per_sec = ifelse(unique(db_id) < 3, 600, fs)),
         FT2=filtfiltnew(FT2, type = 'low', hz = cutoff, 
                         samples_per_sec = ifelse(unique(db_id) < 3, 600, fs)),
         FT3=filtfiltnew(FT3, type = 'low', hz = cutoff, 
                         samples_per_sec = ifelse(unique(db_id) < 3, 600, fs)),
         FT4=filtfiltnew(FT4, type = 'low', hz = cutoff, 
                         samples_per_sec = ifelse(unique(db_id) < 3, 600, fs)),
         FT5=filtfiltnew(FT5, type = 'low', hz = cutoff, 
                         samples_per_sec = ifelse(unique(db_id) < 3, 600, fs)),
         FT6=filtfiltnew(FT6, type = 'low', hz = cutoff, 
                         samples_per_sec = ifelse(unique(db_id) < 3, 600, fs))
  )

# isolate data from trials with 2 second buffer and attach metadata
justtrialdat <- trialdat %>%
  filter(trial_samples == 1) %>%
  mutate(modality = case_when(trialcount %in% c(1:10) ~ 'Pre',
                              trialcount %in% c(41:50) ~ 'Post',
                              trialcount %in% c(11:40) ~ 'Training',
                              TRUE ~ as.character(trialcount))) %>%
  inner_join(partinfo) %>%
  mutate(brick = case_when(order == "S-B" & trialcount %% 2 == 0 ~ 'Large',
                           order == "B-S" & trialcount %% 2 == 1 ~ 'Large',
                           order == "B-S" & trialcount %% 2 == 0 ~ 'Small',
                           order == "S-B" & trialcount %% 2 == 1 ~ 'Small',
                           TRUE ~ as.character(trialcount))) %>%
  inner_join(trialinfo) %>%
  select(db_id, order, group, brick, modality, trial = trialcount, goggles, 
         switch, EMG, FT1:FT6, heaviness = `Heaviness Rating`, comment, 
         miss = `FT3-Switch status`)

# re-Zero force transducer to the onset of each trial and calculate rate of 
# change (first derivative)
justtrialdat <- justtrialdat %>%
  group_by(db_id, trial) %>%
  mutate(FT1 = sqrt((FT1-FT1[1])^2),
         FT1roc = n_point_cen_diff((FT1-FT1[1]),1e-3,9),
         FT2 = sqrt((FT2-FT2[1])^2),
         FT2roc = n_point_cen_diff((FT2-FT2[1]),1e-3,9),
         FT3 = sqrt((FT3-FT3[1])^2),
         FT3roc = n_point_cen_diff((FT3-FT3[1]),1e-3,9),
         switchfilt = ifelse(slide_dbl(switch, ~ sd(.x), .before = 5, 
                                       .after = 5) > 0.01, 1, 0),
         samples = as.double(c(0:(n() - 1))),
         time = case_when(db_id < 3 ~ samples / 600,
                          db_id > 2 ~ samples / fs)
  ) %>%
  select(-samples) %>%
  select(1:5, time, everything())

rm(trialdat)

#### DROP INCOMPLETE DATASETS ####

#see who has all 50 trials with data
trialsum <- justtrialdat%>%
  group_by(db_id) %>%
  summarise(max = max(trial))


#drop bad participants and participants with <50 trials
justtrialdat <- justtrialdat %>%
  filter(! db_id %in% c(5, filter(trialsum, max < 50)$db_id))

#### ISOLATE TRIALMAX, REMOVE BAD TRIALS AND CHECK ILLUSION ####

#get channel max for transducer
trialmax <- justtrialdat %>%
  filter(time < trial_len) %>%
  group_by(db_id, trial) %>%
  summarise(group = unique(group),
            order = unique(order),
            trial = unique(trial),
            brick = unique(brick), 
            modality = unique(modality), 
            rating = unique(heaviness),
            comment = unique(comment),
            miss = unique(miss),
            maxFT1 = max(FT1),
            maxFT2 = max(FT2),
            maxFT3 = max(FT3),
            timetomaxFT1 = which(FT1 == max(FT1))[1],
            timetomaxFT2 = which(FT2 == max(FT2))[1],
            timetomaxFT3 = which(FT3 == max(FT3))[1],
            maxFT1roc = max(FT1roc),
            maxFT2roc = max(FT2roc),
            maxFT3roc = max(FT3roc),
            timetomaxFT1roc = which(FT1roc == max(FT1roc))[1],
            timetomaxFT2roc = which(FT2roc == max(FT2roc))[1],
            timetomaxFT3roc = which(FT3roc == max(FT3roc))[1],
            resolvedforce = sqrt((maxFT1^2 + maxFT2^2))) %>%
  mutate(even_trial = cumsum(trial %% 2 != 0),
         db_id = as.factor(db_id),
         timetomaxFT1 = ifelse(db_id == 1, 
                               timetomaxFT1 * (1 / .6), timetomaxFT1),
         timetomaxFT2 = ifelse(db_id == 1, 
                               timetomaxFT2 * (1 / .6), timetomaxFT2),
         timetomaxFT3 = ifelse(db_id == 1, 
                               timetomaxFT3 * (1 / .6), timetomaxFT3),
         timetomaxFT1roc = ifelse(db_id == 1, 
                                  timetomaxFT1roc * (1 / .6), timetomaxFT1roc),
         timetomaxFT2roc = ifelse(db_id == 1, 
                                  timetomaxFT2roc * (1 / .6), timetomaxFT2roc),
         timetomaxFT3roc=ifelse(db_id == 1, 
                                timetomaxFT3roc * (1 / .6), timetomaxFT3roc))

#view trials with comments
badtrials <- trialmax %>%
  group_by(db_id, trial) %>%
  filter(!is.na(comment) | !is.na(miss))

#drop mistrials
trialmax <- trialmax %>%
  filter(!comment %in% c('mistrial')) %>%
  filter(!miss %in% c('Missed'))

#see if illusion worked
maxdiffs <- trialmax %>%
  select(db_id, group, even_trial, brick, rating, maxFT3)%>%
  group_by(db_id, even_trial) %>% 
  pivot_wider(values_from = c(rating, maxFT3), names_from = brick)%>%
  filter(even_trial == 1)%>%
  summarise(diffFT3 = maxFT3_Large - maxFT3_Small,
            diffRate = rating_Small - rating_Large) 

#identify participants where the illusion didn't work
badmaxdiffs <- maxdiffs %>%
  filter(diffRate < 0 & diffFT3 < 0)

#drop participants where illusion did not work
goodtrialmax <- trialmax %>%
  filter(! db_id %in% c(badmaxdiffs$db_id)) %>%
  ungroup() %>%
  mutate(group = fct_recode(group, OE = "OE",`MI-10` = "MI", `MI-2` = "MI2"))

#### Chronometry ####

MC_timing <- exp_Bdat %>%
  inner_join(partinfo, by = c("id"="db_id")) %>%
  mutate(trials = (block_number-1)*10 + trials,
    movement_type = case_when(group == "MI2" & trials > 2 & trials < 41 ~ "MI",
                                   group =="MI" & trials > 10 & trials < 41 ~ "MI",
                                   TRUE ~ "OE"),
    id = as.factor(id)) %>%
  filter(reaction_time != 9999) %>%
  semi_join(goodtrialmax, by = c("id"="db_id", "trials" = "trial")) %>%
  select(c(id, group, trials, movement_type, go_cue_time, 
           movement_start_time, reaction_time))
