#### Import data for the SWI project ####

#file I/O

path <- "~/Documents/Dalhousie/PhD-Neuroscience/Comps/Comp #1 - David Westwood/SWI_Analysis"
setwd(path)

if(file.exists("Vis")){
  print("Vizualization directory previously created.")
}else{
  dir.create("Vis")
  print("Vizualization directory created.")
}

#### Requirements ####

library(tidyverse)


#### IMPORTS ####

#import demographics
demographics <- read_csv("./_Data/Masterlist.csv",
                         col_names = c("Participant ID", "Sex", "Age", "Group",
                                       "Order", "Cable Orientation", "KVIQ-V",
                                       "KVIQ-K"),
                         col_types = c("n", "f", "n", "f", "f", "n", "n", "n"),
                         skip = 1
) %>%
  mutate(Group = fct_recode(Group,`MI-10` = "MI")) %>%
  select(-`Cable Orientation`)

#import trial data
SWI_datasheet <- read_csv("./_Data/SWI_datasheet.csv") %>%
  mutate(Trial = (as.numeric(Block) - 1) * 10 + Trial) %>%
  filter(!is.na(`Group (from masterlist)`))

trialinfo <- SWI_datasheet %>%
  select(db_id = Participant_ID, trialcount = Trial,
                    `Heaviness Rating`, comment, `FT3-Switch status`)

partinfo <- SWI_datasheet %>%
  group_by(Participant_ID) %>%
  summarise(group = unique(`Group (from masterlist)`),
            order = unique(`Order (from masterlist)`)) %>%
  rename(db_id = Participant_ID)%>%
  mutate(group = ifelse(group == "PP", "OE", group))

#import transducer & EMG data
raw.files <- dir('./_Data/labview', recursive = TRUE,
                 full.names = TRUE, pattern="[[:digit:]]{2}.txt$")

taskdat <- map_df(raw.files, function(f) {
  df <- read_delim(f, col_names = F, delim = ',', col_types = 'dddddddddd') %>%
    mutate(id=as.double(str_sub(basename(f), 2, 4)))
}) %>%
  select(id, X2:X10)

colnames(taskdat) <- c('db_id', 'goggles', 'switch', 'EMG', 'FT1', 'FT2', 'FT3', 
                       'FT4', 'FT5', 'FT6')

#import Experiment Builder Data
exp_B.files <- dir('./_Data/Data R Comp', recursive = TRUE,
                   full.names = TRUE, pattern="RESULTS_FILE.txt$")

exp_Bdat <- map_df(exp_B.files, function(f) {
  df <- read_delim(f, col_names = T, delim = '\t', col_types = 'ddddc') %>%
    mutate(reaction_time=as.double(reaction_time),
           id=as.double(str_sub(dirname(f), -3)), .before=1
    )
}) %>%
  arrange(id,go_cue_time)
