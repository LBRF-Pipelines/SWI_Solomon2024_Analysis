#### Visualize data analysis steps for data cleaning in the SWI project ####

#### TROUBLESHOOTING PLOTS ####

# Plot raw data to see trialcount with trials are denoted by the (goggles) and
# the buffer by (trial_samples).
# If you are clearing variables to save memory, run the function on the object
# created by the first pipe in 1_Analysis.

raw <- ggplot(data = filter(trialdat, db_id == 1), 
              mapping = aes(x = 1:length(goggles))) + 
  geom_line(mapping=aes(y = trial_samples), col = 'red') + 
  geom_line(mapping=aes(y = goggles), col = 'cyan') + 
  geom_line(mapping=aes(y = trial_index), col = 'blue') +
  geom_line(mapping=aes(y = trialcount), col = 'orange')
#show(raw)

#troubleshooting plot trial
trouble <- ggplot(data = filter(justtrialdat, db_id == 53& trial == 1),
                  aes(x = 1:length(EMG)))+
  geom_line(aes(y = EMG / 10), color = "grey") +
  geom_line(aes(y = (switch / 10) - 0.25), color = "red") +
  geom_line(aes(y = goggles), color = "blue") +
  geom_line(aes(y = FT3), color = "darkgreen")
#show(trouble)