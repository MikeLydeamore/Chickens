system("./runSIR")
df <- lapply(1:500, function(x) {
  df <- read.csv(paste("~/Documents/PhD/MarkovChainSimulator/outputs/SIR/run-", x, ".csv", sep=""))
  df$sim_number <- x
  df <- melt(df, id.vars=c("t", "sim_number"))
  return (df)
})

df <- rbindlist(df)
df$variable <- factor(df$variable, levels=c("S","I","R"))

df_deterministic <- read.csv("~/Documents/PhD/MarkovChainSimulator/SIRDeterministic.csv")
df_deterministic <- df_deterministic[df_deterministic$t < max(df$t), ]
df_deterministic <- melt(df_deterministic, id.vars="t")
df_deterministic$variable <- factor(df_deterministic$variable, levels=c("S","I","R"))
p <- ggplot(df) + geom_line(aes(x=t, y=value, group=sim_number), alpha = 0.2, colour="gray") + facet_grid(.~variable) +
  geom_line(data = df_deterministic, aes(x=t, y=value)) + labs("y"="Number of Individuals", "x"="Time")


save_plot("/Users/mlydeamore/cloudstor/JCompSimPaper/SIROutput.png", p, base_aspect_ratio = 2.6)

df <- lapply(1:500, function(x) {
  df <- read.csv(paste("~/Documents/PhD/MarkovChainSimulator/outputs/Scabies/run-", x, ".csv", sep=""))
  df2 <- data.frame("t"=df$t, "S"=df$S+df$S2, "I"=df$IA+df$IAHat+df$IHat+df$I2C+df$I2CHat+df$I2+df$I2Hat, "E"=df$GHat+df$YHat+df$I2CHat+df$I2Hat+df$IHat)
  df2$sim_number <- x
  df2 <- melt(df2, id.vars=c("t", "sim_number"))
  return (df2)
})

df <- rbindlist(df)
df$variable <- factor(df$variable, levels=c("S","I","E"))

df_deterministic <- read.csv("~/Documents/PhD/MarkovChainSimulator/ScabiesDeterministic.csv")
df_deterministic <- with(df_deterministic, data.frame("t"=t, "S"=S+S2, "I"=IA+IAHat+IHat+I2C+I2CHat+I2+I2Hat, "E"=GHat+YHat+I2CHat+I2Hat+IHat))
df_deterministic <- df_deterministic[df_deterministic$t < max(df$t), ]
df_deterministic <- melt(df_deterministic, id.vars="t")
df_deterministic$variable <- factor(df_deterministic$variable, levels=c("S","I","E"))
p <- ggplot(df) + geom_line(aes(x=t, y=value, group=sim_number), alpha = 0.2, colour="gray") + facet_grid(.~variable) +
  geom_line(data = df_deterministic, aes(x=t, y=value)) + labs("y"="Number of Individuals", "x"="Time")


save_plot("/Users/mlydeamore/cloudstor/JCompSimPaper/ScabiesOutput.png", p, base_aspect_ratio = 2.6)
