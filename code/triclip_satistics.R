rm(list=ls(all=TRUE)) 
library(afex)
library(emmeans)
library(lmerTest)
setwd('C:\\Users\\trace\\OneDrive\\Documents\\Rausch Lab\\TriClip Experiment')


TriClipData = read.csv('TriClipXT_Statistics.csv')
TriClipData_Flow = read.csv('TriClipXT_Flow_Statistics.csv')


TriClipData$Intervention = factor(TriClipData$Intervention)
TriClipData$Test = factor(TriClipData$Test)
TriClipData$Pin = factor(TriClipData$Pin)
TriClipData_noDiseased = subset(TriClipData, Intervention != "Diseased")
TriClipData_Flow_noDiseased = subset(TriClipData_Flow, Intervention != "Diseased")




cat("\n----------------------------------------------\n")
cat("\n           Force Mixed Effects                \n")
cat("\n----------------------------------------------\n")

ForceDiff.mixed = mixed(ForceDifference ~ Intervention * Pin + (1|Test), data = TriClipData_noDiseased)
print(ForceDiff.mixed[["anova_table"]])


forcediff_emm = emmeans(ForceDiff.mixed$full_model, ~ Intervention * Pin)
print(summary(forcediff_emm, infer = c(TRUE, TRUE), null = 0))

print(pairs(forcediff_emm, by = "Intervention"))


cat("\n----------------------------------------------\n")
cat("\n           Flow Stats                         \n")
cat("\n----------------------------------------------\n")

FlowDiff.mixed = mixed(FlowDifference ~ Intervention + (1|Test), data = TriClipData_Flow_noDiseased)
print(FlowDiff.mixed[["anova_table"]])


flowdiff_emm = emmeans(FlowDiff.mixed$full_model, ~ Intervention)
print(summary(flowdiff_emm, infer = c(TRUE, TRUE), null = 0))
