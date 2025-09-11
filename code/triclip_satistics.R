rm(list=ls(all=TRUE)) 
library(afex)
library(emmeans)
library(lmerTest)
setwd('C:\\Users\\trace\\OneDrive\\Documents\\Rausch Lab\\TriClip Experiment')


TriClipData = read.csv('TriClipXT_Statistics.csv')
TriClipData_1pin = read.csv('TriClipXT_Flow_Statistics.csv')


TriClipData$Intervention = factor(TriClipData$Intervention)
TriClipData$Test = factor(TriClipData$Test)
TriClipData$Pin = factor(TriClipData$Pin)
TriClipData_noDiseased = subset(TriClipData, Intervention != "Diseased")
TriClipData_1pin_noDiseased = subset(TriClipData_1pin, Intervention != "Diseased")


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

Flow.mixed = mixed(FlowRate ~ Intervention + (1|Test), data = TriClipData_1pin)
print(Flow.mixed[["anova_table"]])


flow_emm = emmeans(Flow.mixed$full_model, ~ Intervention)
flow_contrasts = contrast(flow_emm, "trt.vs.ctrl", ref = "Diseased")
print(summary(flow_contrasts, infer = c(TRUE, TRUE)))


cat("\n----------------------------------------------\n")
cat("\n           Pressure Stats                     \n")
cat("\n----------------------------------------------\n")

Pressure.mixed = mixed(Pressure ~ Intervention + (1|Test), data = TriClipData_1pin)
print(Pressure.mixed[["anova_table"]])


pressure_emm = emmeans(Pressure.mixed$full_model, ~ Intervention)
pressure_contrasts = contrast(pressure_emm, "trt.vs.ctrl", ref = "Diseased")
print(summary(pressure_contrasts, infer = c(TRUE, TRUE)))
