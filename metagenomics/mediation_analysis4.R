# for mediation analysis
# 参考：https://towardsdatascience.com/doing-and-reporting-your-first-mediation-analysis-in-r-2fe423b92171
rm(list = ls())
library(mediation)

# pareparing the datasets for this part of analysis
df <- read.csv("mediation_16S_metabo_serum1.csv", header = T, row.names = 1,
               check.names = F)
# scale the dataframe
df2 <- as.data.frame(scale(df))

# Step #1: The total effect
# The total effect describes the total effect that the independent variable (iv) 
# has on the dependent variable (dv)
fit.totaleffect=lm(AST~Clostridium_sensu_stricto_1,df2)
summary(fit.totaleffect)

# Step #2: The effect of the IV onto the mediator
fit.mediator=lm(`Aminoadipic acid`~Clostridium_sensu_stricto_1,df2)
summary(fit.mediator)

# Step #3: The effect of the mediator on the dependent variable
fit.dv=lm(AST~Clostridium_sensu_stricto_1+`Aminoadipic acid`,df2)
summary(fit.dv)

# Step #4: Causal Mediation Analysis
results = mediate(fit.mediator, fit.dv, treat='Clostridium_sensu_stricto_1', 
                  mediator='Aminoadipic acid', boot=T)
summary(results)

# 部分结果解释
# (1) ACME stands for average causal mediation effects. This is the indirect 
# effect of the IV (Clostridium_sensu_stricto_1) on the DV (TC) 
# that goes through the mediator (Aminoadipic acid)

# (2) ADE stands for average direct effects. It describes the direct effect of 
# the IV on the DV. Again, this is not new information. We have calculated this 
# effect in step #3: the direct effect of the IV on the DV when controlling for 
# the mediator.

# (3) Total Effect stands for the total effect (direct + indirect) of the IV 
# onto the DV.

# (4) Prop. Mediated describes the proportion of the effect of the IV on the 
# DV that goes through the mediator.

# 反向中介作用分析
# 以代谢物作为因变量，看细菌通过表型（如TC）对Aminoadipic acid的作用

# Step #1: The total effect
# The total effect describes the total effect that the independent variable (iv) 
# has on the dependent variable (dv)
fit.totaleffect2=lm(`Aminoadipic acid`~Clostridium_sensu_stricto_1,df2)
summary(fit.totaleffect2)

# Step #2: The effect of the IV onto the mediator
fit.mediator2=lm(AST~Clostridium_sensu_stricto_1,df2)
summary(fit.mediator2)

# Step #3: The effect of the mediator on the dependent variable
fit.dv2=lm(`Aminoadipic acid`~Clostridium_sensu_stricto_1+AST,df2)
summary(fit.dv2)

# Step #4: Causal Mediation Analysis
results2 = mediate(fit.mediator2, fit.dv2, treat='Clostridium_sensu_stricto_1', 
                   mediator='AST', boot=T)
summary(results2)



