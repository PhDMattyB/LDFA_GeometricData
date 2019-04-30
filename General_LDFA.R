##############################
## Linear discriminant function analysis body shape data
##
## Matt Brachmann (PhDMattyB)
##
## 2019-03-13
##
##############################

setwd('~/PhD/Morphometrics/Working Directory')

library(tidyverse)
library(wesanderson)
library(patchwork)
library(janitor)
library(devtools)
library(skimr)

theme_set(theme_bw())

# Other packages to load

## Read in your morphometric data
## I used the filter function to only look at polymorphic populations
Morpho_data = read_csv('AllLakes_AllometryExcluded_PWS_Combined.csv') %>% 
  filter(LaMorph %in% c('G.SB', 'G.PI', 'S.LGB', 'S.PI', 
                        'S.PL', 'T.LGB', 'T.PL', 'T.SB', 
                        'V.BR', 'V.SIL'))
## Needed to basically copy and paste the data for a benthic and pelagic morph
## This allowed my to make Vectors corresponding to different benthic and pelagic 
## morph pairs
Morpho_SLGB = Morpho_data %>% 
  slice(78:103) %>% 
  mutate(LaMorph2 = as.factor(case_when(
    LaMorph == 'S.LGB' ~ 'S.LGB2'))) %>%
  select(-LaMorph) %>% 
  rename(LaMorph = LaMorph2) %>% 
  select(id:BP2, LaMorph, Sex:BPLD1,
         contains('PW'), UNIX:CS)


Morpho_TPL = Morpho_data %>% 
  slice(297:354) %>% 
  mutate(LaMorph2 = as.factor(case_when(
    LaMorph == 'T.PL' ~ 'T.PL2'))) %>%
  select(-LaMorph) %>% 
  rename(LaMorph = LaMorph2) %>% 
  select(id:BP2, LaMorph, Sex:BPLD1,
         contains('PW'), UNIX:CS)

Morpho_cleaned = bind_rows(Morpho_data, Morpho_SLGB, Morpho_TPL) %>% 
  group_by(LaMorph)

## The vectors I mentioned above
Morpho_cleaned = mutate(.data = Morpho_cleaned,
              Vector = as.factor(case_when(
                LaMorph == "G.SB" ~ "GSBPI",
                LaMorph == 'G.PI' ~ 'GSBPI',
                LaMorph == 'S.LGB' ~ 'SLGBPI',
                LaMorph == 'S.PI'~ 'SLGBPI',
                LaMorph == 'S.LGB2' ~ 'SLGBPL',
                LaMorph == 'S.PL' ~ 'SLGBPL',
                LaMorph == 'T.LGB' ~ 'TLGBPL',
                LaMorph == 'T.PL' ~ 'TLGBPL',
                LaMorph == 'T.SB' ~ 'TSBPL', 
                LaMorph == 'T.PL2' ~ 'TSBPL',
                LaMorph == 'V.BR' ~ 'VSILBR', 
                LaMorph == 'V.SIL' ~ 'VSILBR')))  

## Filter for the vector you want to analyze body shape for
Morpho_grouped = Morpho_cleaned %>% 
  arrange(BP) %>% 
  group_by(Vector) %>% 
  filter(Vector == 'VSILBR')

## The LDFA analysis, put in all of your partial warp and uniform component measurements
ldfa = lda(Morpho_grouped$BP ~ Morpho_grouped$PW1X+ Morpho_grouped$PW1Y+ Morpho_grouped$PW2X+ Morpho_grouped$PW2Y+
      Morpho_grouped$PW3X+ Morpho_grouped$PW3Y+ Morpho_grouped$PW4X+ Morpho_grouped$PW4Y+
      Morpho_grouped$PW5X+ Morpho_grouped$PW5Y+ Morpho_grouped$PW6X+ Morpho_grouped$PW6Y+
      Morpho_grouped$PW7X+ Morpho_grouped$PW7Y+ Morpho_grouped$PW8X+ Morpho_grouped$PW8Y+
      Morpho_grouped$PW9X+ Morpho_grouped$PW9Y+ Morpho_grouped$PW10X+ Morpho_grouped$PW10Y+
      Morpho_grouped$PW11X+ Morpho_grouped$PW11Y+ Morpho_grouped$PW12X+ Morpho_grouped$PW12Y+
      Morpho_grouped$PW13X+ Morpho_grouped$PW13Y+ Morpho_grouped$PW14X+ Morpho_grouped$PW14Y+
      Morpho_grouped$PW15X+ Morpho_grouped$PW15Y+ Morpho_grouped$PW16X+ Morpho_grouped$PW16Y+
      Morpho_grouped$PW17X+ Morpho_grouped$PW17Y+ Morpho_grouped$PW18X+ Morpho_grouped$PW18Y+
      Morpho_grouped$PW19X+ Morpho_grouped$PW19Y+ Morpho_grouped$UNIX+ Morpho_grouped$UNIY, CV = T)

lda_predict = predict(ldfa)
apply(lda_predict$posterior, MARGIN = 1, FUN = max)


table = table(Morpho_grouped$BP, lda_predict$class)
sum(table[row(table) == col(table)])/sum(table)

#leave one out cross validation
##CV needs to be true in LDA
table2 = table(Morpho_grouped$BP, ldfa$class)

##Calculate the re-substitution error for the cross validation
sum(table2[row(table2) == col(table2)])/sum(table2)

