set.seed(123)
options(encoding = 'UTF-8')
library(tidyverse)
library(readxl)
library(spatstat)
library(sandwich)
library(lmtest)
MODEL <- '3v1_1p-MN-SV-1000-15000'
path <- paste0('../Results/', MODEL)
grp <- readRDS('../Data/mult_supplementary.RDS')
explanatory <- readRDS('../Data/explanatory.RDS')
data <- readRDS(paste0(path,'/',MODEL,'-regIRFs-chol.RDS')) %>%
  filter(shock == 'exp')
codes <- read_xlsx('../Technical/regioncodes.xlsx')
data <- data %>%
  left_join(codes %>% select(region = Abbrev, name_rus = `Rus Name`, name_eng = Name, OKATO_id))
grid <- read.table('../Technical/rus_grid_reg2.csv', header=TRUE, sep=',', encoding = 'UTF-8') %>%
  mutate(code = ifelse(nchar(code)==10, paste0('0', code), code),
         code = case_when(
           code == '71001000000' ~ '71000000000',
           code == '11001000000' ~ '11000000000',
           TRUE ~ code
         )) %>%
  select(-code_FD, -FD) %>% rename("name" = "name_short")
data <- data %>%
  left_join(grid %>% select(OKATO_id=code, name))  %>%
  relocate(c(OKATO_id, name, name_rus, name_eng), .before = 'variable')
# NORMALIZE IRFs ----
## Own Output Multiplier ---- 
# IRF represents 100 * dlog(y) / dlog(exp) = 100 dy/y * exp/dexp
# Divide by exp to get 100 (dy/y) / dexp - percentage point change in economic activity index 
# Following a unit (rubble per capita) change in exp 
# Then multiply by 0.01 * GRP per capita to normalize shock size to 1% of per capita GRP
multY <- data %>%
  filter(region == source,
         variable == 'y') %>%
  left_join(grp) %>%
  group_by(region) %>%
  mutate(value = Q50,
         sig68 = Q16*Q84>0,
         sig80 = Q10*Q90>0) %>%
  pivot_wider(names_from = 'variable', values_from = 'value') %>%
  mutate(value = y / exp_level * grp * 0.01,
         Q16 = Q16 / exp_level * grp * 0.01,
         Q84 = Q84 / exp_level * grp * 0.01)
## Spill-out Output Multiplier ---- 
# IRF represents 100 * dlog(y) / dlog(exp) = 100 dy/y * exp/dexp
# Divide by exp to get 100 (dy/y) / dexp - percentage point change in economic activity index 
# Following a unit (rubble per capita) change in exp 
# Then multiply by 0.01 * GRP per capita to normalize shock size to 1% of per capita GRP
multspillY <- data %>%
  filter(!region == source,
         variable == 'y') %>%
  left_join(grp %>% 
              left_join(codes %>% select(source = Abbrev, OKATO_id)) %>%
              select(-OKATO_id, grp_source = grp, exp_level_source = exp_level)) %>%
  left_join(grp %>% select(-exp_level)) %>%
  mutate(value = Q50 / exp_level_source * grp_source * 0.01,
         sig68 = Q16*Q84>0) %>%
  group_by(source, h) %>%
  reframe(value = weighted.mean(value, grp), sig = mean(sig68))
## Spill-in Output Multiplier ---- 
# IRF represents 100 * dlog(y) / dlog(exp) = 100 dy/y * exp/dexp
# Divide by exp to get 100 (dy/y) / dexp - percentage point change in economic activity index 
# Following a unit (rubble per capita) change in exp 
# Then multiply by 0.01 * GRP per capita to normalize shock size to 1% of per capita GRP
multspillinY <- data %>%
  filter(!region == source,
         variable == 'y') %>%
  left_join(grp %>% 
              left_join(codes %>% select(source = Abbrev, OKATO_id)) %>%
              select(-OKATO_id, grp_source = grp, exp_level_source = exp_level)) %>%
  left_join(grp %>% select(-exp_level)) %>%
  mutate(value = Q50 / exp_level_source * grp_source * 0.01,
         sig68 = Q16*Q84>0) %>%
  group_by(OKATO_id, name_rus, h) %>%
  reframe(value = weighted.mean(value, grp_source), sig = mean(sig68))
## Own Output Response by Region ----
peakY <- multY %>%
  select(OKATO_id, name_rus, h, value, sig68, grp) %>%
  group_by(OKATO_id, name_rus) %>%
  mutate(absvalue = abs(value)) %>%
  summarize(h = h[which.max(absvalue)], value = value[which.max(absvalue)], sig = sig68[which.max(absvalue)]) %>%
  ungroup() %>%
  mutate(value = paste0(sprintf("%.2f", round(value,2)), ifelse(sig, '*', ''))) %>%
  select(-sig) %>%
  rename_with(~c('OKATO', 'Регион', 'Горизонт', 'Значение'))
horizonY <- multY %>%
  select(OKATO_id, name_rus, h, value, sig68, grp) %>%
  ungroup() %>%
  select(-region) %>%
  filter(h %in% c(0,1,2,3,6,12))
## Output Spillovers and Spill-ins ----
tableregionspillY <- multspillY %>%
  left_join(codes %>% select(source = Abbrev, OKATO_id, OKATO = `Rus Name`)) %>%
  select(-c(source,sig)) %>%
  filter(h %in% c(0,3,6,12)) %>%
  ungroup() %>%
  pivot_wider(names_from = 'h', values_from = 'value') %>%
  left_join(
    multspillinY %>%
      select(-c(name_rus,sig)) %>%
      filter(h %in% c(0,3,6,12)) %>%
      ungroup() %>%
      pivot_wider(names_from = 'h', values_from = 'value', names_prefix = 'SP')
  ) 
# Own Output Response
tableY <- multY %>%
  select(region, OKATO_id, name_rus, h, value, Q16, Q84, sig68, grp) %>%
  group_by(h) %>%
  reframe(o = c(min(value), quantile(value, QUANTILES), max(value), 
                weighted.mean(value, grp), weighted.median(value, grp), mean(sig68),
                mean(Q16), mean(Q84)), 
          type = c('Min', paste0(QUANTILES*100,'%'), 'Max', 'Weighted Mean',
                   'Weighted Median', '% sig', 'CU_l', 'CU_u')) %>%
  pivot_wider(names_from = 'h', values_from = 'o')
YMax <- max(tableY[3,] %>% as.numeric(), na.rm=T) %>% round(.,2)
YMaxHor <- names(tableY)[which.max(tableY[3,] %>% as.numeric())]
YMaxMean <- max(tableY[6,] %>% as.numeric(), na.rm=T) %>% round(.,2)
# Output Spill-Out Response
tableYOut <- multspillY %>%
  select(region = source, h, value, sig) %>%
  group_by(h) %>%
  reframe(o = c(min(value), quantile(value, QUANTILES), max(value), mean(value), mean(sig)), 
          type = c('Min', paste0(QUANTILES*100,'%'), 'Max', 'Mean', '% sig')) %>%
  pivot_wider(names_from = 'h', values_from = 'o') %>%
  select(`Статистика`=type, `0`, `1`, `3`, `6`, `12`, `24`)
# Output Spill-In Response
tableYIn <- multspillinY %>%
  group_by(h) %>%
  reframe(o = c(min(value), quantile(value, QUANTILES), max(value), mean(value), mean(sig)), 
          type = c('Min', paste0(QUANTILES*100,'%'), 'Max', 'Mean', '% sig')) %>%
  pivot_wider(names_from = 'h', values_from = 'o') %>%
  select(`Статистика`=type, `0`, `1`, `3`, `6`, `12`, `24`)
summ <- readRDS(paste0(path,'/',MODEL,'-summary.RDS'))
n <- str_extract(summ$CD, '\\d+\\.\\d+\\%')
n <- as.numeric(substr(n, 1, nchar(n)-1))
conv <- round((1-n/100)*100,1)