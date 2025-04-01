library(tidyverse)
library(rvest)
library(readxl)
library(tsibble)
library(foreach)
library(doParallel)
select <- dplyr::select
options(scipen=9999)
cluster <- makeCluster(detectCores() - 1)
registerDoParallel(cluster)
## REGIONAL DATA ----
# Population
pop <- read_excel('../Data/population.xlsx') %>%
  pivot_longer(-c(OKATO,OKATO_id), names_to = 'date', values_to = 'population') %>%
  mutate(date = as.Date(as.numeric(date),origin = "1899-12-30"))
# Economic indicators
dbSA <- readRDS('../Data/dbSA.RDS')
# get GRP weights
for (i in 2009:2024) {
  df <- read_xlsx('../Data/VRP_OKVED2_s2016.xlsx', sheet=as.character(i))
  if (i < 2016) {
    pos <- which(df[3,] %in% c('Раздел C',
                               'Раздел D',
                               'Раздел E',
                               'Раздел F',
                               'Раздел G',
                               'Раздел K',
                               'Раздел О',
                               'Раздел O'
    ))
    idx <- which(df[,1] == 'name') + 1
    df <- df[idx:nrow(df),c(1,pos)]
    names(df) <- c('region', 'C', 'D', 'E', 'F', 'G', 'K', 'O')
    df <- df %>%
      mutate(across(-region, ~ as.numeric(.))) %>%
      filter(!is.na(region)) %>%
      mutate(year = i,
             industry = C + D + E,
             constr = F,
             rtrade = G,
             services = K + O) %>%
      select(region, year, industry, constr, rtrade, services)
  }
  if (i >= 2016) {
    pos <- which(df[3,] %in% c('Раздел B',
                               'Раздел C',
                               'Раздел D',
                               'Раздел E',
                               'Раздел F',
                               'Раздел G',
                               'Раздел L',
                               'Раздел M',
                               'Раздел R',
                               'Раздел S'
    ))
    idx <- which(df[,1] == 'name') + 1
    df <- df[idx:nrow(df),c(1,pos)]
    names(df) <- c('region', 'B', 'C', 'D', 'E', 'F', 'G', 'L', 'M', 'R', 'S')
    df <- df %>%
      mutate(across(-region, ~ as.numeric(.))) %>%
      filter(!is.na(region)) %>%
      mutate(year = i,
             industry = B + C + D + E,
             constr = F,
             rtrade = G,
             services = L + M + R + S) %>%
      select(region, year, industry, constr, rtrade, services)
  }
  if (i == 2009) {
    result <- df
  } else {
    result <- rbind(result, df)
  }
}
result <- result %>%
  mutate(sum = industry + constr + rtrade + services,
         industry = industry / sum,
         constr = constr / sum,
         rtrade = rtrade/sum,
         services = services/sum) %>%
  select(-sum, OKATO = region, year) %>%
  left_join(dbSA %>% distinct(OKATO, OKATO_id)) %>%
  pivot_longer(-c(OKATO,OKATO_id,year), values_to = 'w', names_to = 'variable')
# Calculate business activity indicator
activity <- dbSA %>%
  filter(variable %in% c('industry', 'constr', 'rtrade', 'services')) %>%
  filter(type == 'индекс SA (2016 = 100)') %>%
  filter(date >= as.Date('2009-01-01')) %>%
  mutate(year = year(date)) %>%
  left_join(result)
activity <- activity %>%
  mutate(index = w * value) %>%
  select(-c(w,year,type)) 
dd <- activity %>%
  select(OKATO, OKATO_id, date, index) %>%
  pivot_longer(-c(OKATO, OKATO_id, date), names_to = 'variable', values_to = 'value') %>%
  group_by(OKATO, OKATO_id, date, variable) %>%
  summarize(value = sum(value)) 
activity <- rbind(activity %>% select(-index), dd) %>%
  arrange(OKATO, date, variable)
# writexl::write_xlsx(activity %>% pivot_wider(names_from='date'), '123.xlsx')
activity %>%
  filter(variable == 'index') %>%
  ggplot() +
  geom_line(aes(x=date, y=value)) +
  facet_wrap(~OKATO_id, scales='free_y') +
  theme(legend.position = 'none')
activity <- activity %>%
  filter(variable == 'index')
# Region codes for Matlab
codes <- read_xlsx('../Technical/regioncodes.xlsx')
# Inflation
infl <- dbSA %>%
  filter(variable == 'inflation' & type == 'индекс SA (2016 = 100)')
# Budget expenditure & revenue data + adjustment for inflation
budgetData <- readRDS('../Data/filteredFullSA.RDS') %>%
  filter(type %in% c('exp', 'rev'),
         budget == 'consolMSA') %>%
  ungroup() %>%
  select(date, region, OKATO_id, variable = type, value) %>%
  # filter(date >= as.Date('2014-12-01')) %>%
  left_join(infl %>% select(OKATO_id, date, cpi = value)) %>%
  left_join(pop %>% select(-OKATO)) %>%
  fill(population) %>%
  mutate(value = value/cpi*100, # million 2016 rubles
         population = population/10^6, # million people,
         value = value / population # 2016 rubles per person
  ) %>%
  left_join(codes %>% select(`Rus Name`, OKATO_id)) %>%
  select(OKATO = `Rus Name`, OKATO_id, date, variable, value) 
# calculate budget deficit
budgetData <- budgetData %>%
  spread(variable, value) %>%
  mutate(deficit = exp - rev) %>%
  gather(variable, value, -c(OKATO,OKATO_id,date))
# Calculate regional variables
vars <- c('inflation', 'unemprate', 'wager')
result <- data.frame()
for (i in vars) {
  data <- dbSA %>%
    filter(variable %in% i) %>%
    filter(type %in% c('индекс SA (2016 = 100)', 'тыс. чел. SA', 'руб. 2016 г. SA', '% SA')) %>%
    group_by(OKATO_id)
  if (i == 'inflation') {
    data <- data %>%
      mutate(
        value = log(value/lag(value))*100
      )
  } else if (i == 'unemprate') {
    data <- data
  } else if (i == 'wager') {
    data <- data %>%
      mutate(value = log(value/lag(value))*100)
  }
  result <- rbind(result, data)
}
result %>%
  left_join(codes %>% select(`Short Name`, OKATO_id)) %>%
  filter(!is.na(`Short Name`)) %>%
  ggplot() +
  geom_line(aes(x=date, y = value, color = OKATO)) +
  facet_wrap(~variable, scale='free_y') +
  theme(legend.position = 'none')
result <- result %>%
  select(OKATO, OKATO_id, date, variable, value)
# Plug back in & add budget
econDataOriginal <- result %>%
  rbind(., budgetData) %>%
  rbind(., activity %>% mutate(value = log(value)*100)) %>%
  left_join(codes %>% select(`Short Name`, OKATO_id)) %>%
  filter(!is.na(`Short Name`))
econDataOriginal %>%
  ggplot() +
  geom_line(aes(x=date, y=value, color=OKATO_id)) +
  facet_wrap(~variable, scales='free_y') +
  theme(legend.position = 'none')
## FEDERAL DATA ----
# RUONIA interest rate
ruonia <- read_excel('../Data/ruonia.xlsx') %>%
  mutate(date = as.Date(DT)) %>%
  select(date, ruonia=ruo) %>%
  group_by(date = yearmonth(date) %>% as.Date(.)) %>%
  summarize(ruo = mean(ruonia, na.rm=T))
# MIACR
miacr <- read_excel('../Data/miacr.xlsx')
# USD/RUB
exrate <- read_excel('../Data/exrate.xls')
# Federal Budget
fedbudget <- read_excel('../Data/fedbudget.xlsx') %>%
  left_join(pop %>% filter(OKATO_id == '643') %>% select(date, population))  %>%
  fill(population) %>%
  left_join(infl %>% select(OKATO_id, date, cpi = value) %>% filter(OKATO_id == '643')) %>%
  mutate(fedexpR = fedbudget/cpi*100, # million 2016 rubles
         fedrevR = revenue/cpi*100,   # million 2016 rubles
         population = population/10^6, # million people,
         fedexpR = fedexpR / population, # 2016 rubles per person
         fedrevR = fedrevR / population # 2016 rubles per person
  )
unsa <- ts(fedbudget$fedexpR, frequency = 12, start = c(year(first(fedbudget$date)), month(first(fedbudget$date))))
sa <- RJDemetra::x13(unsa)
fedbudget$fedexpR <- sa$final$series[,2]
unsa <- ts(fedbudget$fedrevR, frequency = 12, start = c(year(first(fedbudget$date)), month(first(fedbudget$date))))
sa <- RJDemetra::x13(unsa)
fedbudget$fedrevR <- sa$final$series[,2]
# Monthly GDP proxy
ivbo <- read_excel('../Data/ivbo.xlsx')
unsa <- ts(ivbo$base, frequency = 12, start = c(year(first(ivbo$date)), month(first(ivbo$date))))
sa <- RJDemetra::x13(unsa)
ivbo$ivbo <- sa$final$series[,2]
# All global variables
fedData <- exrate %>%
  left_join(ruonia) %>%
  left_join(miacr) %>%
  left_join(fedbudget %>% select(date, fedexpR, fedrevR)) %>%
  left_join(ivbo %>% select(date, ivbo))
fedData <- fedData %>%
  pivot_longer(-date, names_to = 'variable')
fedData %>%
  ggplot() +
  ggtitle('Российская Федерация') +
  geom_line(aes(x=date, y=value)) +
  facet_wrap(~variable, scales='free_y') +
  theme(legend.position = 'none',
        text = element_text(size=25))
# ggsave('../Figures/RawData/Rus.png', scale=2)
## OIL PRICE MODEL ----
# Urals, FFR and global activity
urals <- read_excel('../Data/urals.xlsx')
urals <- urals %>%
  mutate(poil = log(urals)) %>%
  select(date, poil, ffr = sffr, yrow)
## PACK TOGETHER ----
# ORIGINAL LEVELS
SAMPLESTART <- as.Date('2009-11-01')
SAMPLEEND <- as.Date('2024-06-01')
l <- list()
econData2 <- econDataOriginal %>%
  left_join(codes %>% select(`Short Name`, Abbrev)) %>%
  arrange(`Short Name`) %>%
  select(-`Short Name`) %>%
  mutate(variable = case_when(
    variable == 'unemprate' ~ 'u',
    variable == 'wager' ~ 'w',
    variable == 'inflation' ~ 'pi',
    variable == 'index' ~ 'y',
    TRUE ~ variable
  )) %>%
  filter(date >= as.Date('2009-01-01'))
# corrections for negatives (2 values over 2009-2024)
ff <- econData2 %>% filter(variable %in% c('exp','rev') & value < 0)
for (i in 1:nrow(ff)) {
  idx <- econData2 %>% with(which(OKATO_id == ff$OKATO_id[i] & date == ff$date[i] & variable == ff$variable[i]))
  econData2[idx,'value'] <- (econData2[idx-1,'value'] + econData2[idx+1,'value'])/2
}
for (j in unique(econData2$Abbrev)) {
  sub <- econData2 %>%
    filter(Abbrev == j) %>%
    mutate(value = case_when(
      variable %in% c('exp','rev') ~ log(value),
      TRUE ~ value
    )) %>%
    pivot_wider(names_from = 'variable', values_from = 'value') %>%
    filter(date >= SAMPLESTART & date<= SAMPLEEND) %>%
    ungroup() %>%
    select(-c(Abbrev, OKATO, OKATO_id)) %>%
    arrange(date)
  df_ts <- ts(sub[,-1], start=c(year(SAMPLESTART), month(SAMPLESTART)), freq=12)  
  l[[j]] <- df_ts
}
sub <- fedData %>% 
  filter(!variable == 'ruo') %>%
  group_by(variable) %>%
  mutate(value = case_when(
    variable %in% c('ruo', 'miacr') ~ value,
    TRUE ~ log(value)
  )) %>% 
  mutate(variable = case_when(
    variable == 'ruo' ~ 'i',
    variable == 'miacr' ~ 'i',
    variable == 'urals' ~ 'poil',
    variable == 'fedexpR' ~ 'expfed',
    variable == 'fedrevR' ~ 'revfed',
    variable == 'usdrub' ~ 'e',
    variable == 'ivbo' ~ 'gdp'
  )) %>%
  pivot_wider(names_from = 'variable', values_from = 'value') %>%
  filter(date >= SAMPLESTART & date<= SAMPLEEND)
df_ts <- ts(sub[,-1], start=c(year(SAMPLESTART), month(SAMPLESTART)), freq=12)  
l[['RU']] <- df_ts
plot(l$RU)
urals <- urals %>% 
  filter(date >= SAMPLESTART & date<= SAMPLEEND)
df_ts <- ts(urals[,-1], start=c(year(SAMPLESTART), month(SAMPLESTART)), freq=12)  
l[['OP']] <- df_ts
plot(l$OP)
saveRDS(l, '../Data/origData.RDS')
# Compute Weight Matrix
codes <- read_xlsx('../Technical/regioncodes.xlsx')
## Load data ----
w <- read_excel('../Data/W.xls')
n1 <- w[8:94,2] %>% pull()
n2 <- w[5,4:90] %>% as.character()
regs <- w[7:94,1:2] %>%
  rename_with(~c('region', 'code'))
## Remove Unreported and International ----
w2 <- w[8:94,4:90] %>% as.matrix() %>% apply(., 2, as.numeric)
colnames(w2) <- n2
rownames(w2) <- n1
## Rearrange ----
# Ex to COL, Im to ROW
# Excluded regions:
# 11100000000 - Ненецкий АО
# 11001000000 - Архангельская область без АО
# 71100000000 - Ханты-Манскийский АО
# 71140000000 - Ямало-Ненецкий АО
# 71001000000 - Тюменская область без АО
# 77000000000 - Чукотский АО
# 35000000000 - Республика Крым
# 67000000000 - Севастополь
except <- c('11100000000', '11001000000',
            '71100000000', '71140000000', '71001000000', 
            '77000000000', 
            '35000000000',
            '67000000000')
w2 <- w2[!rownames(w2) %in% except,
         !colnames(w2) %in% except]
## Get Trade Shares ----
ExSh <- (colSums(w2)-diag(w2))/diag(w2)
ImSh <- (rowSums(w2)-diag(w2))/diag(w2)
trade <- data.frame(exsh = ExSh,
                    imsh = ImSh) %>%
  rownames_to_column(var = 'OKATO_id')
## Add Codes ----
rc <- codes %>%
  select(OKATO_id, `Short Name`) %>%
  deframe()
colnames(w2) <- str_replace_all(colnames(w2), rc)
rownames(w2) <- str_replace_all(rownames(w2), rc)
w2 <- w2[order(rownames(w2)), order(colnames(w2))]
## Compute W ----
m <- w2+t(w2) - 2*diag(diag(w2), nrow=nrow(w2), ncol=ncol(w2))
rc2 <- codes %>%
  filter(!is.na(Abbrev)) %>%
  arrange(`Short Name`) %>%
  pull(Abbrev)
W <- prop.table(m, 1)
colnames(W) <- rc2
rownames(W) <- rc2
saveRDS(W, '../Data/W.RDS')
# Calculate ratios for multiplier computation
rr <- econData2 %>%
  filter(date >= SAMPLESTART & date<= SAMPLEEND) %>%
  filter(variable == 'exp') %>%
  pivot_wider(names_from = 'variable', values_from = 'value') %>%
  mutate(exp_level=exp) %>%
  group_by(OKATO_id) %>%
  summarize(exp_level = mean(exp_level))
grp <- read_xlsx('../Data/VRP.xlsx') %>%
  pivot_longer(-OKATO, names_to = 'year') %>%
  mutate(year = as.numeric(year),
         value = as.numeric(value)) %>%
  filter(!is.na(OKATO)) %>%
  left_join(codes %>% select(OKATO = `Rus Name`, OKATO_id)) %>%
  left_join(infl %>%
              select(OKATO_id, date, value) %>%
              group_by(OKATO_id, year = year(date)) %>%
              summarize(cpi = mean(value))) %>%
  mutate(grp = value / cpi * 100) %>%
  filter(!is.na(cpi)) %>%
  select(OKATO_id, year, grp) %>%
  filter(year >= 2009) %>%
  group_by(OKATO_id) %>%
  summarize(grp = mean(grp))
rr <- rr %>%
  left_join(grp)
saveRDS(rr, '../Data/mult_supplementary.RDS')  
# Compute explanatory variables
data <- readRDS('../Data/dbSA.RDS')
ltight <- data %>%
  filter(variable %in% c('unempl', 'ldemand', 'unemprate'),
         type %in% c('тыс. чел.', '%')) %>%
  group_by(OKATO, OKATO_id, year = year(date), variable) %>%
  summarize(value = mean(value)) %>%
  pivot_wider(names_from = 'variable', values_from = 'value') %>%
  mutate(tightness = unempl/ldemand) %>%
  filter(year >= 2017 & year <= 2023) %>%
  group_by(OKATO, OKATO_id) %>%
  summarize(tightness = mean(tightness), unemprate = mean(unemprate))
debtData <- readRDS('../Data/debtFull.RDS') %>%
  group_by(OKATO_id, year = year(date)) %>%
  summarize(debt = mean(debt, na.rm=T), debtM = mean(debt - budget, na.rm=T))
grp <- read_xlsx('../Data/VRP.xlsx', sheet = 'Total') %>%
  pivot_longer(-OKATO, names_to = 'year') %>%
  mutate(year = as.numeric(year),
         grp = as.numeric(value)/10^3) %>%
  filter(!is.na(OKATO)) %>%
  left_join(codes %>% select(OKATO = `Rus Name`, OKATO_id)) %>%
  select(-value)
debtData <- debtData %>%
  left_join(grp) %>%
  mutate(debtgrp = debtM/grp * 100) %>%
  group_by(OKATO) %>%
  summarize(debtgrp = mean(debtgrp, na.rm=T))
shares <- readRDS('../Data/spendingShares.RDS') %>%
  group_by(OKATO_id, name) %>%
  summarize(value = mean(consShare, na.rm = T)*100) %>%
  pivot_wider(names_from = 'name') %>%
  rename_with(~c('OKATO_id', 'utilities', 'healthcare', 'culture', 'security', 'defense', 'economy', 
                 'education', 'state', 'social', 'sport'))
# HTM calculation
data <- read.csv('../Data/other/5th_wave_hh_230323.csv', fileEncoding = 'windows-1251')
data2 <- read.csv('../Data/other/5thwave_ind_240323.csv', fileEncoding = 'windows-1251')
filt <- data %>%
  select(id_h, psu, hhwgt,
         labinc = h37_n_s_r, reginc = h47_n_s_r, totinc = h48_n_s_r, 
         saving = h34_n_s_r,
         propvalue1 = a46a_s_r, mortgage1 = a66a_s_r, # основное место жительства
         propvalue2 = b136a1_s_r, propvalue3 = b136a2_s_r, propvalue4 = b136a3_s_r, mortgage2 = b1_47a_s_r,  # квартиры
         propvalue5 = b246a1_s_r, propvalue6 = b246a2_s_r, propvalue7 = b246a3_s_r, mortgage3 = b2_57a_s_r,  # дома
         propvalue8 = b337a1_s_r, propvalue9 = b337a2_s_r, propvalue10 = b337a3_s_r, mortgage4 = b3_48a_s_r, # земля
         propvalue11 = b4_35a_s_r # гаражи
  ) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(prop = select(., starts_with("propvalue")) %>% rowSums(),
         mortgage = select(., starts_with("mortgage")) %>% rowSums()) %>%
  select(-matches('mortgage\\d')) %>%
  select(-matches(regex('propvalue\\d'))) 
filt2 <- data2 %>%
  select(id_h, idind, 
         cc1 = p615a1_s_r, cc2 = p615a2_s_r, cc3 = p615a3_s_r, # кредитная карта
         climit = c2_4,
         # clim1 = p610a1_s_r, clim2 = p610a2_s_r, clim3 = p610a3_s_r, # лимиты 3 наиболее используемых КК
         microfin = c3_3n_s_r, # микрофинансовые организации, ломбарды
         potreb1 = p911a1_s_r, potreb2 = p911a2_s_r, potreb3 = p911a3_s_r, # потребительские кредиты
         selfconstr = c1_29, # не планирует брать, т. к. уверен, что откажут в выдаче кредита
         refusal = c1_12 # отказы в выдаче кредита в последние 2 года
  ) %>%
  mutate(across(starts_with('cc'), ~ as.numeric(.)),
         across(starts_with('clim'), ~ as.numeric(.)),
         climit = as.numeric(climit),
         selfconstr = selfconstr == 'УВЕРЕН(А), ЧТО МНЕ ОТКАЖУТ',
         refusal = refusal == 'СЛУЧАЛОСЬ'
  ) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(cc = cc1 + cc2 + cc3,
         potreb = potreb1 + potreb2 + potreb3,
         selfconstr = as.numeric(selfconstr)*1,
         refusal = as.numeric(refusal)*1) %>%
  select(id_h, idind, cc, climit, potreb, microfin, selfconstr, refusal) %>%
  group_by(id_h) %>%
  summarize(cc = sum(cc),
            potreb = sum(potreb),
            microfin = sum(microfin),
            climit = sum(climit),
            selfconstr = sum(selfconstr) >= 1,
            refusal = sum(refusal) >= 1)
res <- filt %>%
  left_join(filt2, by = 'id_h') %>%
  mutate(lwealth = saving - cc - potreb - microfin,
         illwealth = prop - mortgage) %>%
  mutate(htm = (lwealth>=0 & lwealth <= totinc/2) | (lwealth<=0 & lwealth <= totinc/2 - climit),
         whtm = htm & illwealth > 0,
         phtm = htm & illwealth <=0
  ) 
sw <- res %>%
  mutate(psu = case_when(
    psu == ' Краснодарский край ' ~ 'Краснодарский край',
    TRUE ~ psu
  )) %>%
  group_by(OKATO=psu) %>%
  mutate(w = hhwgt/sum(hhwgt)) %>%
  summarize(htm = sum(htm)/n(),
            whtm = sum(whtm)/n(),
            phtm = sum(phtm)/n()) %>%
  mutate(OKATO = case_when(
    OKATO == 'Республика Кабардино-Балкария' ~ 'Кабардино-Балкарская Республика',
    OKATO == 'Республика Татарстан' ~ 'Республика Татарстан (Татарстан)',
    OKATO == 'Республика Удмуртия' ~ 'Удмуртская Республика',
    OKATO == 'Республика Чувашия' ~ 'Чувашская Республика - Чувашия',
    OKATO == 'г. Москва' ~ 'Город Москва столица Российской Федерации город федерального значения',
    OKATO == 'г. Санкт-Петербург' ~ 'Город Санкт-Петербург город федерального значения',
    TRUE ~ OKATO    
  ))
explanatory <- read_xlsx('../Data/explanatory.xlsx') %>%
  left_join(ltight) %>%
  select(-region) %>%
  left_join(trade) %>%
  left_join(debtData) %>%
  left_join(shares) %>%
  left_join(sw)
saveRDS(explanatory, '../Data/explanatory.RDS')
