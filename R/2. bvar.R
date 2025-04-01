set.seed(123)
library(readxl)
library(tidyverse)
library(BGVAR)
options(scipen=9999)
saver <- function(object, folder, name) {
  saveRDS(object, paste0(folder, '/', name), compress = FALSE)
}
source('irf.bgvar.R')
environment(irf.bgvar2) <- asNamespace('BGVAR')
assignInNamespace("irf.bgvar", irf.bgvar2, ns = "BGVAR")

# WEIGHT MATRIX ----
W <- readRDS('../Data/W.RDS')
## GRP Weights for Aggregation ----
aggW <- read_excel('../Data/aggregationWeights.xlsx')[1:79,c(2,3)]
# ESTIMATION ----
## Load Data ----
data <- readRDS('../Data/origData.RDS')
variables <- c('exp', 'y', 'pi')
o <- lapply(head(data,-2), function(x) x[,variables])
## Check Order  ----
all(colnames(W) == names(o))
all(colnames(W) == aggW[,1])
## Dominant Unit Settings ----
OIL_COUNTRY <- FALSE # model oil country as separate?
if (OIL_COUNTRY) {
  # Federal block
  RU.weights<-list()
  RU.weights$weights <- aggW %>% pull(w)
  names(RU.weights$weights) <- aggW %>% pull(abb) 
  RU.weights$variables <- c('e', 'expfed', 'i', 'revfed', 'pi', 'y') 
  RU.weights$exo <- c('e', 'expfed', 'i', 'revfed')
  # Oil block
  OC.weights<-list()
  OC.weights$weights <- rep(0,79)
  names(OC.weights$weights) <- aggW %>% pull(abb) 
  OC.weights$variables <- c('poil', 'yrow', 'ffr')
  OC.weights$exo <- c('poil')
  OE.weights <- list(RU=RU.weights, OC=OC.weights)
  o$RU <- data$RU[,c('e', 'expfed', 'i', 'revfed')]
  o$OC <- data$OP[,c('poil','yrow','ffr')]
} else {
  # Federal block
  RU.weights<-list()
  RU.weights$weights <- aggW %>% pull(w)
  names(RU.weights$weights) <- aggW %>% pull(abb) 
  RU.weights$variables <- c('poil', 'expfed', 'e', 'revfed', 'i', 'pi', 'y') 
  RU.weights$exo <- c('poil', 'expfed', 'e', 'revfed', 'i')
  OE.weights <- list(RU=RU.weights)
  o$RU <- cbind(data$OP[,'poil'], data$RU[,c('expfed', 'e', 'revfed', 'i')])
  colnames(o$RU) <- c('poil', 'expfed', 'e', 'revfed', 'i')
}
## Estimate model ----
folder <- '../Results'
nvar <- length(variables)
draws <- 1000
burnin <- 15000
plag <- c(1,1)
prior <- 'MN'
SV <- TRUE
mname <- paste0(nvar, 'v', paste0(plag,collapse='_'), 'p-', prior, '-', ifelse(SV, 'SV', 'NSV'), '-', draws, '-', burnin)
path <- paste0(folder,'/',mname)
dir.create(path, showWarnings = FALSE)
# IRF ANALYSIS ----
if (file.exists(paste0(path,'/',mname,'.RDS')) & !exists('model.1')) {
  model.1 <- readRDS(paste0(path,'/',mname,'.RDS'))
} else if (!file.exists(paste0(path,'/',mname,'.RDS'))) {
  model.1 <- bgvar(Data=o,
                   W = W,
                   draws = draws,
                   burnin = burnin,
                   plag = plag,
                   prior = prior,
                   hyperpara = NULL,
                   SV = T,
                   thin = 1,
                   trend = T,
                   hold.out = 0,
                   eigen = 1.05,
                   expert = list(OE.weights=OE.weights,
                                 cores = 30
                   )
  )
  saver(model.1, path, paste0(mname, '.RDS'))
  msumm <- summary(model.1)
  saver(list(
    CD = msumm$CD$perc,
    res = msumm$res$p.res,
    crosscorr = msumm$cross.corr$res.res,
    crosscorrdata = msumm$cross.corr$dat.res
  ), path, paste0(mname, '-summary.RDS'))
}
## Federal Expenditure Shock ----
### Generalized IRF ----
ident <- 'girf'
shockinfo <- get_shockinfo(ident)
shockinfo$shock <- "RU.expfed"
shockinfo$scale <- 1
i3 <- irf(model.1, n.ahead=24, shockinfo=shockinfo, 
          expert = list(use_R = T, cores = 30))
plot(i3, cumulative = F, resp = 'y')
saver(i3, path, paste0(mname, '-', shockinfo$shock, '_', ident, '.RDS'))
### Recursive Identification ----
ident <- 'chol'
shockinfo <- get_shockinfo(ident)
shockinfo$shock <- "RU.expfed"
shockinfo$scale <- 1
i3 <- irf(model.1, n.ahead=24, shockinfo=shockinfo, 
          expert = list(use_R = T, cores = 30))
plot(i3, cumulative = F, resp = 'y')
saver(i3, path, paste0(mname, '-', shockinfo$shock, '_', ident, '.RDS'))
### Sign Restrictions ----
ident <- 'sign'
shockinfo <- get_shockinfo(ident)
cc <- names(model.1$args$Data)[1:79]
shockinfo <- add_shockinfo(shockinfo, shock="RU.expfed", 
                           restriction=paste0(cc,'.y'), sign=rep(">",length(cc)),
                           horizon=3, prob=0.75, scale=1)
shockinfo <- add_shockinfo(shockinfo, shock="RU.expfed", 
                           restriction='RU.poil', sign='0',
                           horizon=1, prob=1, scale=1)
shockinfo <- add_shockinfo(shockinfo, shock="RU.expfed", 
                           restriction='RU.i', sign='>',
                           horizon=1, prob=1, scale=1)
i3 <- irf(model.1, n.ahead=12, shockinfo=shockinfo, expert = list(MaxTries = 4000,
                                                                  cores = 30,
                                                                  use_R = T))
plot(i3, cumulative = F, resp = 'RU')
saver(i3, path, paste0(mname, '-', shockinfo$shock[1], '_', ident, '.RDS'))
## Regional Expenditure Shocks ----
### Calculate IRFs (Cholesky) ----
shocks <- paste0(names(data) %>% head(-2), '.exp')
ident <- 'chol'
shocks <- shocks[!shocks %in% str_sub(list.files(paste0(path,'/regions_chol')), -15, -10)]
dir.create(paste0(path, '/regions_chol'), showWarnings = FALSE)
for (i in shocks) {
  print(i)
  shockinfo <- get_shockinfo(ident)
  shockinfo$shock <- i
  shockinfo$scale <- 1
  ii <- irf(model.1, n.ahead=24, shockinfo=shockinfo, expert = list(cores = 30, use_R = T))
  gc()
  saver(ii, paste0(path, '/regions_chol'), paste0(mname, '-', shockinfo$shock, '_', ident, '.RDS'))
  rm(ii)
}
f <- list.files(paste0(path, '/regions_chol'), full.names = T)
for (i in 1:length(f)) {
  est <- readRDS(f[i])
  effects <- plyr::adply(est$posterior, .margins = 1:3, .id = c('response','h','shock')) %>%
    separate(response, into = c('region', 'variable'), sep = '\\.') %>%
    group_by(region, variable) %>%
    arrange(h)
  if (i == 1) result <- effects else result <- rbind(result, effects)
}
result <- result %>%
  separate(shock, c("source", "shock"))
saver(result, path, paste0(mname,'-regIRFs-chol.RDS'))
### Calculate IRFs (Sign Restrictions) ----
shocks <- paste0(names(data) %>% head(-2), '.exp')
ident <- 'sign'
shocks <- shocks[!shocks %in% str_sub(list.files(paste0(path,'/regions_sign')), -15, -10)]
dir.create(paste0(path, '/regions_sign'), showWarnings = FALSE)
for (i in shocks) {
  print(i)
  reg <- substr(i, 1, 2)
  shockinfo <- get_shockinfo(ident)
  # # AD shock
  shockinfo <- add_shockinfo(shockinfo, shock=paste0(reg, '.y'),
                           restriction=c(paste0(reg,'.pi'), paste0(reg, '.exp')),
                           sign=c(">",0), horizon=1, scale=1)
  # AS shock
  shockinfo <- add_shockinfo(shockinfo, shock=paste0(reg,'.pi'),
                             restriction=c(paste0(reg,'.y'), paste0(reg, '.exp')),
                             sign=c("<",0), horizon=1, scale=1)
  # Fiscal shock
  shockinfo <- add_shockinfo(shockinfo, shock=i,
                             restriction=paste0(reg, '.y'), sign=">", horizon=3, scale=1)
  ii <- irf(model.1, n.ahead=12, shockinfo=shockinfo, verbose = T, 
            expert = list(cores = 30, use_R = T, MaxTries = 4000))
  gc()
  saver(ii, paste0(path, '/regions_sign'), paste0(mname, '-', i, '_', ident, '.RDS'))
  rm(ii)
}
f <- list.files(paste0(path, '/regions_sign'), full.names = T)
for (i in 1:length(f)) {
  est <- readRDS(f[i])
  effects <- plyr::adply(est$posterior, .margins = 1:3, .id = c('response','h','shock')) %>%
    separate(response, into = c('region', 'variable'), sep = '\\.') %>%
    group_by(region, variable) %>%
    arrange(h)
  if (i == 1) result <- effects else result <- rbind(result, effects)
}
result <- result %>%
  separate(shock, c("source", "shock"))
saver(result, path, paste0(mname,'-regIRFs-sign.RDS'))
## Joint Fiscal Expansion ----
### Top-10 by GRP ----
codes <- read_xlsx('../Technical/regioncodes.xlsx')
explanatory <- readRDS('../Data/explanatory.RDS') %>%
  select(OKATO, OKATO_id, GRP) %>%
  arrange(GRP) %>%
  left_join(codes %>% select(region = Abbrev, name_rus = `Rus Name`, name_eng = Name, OKATO_id)) %>%
  filter(!is.na(region)) %>%
  arrange(GRP) %>%
  left_join(readRDS('../Data/mult_supplementary.RDS') %>% mutate(size = grp * 0.01/exp_level))
poor10 <- explanatory$region[1:10]
poor10s <- explanatory$size[1:10]
rich10 <- explanatory$region[70:79]
rich10s <- explanatory$size[70:79]
ident <- 'chol'
# smallest
shockinfo <- data.frame(shock = paste0(poor10, '.exp'), scale = 1, global = TRUE)
ii <- irf(model.1, n.ahead=24, shockinfo=shockinfo, expert = list(cores = 30, use_R = T))
effects <- plyr::adply(ii$posterior, .margins = 1:3, .id = c('response','h','shock')) %>%
  separate(response, into = c('region', 'variable'), sep = '\\.') %>%
  group_by(region, variable) %>%
  arrange(h)
saver(effects, path, paste0(mname, '-', 'joint_10poorest_', ident, '.RDS'))
# largest
shockinfo <- data.frame(shock = paste0(rich10, '.exp'), scale = 1, global = TRUE)
ii <- irf(model.1, n.ahead=24, shockinfo=shockinfo, expert = list(cores = 30, use_R = T))
effects <- plyr::adply(ii$posterior, .margins = 1:3, .id = c('response','h','shock')) %>%
  separate(response, into = c('region', 'variable'), sep = '\\.') %>%
  group_by(region, variable) %>%
  arrange(h)
saver(effects, path, paste0(mname, '-', 'joint_10richest_', ident, '.RDS'))


