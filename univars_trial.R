library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
# Set random seed to ensure reproducible results
rand_seed = 9430568

baci_diffs <- function(mat){
  data.frame(delta_baci1_low = 
               mat[,combos$control_lowi[1]] - mat[,combos$impact_lowi[1]] -
               (mat[,combos$control_lowi[2]] - mat[,combos$impact_lowi[2]]),
             delta_baci2_low = 
               mat[,combos$control_lowi[1]] - mat[,combos$impact_lowi[1]] -
               (mat[,combos$control_lowi[3]] - mat[,combos$impact_lowi[3]]),
             delta_baci1_hi =
               mat[,combos$control_hii[1]] - mat[,combos$impact_hii[1]] -
               (mat[,combos$control_hii[2]] - mat[,combos$impact_hii[2]]),
             delta_baci2_hi = mat[,combos$control_hii[1]] - mat[,combos$impact_hii[1]] -
               (mat[,combos$control_hii[3]] - mat[,combos$impact_hii[3]]))
}

## Load data: ultimately from OSF repository
# library(osfr); library(dplyr)
# if(!"data" %in% dir()){system("mkdir data")}
# if(!"wq_data_compiled.xlsx" %in% dir("data")){
# wq_files <- osf_retrieve_node("4ywvq") %>% osf_ls_files()
# osf_download(wq_files[wq_files$name == "wq_data_compiled.xlsx",], path = "data")
# }
# data_for_model.xlsx compiled in urban_riffle_exp_data_compilation.R
dat_file <- 
  "data/urban_riffle_experiment_data_for_model.xlsx"
sites <- data.frame(readxl::read_excel(dat_file, sheet = "sites"), 
                    stringsAsFactors = FALSE)
samples <- data.frame(readxl::read_excel(dat_file, sheet = "samples"), 
                      stringsAsFactors = FALSE)
# one field didn't translate as numeric
samples$vel_m_s_mean <- as.numeric(samples$vel_m_s_mean)
biota <- data.frame(readxl::read_excel(dat_file, sheet = "biota"), 
                    stringsAsFactors = FALSE)
taxa <- data.frame(readxl::read_excel(dat_file, sheet = "taxa"), 
                   stringsAsFactors = FALSE)
higher_taxa <- data.frame(readxl::read_excel(dat_file, sheet = "higher_taxa"), 
                          stringsAsFactors = FALSE)

samples$seg <- substr(samples$old_samplecode,nchar(samples$old_samplecode)-1,
                      nchar(samples$old_samplecode)-1)
# The dataset contains samples from the segments upstream and downstream of the 
# riffle (putative or real) in each site for trip 1-4 (out of 6).  Given the 
# small effects in the 'riffle' segments (M), theu U and L segments were not 
# included in the analysis.
samples <- samples[samples$seg == "M",]

sites <- sites[order(sites$exp_treatment,sites$ai),]
sites$site_no <- 1:nrow(sites)
samples$site_no <- sites$site_no[match(samples$sitecode,sites$sitecode)]
samples$sample <- substr(samples$old_samplecode,1,nchar(samples$old_samplecode)-1)
sample_nos <- data.frame(sample = unique(samples$sample))
sample_nos$sample_no <- 1:nrow(sample_nos)
samples$sample_no <- sample_nos$sample_no[match(samples$sample, sample_nos$sample)]
samples$t <- as.numeric(substr(samples$old_samplecode, 1,1))

# Assemble fixed predictors (a1,a2, ci, ba1ci, ba2ci, i, ba1cii, ba2cii, spring) into a matrix, u
samples$ba <- as.numeric(as.numeric(substr(samples$old_samplecode,1,1)) %in% c(3,4)) 
samples$ba[as.numeric(as.numeric(substr(samples$old_samplecode,1,1))) %in% c(5,6)] <- 2
samples$ba <- factor(samples$ba)
# a1 and a2 are the two after periods treated as categories with b as a reference
samples$ci <- as.numeric(sites$exp_treatment[match(samples$sitecode,sites$sitecode)] == "riffle") 
# 0 = control, 1 = riffle
samples$ai <- sites$ai[match(samples$sitecode, sites$sitecode)]
i_scaled <- scale(log10(samples$ai*100 + 0.1))
samples$i <- as.vector(i_scaled)
samples$range_phi <- samples$max_phi - samples$min_phi

unscale <- function(y_scaled_transformed, x, scaled = TRUE, log = TRUE, log_add = 0){  
  #y_scaled_transformed is modelled output of variable x calculated on scaled 
  # (if scale = TRUE), and logged (if log = TRUE) data, x
  #x is the original vector of raw modelled data before transformation (if log = TRUE) and scaling
  # log_add is the amount added before logging to avoid log(0)
  if(log){
  x_transformed <- log(x + log_add)
  }else{
      x_transformed <- x
  }
  if(scaled){
  scale_x <- scale(x_transformed) 
  y_transformed <- y_scaled_transformed * attr(scale_x, 'scaled:scale') + 
                                  attr(scale_x, 'scaled:center')
  }else{
    y_transformed <- y_scaled_transformed
  }
  if(log){
    y <- exp(y_transformed) - log_add
  }else{
    y <- y_transformed
  }
  list(y_transformed = y_transformed, 
       y = y)
}

ys <- list(depth_m_var = as.vector(scale(log(samples$depth_m_var + 1.75e-04))),
           # min(samples$depth_m_var[!is.na(samples$depth_m_var) & samples$depth_m_var > 0]) # 0.0007442
           depth_m_mean = log(samples$depth_m_mean),
           vel_m_s_var = as.vector(scale(log(samples$vel_m_s_var + 0.0007442))),
           # min(samples$vel_m_s_var[!is.na(samples$vel_m_s_var) & samples$vel_m_s_var > 0]) # 0.0007442
           vel_m_s_mean = as.vector(scale(log(samples$vel_m_s_mean))), 
           cpom_g = as.vector(scale(log(samples$cpom_g + 0.001))),
           # min(samples$cpom_g[!is.na(samples$cpom_g) & samples$cpom_g > 0]) # 0.001
           cpomw_g = as.vector(scale(log(samples$cpom_g + samples$wood_g + 0.001))),
           fpom_g = as.vector(scale(log(samples$fpom_g  + 0.0094382))),
           # min(samples$fpom_g[!is.na(samples$fpom_g) & samples$fpom_g > 0]) # 0.0094382
           algae_g = as.vector(scale(log(samples$algae_g + 0.001))), 
           # min(samples$algae_g[!is.na(samples$algae_g) & samples$algae_g > 0]) # 0.001
           macrophyte_g = as.vector(scale(log(samples$macrophyte_g + 5e-04))),
           # min(samples$macrophyte_g[!is.na(samples$macrophyte_g) & samples$macrophyte_g > 0]) # 5e-04
           min_phi = as.vector(scale(samples$min_phi)),
           max_phi = as.vector(scale(samples$max_phi)),
           med_phi = as.vector(scale(samples$med_phi)))

y_transform_pars <- data.frame(stat = names(ys)[1:12],
                               scaled = c(TRUE,FALSE,rep(TRUE,10)), 
                               log = c(rep(TRUE,9),rep(FALSE,3)), 
                               log_add = c(1.75e-04,0,0.0007442,0,0.001,
                                           0.001,0.0094382,0.001,5e-04,0,0,0))

mod <- cmdstan_model("normal_rand_sa_si_fixedmatrix_baci.stan", pedantic = TRUE) 

univar_beta_summary <- 
  data.frame(model = NA, param = NA, variable = NA, mean = NA, median = NA,
             sd = NA, mad = NA, q5 = NA, q95 = NA, rhat = NA, ess_bulk = NA, 
             ess_tail = NA)[0,]
univar_diagnostics <-
  data.frame(model = NA,divergences = NA,max_treedepth_exceedences = NA,
             min_bfmi = NA,max_bfmi = NA,max_rhat = NA,min_ess_bulk = NA,
             min_ess_tail = NA)[0,]
gofs <-data.frame(model = NA, oe_r_squ = NA, oe_slope = NA, pred_95 = NA)[0,]
param_summary <- list()
cfs <- list()

for(i in 1:12){
# remove NAs, and check if data covers trips 5 and 6
  yi <- names(ys)[i]
samplesi <- samples[!is.na(samples[ifelse(yi == "cpomw_g","cpom_g",yi)]),]
if(length(unique(samplesi$sitecode)) != nrow(sites)) stop()
# recode ba for cases where there are no data for ba 2 (trips 5 and 6)
samplesi$ba <- factor(samplesi$ba)
# OM data for trips 5-6 unreliable
if(grepl("_g",yi)){
  samplesi <- samplesi[samplesi$t < 5,]
}
u <- model.matrix(~ ba + ci + ba:ci + i + ba:ci:i, data = samplesi)
# and for sample no (as there a variable missing samples)
sample_nos <- data.frame(sample = unique(samplesi$sample))
sample_nos$sample_no <- 1:nrow(sample_nos)
samplesi$sample_no <- sample_nos$sample_no[match(samplesi$sample, sample_nos$sample)]

# Data list for Stan
sdata <- list(n_obs = nrow(samplesi),
              n_site = nrow(sites),
              n_sample = nrow(sample_nos),
              n_pred = ncol(u),
 #             n_t = max(samplesi$t),
              site_no = samplesi$site_no,
              samp_no = samplesi$sample_no,
#              t_no = samplesi$t,
              u = u,
              y = ys[[i]][match(samplesi$smpcode,samples$smpcode)]
)

nt <- 4; nc <- 4
ni <- rep(2000,12) # c(6000,2000,6000,6000,6000,3000,3000,3000,3000,4000,4000,4000)
nb <- rep(1000,12) # c(3000,1000,3000,3000,3000,2000,2000,2000,2000,2500,2500,2500)
# even with standardized y, 6000 iters, vel_m_s_var has low BFMI (0.035-0.076) and 0% divergences    
# sigma_sa/a_sa, seemingly the main problem: even if prior for sigma_sa is set to exponential(3) (as it seems to sit near 0)
# bfmi remains low.  All main effects are well sampled
# similar story for vel_m_s_mean

ad <- rep(0.8,12) # c(0.9,0.8,0.99,0.95,0.95,0.99,0.99,0.95,0.98,0.85,0.95,0.95)
stanfit_i <- mod$sample(data = sdata,
                        seed = rand_seed, chains = nc,
                        parallel_chains = nc, iter_warmup = nb[i],
                        adapt_delta = ad[i],
                        max_treedepth = ifelse(i %in% c(1,2,3,5,6,7,9,11,12),15,10),
                        iter_sampling = (ni - nb)[i], refresh = 100)
# #  save csv files rather than the model object to use less RAM
stanfit_i$save_output_files(
  dir = "~/uomShare/wergStaff/ChrisW/git-data/urban_riffle_experiment/model_fits/",
  basename = paste0("fit_riffle_baci_", names(ys)[i]), timestamp = FALSE, random = FALSE)
stanfit_i$sampler_diagnostics()
saveRDS(stanfit_i, file = paste0("~/uomShare/wergStaff/ChrisW/git-data/urban_",
                                 "riffle_experiment/model_fits/fit_riffle_baci_",
                                 names(ys)[i], ".rds"))
stanfit_i$diagnostic_summary()
# EBFMI 0.273, 0.289, 0.290, 0.309, zero divergences, zero max treedepth reached.
# # The above three steps required < 500 Mb RAM
summ <- stanfit_i$summary() # This took ~2h and needed >40 Gb RAM
min(summ$ess_bulk,na.rm=TRUE) 
min(summ$ess_tail,na.rm=TRUE) 
max(summ$rhat,na.rm=TRUE)  

univar_beta_summary <- rbind(univar_beta_summary, 
                             data.frame(model = names(ys)[i], 
                                  param = names(as.data.frame(u)),
                                  as.data.frame(summ[grep("gamma",summ$variable),])))
univar_diagnostics <- rbind(univar_diagnostics, data.frame(model = names(ys)[i],
                                 divergences = sum(stanfit_i$diagnostic_summary()$num_divergent),
                                 max_treedepth_exceedences = sum(stanfit_i$diagnostic_summary()$num_max_treedepth),
                                 min_bfmi = min(stanfit_i$diagnostic_summary()$ebfmi),
                                 max_bfmi = min(stanfit_i$diagnostic_summary()$ebfmi),
                                 max_rhat = max(summ$rhat,na.rm=TRUE),
                                 min_ess_bulk = min(summ$ess_bulk,na.rm=TRUE),
                                 min_ess_tail = min(summ$ess_tail,na.rm=TRUE)))

drawsi <- as.data.frame(stanfit_i$draws(format = "df", variables = c("a_si","a_sa","gamma",  #"a_t",
                                                                     "sigma_si","sigma_sa","sigma")))  #,"sigma_t"
predx <- unique(cbind(samplesi[c("sample","site_no","sample_no","t")],u))
temp <- aggregate(sdata$y, by = list(sample = samplesi$sample),FUN = mean)
predx$x <- temp$x[match(predx$sample, temp$sample)]
  if("ba2" %in% names(as.data.frame(u))){
    predy_draws <-
      drawsi[grep("a_si",names(drawsi))][match(predx$site_no, 1:sdata$n_site)] +
      drawsi[grep("a_sa",names(drawsi))][match(predx$sample_no, 1:sdata$n_sample)] +
#      drawsi[grep("a_t",names(drawsi))][match(predx$t, 1:sdata$n_t)] +
      drawsi[,"gamma[1]"] %*% t(predx$`(Intercept)`) +
      drawsi[,"gamma[2]"] %*% t(predx$ba1) +
      drawsi[,"gamma[3]"] %*% t(predx$ba2) +
      drawsi[,"gamma[4]"] %*% t(predx$ci) +
      drawsi[,"gamma[5]"] %*% t(predx$i) + 
      drawsi[,"gamma[6]"] %*% t(predx$`ba1:ci`) +
      drawsi[,"gamma[7]"] %*% t(predx$`ba2:ci`) +
      drawsi[,"gamma[8]"] %*% t(predx$`ba0:ci:i`) + 
      drawsi[,"gamma[9]"] %*% t(predx$`ba1:ci:i`) +
      drawsi[,"gamma[10]"] %*% t(predx$`ba2:ci:i`) 
  }else{
    predy_draws <-
      drawsi[grep("a_si",names(drawsi))][match(predx$site_no, 1:sdata$n_site)] +
      drawsi[grep("a_sa",names(drawsi))][match(predx$sample_no, 1:sdata$n_sample)] +
#      drawsi[grep("a_t",names(drawsi))][match(predx$t, 1:sdata$n_t)] +
      drawsi[,"gamma[1]"] %*% t(predx$`(Intercept)`) +
      drawsi[,"gamma[2]"] %*% t(predx$ba1) +
      drawsi[,"gamma[3]"] %*% t(predx$ci) +
      drawsi[,"gamma[4]"] %*% t(predx$i) + 
      drawsi[,"gamma[5]"] %*% t(predx$`ba1:ci`) +
      drawsi[,"gamma[6]"] %*% t(predx$`ba0:ci:i`) +
      drawsi[,"gamma[7]"] %*% t(predx$`ba1:ci:i`)
   }
names(predy_draws) <- predx$sample

param_summary[[i]] <- 
  data.frame(param = names(as.data.frame(u)),
             mean = unlist(lapply(drawsi[grep("gamma",names(drawsi))], mean)),
             t(as.data.frame(lapply(drawsi[grep("gamma",names(drawsi))], 
                                    quantile, probs = c(0.025,0.05,0.5,0.95,0.975)))))
names(param_summary)[i] <- names(ys)[i]

meansi <- apply(predy_draws, 2, FUN = mean)
cisi <- apply(predy_draws, 2, quantile, probs = c(0.025,0.975))
lmi <- lm(meansi ~ predx$x)
# plot(predx$x,meansi, xlab="Observed",ylab="Predicted",
#      ylim = range(cisi))
# abline(0,1); abline(lmi,lty = 3)
# for(i in 1:nrow(predx)){
#   lines(rep(predx$x[i],2), cisi[,i])
# }
gofsi <- data.frame(model = names(ys)[i],
                   oe_r_squ = summary(lmi)$adj.r.squared,
                   oe_slope = lmi$coefficients[2],
                   pred_95 = 0)
for(j in 1:nrow(predx)){
  gofsi$pred_95 <- gofsi$pred_95 + 
    as.numeric(findInterval(predx$x[j], cisi[,j]) == 1)
  #This adds one if the observed value is between the two CIs
}
gofsi$pred_95 <- gofsi$pred_95/nrow(predx)
gofs <- rbind(gofs, gofsi)

if("ba2" %in% names(as.data.frame(u))){
  predx_cf <- expand.grid(intercept = 1, ba = factor(c(0,1,2)), ci = c(0,1), 
                        i = seq(min(u[,5]),max(u[,5]),length=10))
}else{
  predx_cf <- expand.grid(intercept = 1, ba = factor(c(0,1)), ci = c(0,1), 
                          i = seq(min(u[,5]),max(u[,5]),length=10))
}
predx_cf <- data.frame(model.matrix(~ ba + ci + ba:ci + i + ba:ci:i, 
                                    data = predx_cf))
names(predx_cf) <- names(as.data.frame(u))
predx_cf$ai <- 10^(predx_cf$i * attr(i_scaled, 'scaled:scale') + 
                     attr(i_scaled, 'scaled:center')) - 0.1
  if("ba2" %in% names(as.data.frame(u))){
    cf_drawsi <- 
      drawsi[,"gamma[1]"] %*% t(predx_cf$`(Intercept)`) +
      drawsi[,"gamma[2]"] %*% t(predx_cf$ba1) +
      drawsi[,"gamma[3]"] %*% t(predx_cf$ba2) +
      drawsi[,"gamma[4]"] %*% t(predx_cf$ci) +
      drawsi[,"gamma[5]"] %*% t(predx_cf$i) + 
      drawsi[,"gamma[6]"] %*% t(predx_cf$`ba1:ci`) +
      drawsi[,"gamma[7]"] %*% t(predx_cf$`ba2:ci`) +
      drawsi[,"gamma[8]"] %*% t(predx_cf$`ba0:ci:i`) + 
      drawsi[,"gamma[9]"] %*% t(predx_cf$`ba1:ci:i`) +
      drawsi[,"gamma[10]"] %*% t(predx_cf$`ba2:ci:i`) 
  }else{
    cf_drawsi <- 
      drawsi[,"gamma[1]"] %*% t(predx_cf$`(Intercept)`) +
      drawsi[,"gamma[2]"] %*% t(predx_cf$ba1) +
      drawsi[,"gamma[3]"] %*% t(predx_cf$ci) +
      drawsi[,"gamma[4]"] %*% t(predx_cf$i) + 
      drawsi[,"gamma[5]"] %*% t(predx_cf$`ba1:ci`) +
      drawsi[,"gamma[6]"] %*% t(predx_cf$`ba0:ci:i`) +
      drawsi[,"gamma[7]"] %*% t(predx_cf$`ba1:ci:i`)
  }
cfs[[i]] <- 
  data.frame(predx_cf, 
             mean = apply(cf_drawsi,2, mean),
             t(as.data.frame(apply(cf_drawsi, 2,
                                    quantile, probs = c(0.025,0.05,0.5,0.95,0.975)))))
names(cfs)[i] <- names(ys)[i]
}

save(cfs, file = "small_data/env_var_predy_cf.rda")
save(univar_beta_summary, univar_diagnostics, file = "small_data/env_var_model_summaries.rda")

# For inspecting parameter spreads across zero
# tidyr::as_tibble(param_summary[[12]])

ylabs = c("Variance of depth (cm)", "Mean depth (cm)","Variance of velocity (m/s)", 
          "Mean velocity (m/s)", "CPOM (g)", "CPOM with wood (g)", "FPOM (g)",
          "Algae (g)", "Macrophyte (g)", "Minimum Phi","Maximum Phi",
          "Median Phi")
cols <- c( "#CD534CFF","#EFC000FF","#868686FF","#0073C2FF")

par(mar = c(2,4,1,1))
par(mfrow = c(2,4))
for(i in 1:length(cfs)){
  cfi <- cfs[[i]]
  if(grepl("_g",names(cfs)[i])){
    cfi <- cfi[cfi$ba2 == 0,]
  }
  combos <- list(control_hii = which(cfi$ci == 0 & cfi$i == max(cfi$i)),
                 impact_hii = which(cfi$ci == 1 & cfi$i == max(cfi$i)),
                 control_lowi = which(cfi$ci == 0 & cfi$i == min(cfi$i)),
                 impact_lowi = which(cfi$ci == 1 & cfi$i == min(cfi$i)))
  x_adjs <- c(-0.03, -0.01, 0.01, 0.03)
  ylabi <- ylabs[i]
  names(cfi)[match(c("mean","X2.5.","X97.5."),names(cfi))] <- c("pred","lo","hi")
  if(y_transform_pars$stat[i] != "cpomw_g")
    cfi[c("pred","hi","lo")] <- unscale(cfi[c("pred","hi","lo")],
                                        x = samples[,y_transform_pars$stat[i]],
                                        scaled = y_transform_pars$scaled[i],
                                        log = y_transform_pars$log[i],
                                        log_add = y_transform_pars$log_add[i])$y_transformed
  cfi1 <- cfi[combos$control_hii,]
       plot(0:ifelse(grepl("_g",names(cfs)[i]), 1, 2) + x_adjs[1], 
             cfi1$pred, ylim = c(min(cfi$lo),max(cfi$hi)), 
            type = 'b', col = cols[1], pch = 16, axes = FALSE, xlim = c(0,2),
            xlab = "", ylab = ylabi)
  lines(0:ifelse(grepl("_g",names(cfs)[i]), 1,2) + x_adjs[2], cfi$pred[combos[[2]]],
                                col = cols[2], pch = 16, type = 'b')
  lines(0:ifelse(grepl("_g",names(cfs)[i]), 1,2) + x_adjs[3], cfi$pred[combos[[3]]],
                                col = cols[3], pch = 16, type = 'b')
  lines(0:ifelse(grepl("_g",names(cfs)[i]), 1,2) + x_adjs[4], cfi$pred[combos[[4]]],
                                col = cols[4], pch = 16, type = 'b')
  for(k in 1:ifelse(grepl("_g",names(cfs)[i]), 2,3)){
    for(j in 1:4){
      cfij <- cfi[combos[[j]],]
      lines(rep(k-1,2) + x_adjs[j], c(cfij$lo[k],cfij$hi[k]), col = cols[j])
    }
  }
  if(i > 2){
    axis(1, at = 0:2, labels = c("B","A1","A2"))
  }else{
    axis(1, at = 0:2, labels = c("","",""))
  }
  axis(2, las = 1); 
  box(bty = "l")
  title(main = paste0(LETTERS[i],"."), adj = 0)
  if(i == 2){
    legend("topright", pch = 16, col = cols,
           legend = c("Control, 32% EI", "Impact, 32% EI", "Control 4% EI","Impact 4% EI"),
           cex = 1)
  }
}

env_diffs <- list()
for(i in 1:12){
stanfit_i <- readRDS(stanfit_i, file = paste0("~/uomShare/wergStaff/ChrisW/git-data/urban_",
                                   "riffle_experiment/model_fits/fit_riffle_baci_",
                                   names(ys)[i], ".rds"))
drawsi <- as.data.frame(stanfit_i$draws(format = "df", variables = c("a_si","a_sa","gamma",  #"a_t",
                                                                     "sigma_si","sigma_sa","sigma"))) #"sigma_t",
cf_drawsi <- 
  drawsi[,"gamma[1]"] %*% t(predx_cf$`(Intercept)`) +
  drawsi[,"gamma[2]"] %*% t(predx_cf$ba1) +
  drawsi[,"gamma[3]"] %*% t(predx_cf$ba2) +
  drawsi[,"gamma[4]"] %*% t(predx_cf$ci) +
  drawsi[,"gamma[5]"] %*% t(predx_cf$i) + 
  drawsi[,"gamma[6]"] %*% t(predx_cf$`ba1:ci`) +
  drawsi[,"gamma[7]"] %*% t(predx_cf$`ba2:ci`) +
  drawsi[,"gamma[8]"] %*% t(predx_cf$`ba0:ci:i`) + 
  drawsi[,"gamma[9]"] %*% t(predx_cf$`ba1:ci:i`) +
  drawsi[,"gamma[10]"] %*% t(predx_cf$`ba2:ci:i`) 

env_diffs[[i]] <- baci_diffs(cf_drawsi)
meani = apply(env_diffs[[i]], 2, FUN = mean)
qlsi <- as.data.frame(t(apply(env_diffs[[i]], 2, FUN = quantile, 
      probs = c(0.025,0.05,0.125,0.5,0.875,0.95,0.975))))
names(qlsi) = c("lo95","lo90","lo75","median","hi75","hi90","hi95")
meani <- apply(env_diffs[[i]],2, FUN = mean)
if(i == 1){
  env_diff_summs <- data.frame(stat = names(cfs)[i],
                                  diff = row.names(qlsi), mean = meani, qlsi)
}else{
  env_diff_summs <- rbind(env_diff_summs,
                          data.frame(stat = names(cfs)[i],
                                     diff = row.names(qlsi), mean = meani, qlsi))
}
}

save(env_diff_summs, file = "small_data/env_diff_summs.rda")

xlabs = names(cfs)
stat_names <- ylabs
lo <- layout(matrix(c(1:8),2,4,byrow = TRUE),
             widths = c(14,10,10,10),heights = c(10,10))
diffs <- c("delta_baci1_hi","delta_baci2_hi","delta_baci1_low","delta_baci2_low")
nms <- c( bquote(Delta~"A1 4% EI"), bquote(Delta~"A2 4% EI"), 
          bquote(Delta~"A1 32% EI"), bquote(Delta~"A2 32% EI"))
miny <- 0.5; maxy <- 4.5
for(j in 8:12){
  summp <- env_diff_summs[env_diff_summs$stat == names(cfs)[j],]
  summp <- summp[match(diffs,summp$diff),]
  par(mar = c(4,ifelse(j %in% c(1,5), 5, 1),1,0))
  xrange <- c(min(summp$lo95),max(summp$hi95))
  plot(xrange,c(miny,maxy),type = 'n', axes = FALSE, xlab = xlabs[j],ylab = "")
  for(i in 1:4){
    lines(c(summp$lo95[i],summp$hi95[i]),rep(miny+maxy-i,2),lend = 2,
          col = ifelse(summp$hi90[i] < 0,"red",ifelse(summp$lo90[i] > 0, "blue","grey")))
    lines(c(summp$lo90[i],summp$hi90[i]),rep(miny+maxy-i,2), lwd = 2,lend = 2,
          col = ifelse(summp$hi90[i] < 0,"red",ifelse(summp$lo90[i] > 0, "blue","grey")))
    lines(c(summp$lo75[i],summp$hi75[i]),rep(miny+maxy-i,2), lwd = 4,lend = 2,
          col = ifelse(summp$hi90[i] < 0,"red",ifelse(summp$lo90[i] > 0, "blue","grey")))
    points(summp$mean[i],miny+maxy-i,pch = 21, 
           bg = ifelse(summp$hi90[i] < 0,"red",ifelse(summp$lo90[i] > 0, "blue","grey")))
  }
  axis(1); 
  if(j %in% c(1,5)){
    ylabs <- c("A1,32%EI","A2,32%EI","A1,4%EI","A2,4%EI")
  }else{
    ylabs <- rep("",length=length(nms))
  }
  axis(2, at = 1:4, 
       labels = ylabs, las = 1) 
  box(bty = 'l')
  abline(v = 0, lty = 3)
  title(main = paste0(LETTERS[j],". "), adj = 0)
}
par(mar = c(0,0,0,0))

