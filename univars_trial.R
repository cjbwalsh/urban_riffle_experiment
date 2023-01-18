library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
# Set random seed to ensure reproducible results
rand_seed = 9430568

## Load data: ultimately from OSF repository
# library(osfr); library(dplyr)
# if(!"data" %in% dir()){system("mkdir data")}
# if(!"wq_data_compiled.xlsx" %in% dir("data")){
# wq_files <- osf_retrieve_node("4ywvq") %>% osf_ls_files()
# osf_download(wq_files[wq_files$name == "wq_data_compiled.xlsx",], path = "data")
# }
# data_for_model.xlsx compiled in urban_riffle_exp_data_compilation.R
dat_file <- 
  "~/uomShare/wergStaff/ChrisW/git-data/urban_riffle_experiment/data_for_model.xlsx"
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
ys <- list(depth_m_var = log(samples$depth_m_var + 1.75e-04),
           # min(samples$depth_m_var[!is.na(samples$depth_m_var) & samples$depth_m_var > 0]) # 0.0007442
           depth_m_mean = log(samples$depth_m_mean),
           vel_m_s_var = log(samples$vel_m_s_var + 0.0007442),
           # min(samples$vel_m_s_var[!is.na(samples$vel_m_s_var) & samples$vel_m_s_var > 0]) # 0.0007442
           vel_m_s_mean = log(samples$vel_m_s_mean), 
           cpom_g = log(samples$cpom_g + 0.001),
           # min(samples$cpom_g[!is.na(samples$cpom_g) & samples$cpom_g > 0]) # 0.001
           cpomw_g = log(samples$cpom_g + samples$wood_g + 0.001),
           fpom_g = log(samples$fpom_g  + 0.0094382),
           # min(samples$fpom_g[!is.na(samples$fpom_g) & samples$fpom_g > 0]) # 0.0094382
           algae_g = log(samples$algae_g + 0.001), 
           # min(samples$algae_g[!is.na(samples$algae_g) & samples$algae_g > 0]) # 0.001
           macrophyte_g = log(samples$macrophyte_g + 5e-04),
           # min(samples$macrophyte_g[!is.na(samples$macrophyte_g) & samples$macrophyte_g > 0]) # 5e-04
           min_phi = samples$min_phi,
           max_phi = samples$max_phi,
           med_phi = samples$med_phi,
           range_phi = samples$range_phi)

mod <- cmdstan_model("normal_rand_sa_si_t_fixedmatrix_baci.stan", pedantic = TRUE) 

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

for(i in 1:13){
# remove NAs, and check if data covers trips 5 and 6
  yi <- names(ys)[i]
samplesi <- samples[!is.na(samples[ifelse(yi == "cpomw_g","cpom_g",yi)]),]
if(length(unique(samplesi$sitecode)) != nrow(sites)) stop()
# no om data for EUM...need to check this...for now...
if(i %in% 5:8){
samplesi <- samplesi[-grep("EUM",samples$sitecode),]
sitesi <- sites[-grep("EUM",sites$sitecode),]
sitesi$site_no <- 1:nrow(sitesi)
samplesi$site_no <- sitesi$site_no[match(samplesi$sitecode, sitesi$sitecode)]
}
# recode ba for cases where there are no data for ba 2 (trips 5 and 6)
samplesi$ba <- factor(samplesi$ba)
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
              n_t = max(samplesi$t),
              site_no = samplesi$site_no,
              samp_no = samplesi$sample_no,
              t_no = samplesi$t,
              u = u,
              y = ys[[i]][match(samplesi$smpcode,samples$smpcode)]
)

ni <- ifelse(i %in% c(3,5,6), 9000, 3500); nt <- 4; 
nb <- ifelse(i %in% c(3,5,6), 3000, 1500); nc <- 4
stanfit_i <- mod$sample(data = sdata,
                        seed = rand_seed, chains = nc,
                        parallel_chains = nc, iter_warmup = nb,
                        adapt_delta = ifelse(i %in% c(5,8,9,11),0.95,
                                             ifelse(i == 1, 0.85,0.8)),  
                        max_treedepth = ifelse(i %in% c(2,6),15,10),
                        iter_sampling = ni - nb, refresh = 100)
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
min(summ$ess_bulk,na.rm=TRUE) # 657 from 6500...rerun this with 5500 tonight!
min(summ$ess_tail,na.rm=TRUE) # 1333
max(summ$rhat,na.rm=TRUE)  # 1.005

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

drawsi <- as.data.frame(stanfit_i$draws(format = "df", variables = c("a_si","a_sa","a_t","gamma")))
predx <- unique(cbind(samplesi[c("sample","site_no","sample_no","t")],u))
temp <- aggregate(sdata$y, by = list(sample = samplesi$sample),FUN = mean)
predx$x <- temp$x[match(predx$sample, temp$sample)]
  if("ba2" %in% names(as.data.frame(u))){
    predy_draws <-
      drawsi[grep("a_si",names(drawsi))][match(predx$site_no, 1:sdata$n_site)] +
      drawsi[grep("a_sa",names(drawsi))][match(predx$sample_no, 1:sdata$n_sample)] +
      drawsi[grep("a_t",names(drawsi))][match(predx$t, 1:sdata$n_t)] +
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
      drawsi[grep("a_t",names(drawsi))][match(predx$t, 1:sdata$n_t)] +
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

# For inspecting parameter spreads across zero
# tidyr::as_tibble(param_summary[[12]])

ylabs = c("Variance of depth (m)", "Mean depth (m)","Variance of velocity (m/s)", 
          "Mean velocity (m/s)", "CPOM (g)", "FPOM (g)",
          "Algae (g)", "Macrophyte (g)", "Minimum Phi","Maximum Phi",
          "Median Phi","Range Phi")
cols <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")

par(mar = c(2,4,1,1))
par(mfrow = c(2,4))
for(i in 1:length(cfs)){
  cfi <- cfs[[i]]
  combos <- list(control_lowi = which(cfi$ci == 0 & cfi$i == min(cfi$i)),
                 control_hii = which(cfi$ci == 0 & cfi$i == max(cfi$i)),
                 impact_lowi = which(cfi$ci == 1 & cfi$i == min(cfi$i)),
                 impact_hii = which(cfi$ci == 1 & cfi$i == max(cfi$i)))
  x_adjs <- c(-0.03, -0.01, 0.01, 0.03)
  ylabi <- ylabs[i]
  names(cfi)[match(c("mean","X2.5.","X97.5."),names(cfi))] <- c("pred","lo","hi")
  with(cfi[combos$control_lowi,], 
       plot((0:ifelse(length(combos$control_lowi) == 2, 1,2)) + x_adjs[1], pred, 
             ylim = c(min(cfi$lo),max(cfi$hi)), 
            type = 'b', col = cols[1], pch = 16, axes = FALSE, xlim = c(0,2),
            xlab = "", ylab = ylabi))
  with(cfi[combos[[2]],], lines((0:ifelse(length(combos$control_lowi) == 2, 1,2))+ x_adjs[2], pred,
                                col = cols[2], pch = 16, type = 'b'))
  with(cfi[combos[[3]],], lines((0:ifelse(length(combos$control_lowi) == 2, 1,2)) + x_adjs[3], pred,
                                col = cols[3], pch = 16, type = 'b'))
  with(cfi[combos[[4]],], lines((0:ifelse(length(combos$control_lowi) == 2, 1,2)) + x_adjs[4], pred,
                                col = cols[4], pch = 16, type = 'b'))
  for(k in 1:ifelse(length(combos$control_lowi) == 2, 2,3)){
    for(j in 1:4){
      cfij <- cfi[combos[[j]],]
      lines(rep(k-1,2) + x_adjs[j], c(cfij$lo[k],cfij$hi[k]), col = cols[j])
    }
  }
  if(i > 2){
    axis(1, at = 0:2, labels = c("Before","2 y after","5 y after"))
  }else{
    axis(1, at = 0:2, labels = c("","",""))
  }
  axis(2, las = 1); 
  box(bty = "l")
  title(main = paste0(LETTERS[i],"."), adj = 0)
  if(i == 2){
    legend("topright", pch = 16, col = cols,
         legend = c("Ctl, low EI", "Ctl high EI","Imp, low EI", "Imp high EI"),
         cex = 0.9)
}
