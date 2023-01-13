### This trial demonstrated that a model with no prior hyper-distribution for 
### sigma_sa, sigma_si, and sigma_t produced weaker predictions despite sampling
### more quickly with higher BFMI.  This got me thinking about the desirability 
### of setting a sigma_sa etc. value for each taxon individually (trialling it 
### on the simpler 46-site data as I type.)

source("misc_functions.R")
# load and check cmdstanr package
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
# Set random seed to ensure reproducible results
rand_seed = 9430568 #9 #70 #1 #2

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
# Reduce biota table to match reduced samples table
biota <- biota[biota$smpcode %in% samples$smpcode,]
# Convert long-form biota table into a wide taxon-by-sample table of counts 
# in subsamples
biota_ct <- as.data.frame(ct(biota$smpcode, biota$shortcode, biota$count))
biota_ct <- biota_ct[match(samples$smpcode,row.names(biota_ct)),]
# Create a table of the same dimensions as biota_ct, with the subsample proportion
# for each observation (for coarsepick specimens, subsample ppn = 1)
ss_ct <- biota_ct
for(i in 1:nrow(samples)){
  ss_ct[i,] <- samples$subsample_perc[i]/100
}
for(i in which(biota$coarsepick == 1)){
  ss_ct[row.names(ss_ct) == biota$smpcode[i], biota$shortcode[i]] <- 1
}

# Assemble random predictors (site_no, sample_no, t)
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
samples$spring <- as.integer(substr(samples$old_samplecode,1,1) %in% c(1,3,5))
u <- model.matrix(~ ba + ci + ba:ci + i + ba:ci:i + spring, data = samples)

# Data list for Stan
sdata <- list(n_obs = nrow(biota_ct),
              n_taxa = ncol(biota_ct),
              n_site = nrow(sites),
              n_sample = nrow(sample_nos),
              n_pred = ncol(u),
              n_t = max(samples$t),
              site_no = samples$site_no,
              samp_no = samples$sample_no,
              t_no = samples$t,
              u = u,
              c = as.matrix(biota_ct),
              s = as.matrix(ss_ct)
)

mod <- cmdstan_model("nbinom_me_no_sigma_rand_sa_si_t_fixedmatrix_baci.stan", pedantic = TRUE) 

ni <- 5000; nt <- 4; nb <- 1500; nc <- 4
stanfit_i <- mod$sample(data = sdata,
                        seed = rand_seed, chains = nc,
                        parallel_chains = nc, iter_warmup = nb,
                        iter_sampling = ni - nb, refresh = 100)
# 6000 iterations 3.6 h
# #  save csv files rather than the model object to use less RAM
stanfit_i$save_output_files(
  dir = "~/uomShare/wergStaff/ChrisW/git-data/urban_riffle_experiment/model_fits/",
  basename = "fit_riffle_baci_mt_no_sigma", timestamp = FALSE, random = FALSE)
stanfit_i$sampler_diagnostics()
saveRDS(stanfit_i, file = "~/uomShare/wergStaff/ChrisW/git-data/urban_riffle_experiment/model_fits/fit_riffle_baci_mt_no_sigma.rds")
stanfit_i$diagnostic_summary()
# EBFMI 0.64, 0.65, 0.64, 0.65, zero divergences, zero max treedepth reached.
# # The above three steps required < 500 Mb RAM
summ <- stanfit_i$summary() # This took ~2h and needed >40 Gb RAM
min(summ$ess_bulk,na.rm=TRUE) # 414
min(summ$ess_tail,na.rm=TRUE) # 774
max(summ$rhat,na.rm=TRUE)  # 1.011

####

mod_draws <- as.data.frame(stanfit_i$draws(format = "df", variables = c("a_si","a_sa","a_t","gamma","phi")))
save(mod_draws, 
     file = "~/uomShare/wergStaff/ChrisW/git-data/urban_riffle_experiment/model_fits/fit_no_sigma_riffle_baci_draws.rda")

####

predx <- unique(data.frame(site_no = samples$site_no,
                           sample_no = samples$sample_no,
                           t = samples$t,u))
predy_draws <- list()
system.time({
  for(i in 1:sdata$n_taxa){
    drawsi <- mod_draws[grep(paste0(",",i,"]"), names(mod_draws))]
    predy_draws[[i]] <-
      drawsi[grep("a_si",names(drawsi))][match(predx$site_no, 1:sdata$n_site)] +
      drawsi[grep("a_sa",names(drawsi))][match(predx$sample_no, 1:sdata$n_sample)] +
      drawsi[grep("a_t",names(drawsi))][match(predx$t, 1:sdata$n_t)] +
      drawsi[,paste0("gamma[1,",i,"]")] %*% t(predx$X.Intercept.) +
      drawsi[,paste0("gamma[2,",i,"]")] %*% t(predx$ba1) +
      drawsi[,paste0("gamma[3,",i,"]")] %*% t(predx$ba2) +
      drawsi[,paste0("gamma[4,",i,"]")] %*% t(predx$ci) +
      drawsi[,paste0("gamma[5,",i,"]")] %*% t(predx$i) + 
      drawsi[,paste0("gamma[6,",i,"]")] %*% t(predx$spring) +
      drawsi[,paste0("gamma[7,",i,"]")] %*% t(predx$ba1.ci) +
      drawsi[,paste0("gamma[8,",i,"]")] %*% t(predx$ba2.ci) +
      drawsi[,paste0("gamma[9,",i,"]")] %*% t(predx$ba0.ci.i) + 
      drawsi[,paste0("gamma[10,",i,"]")] %*% t(predx$ba1.ci.i) +
      drawsi[,paste0("gamma[11,",i,"]")] %*% t(predx$ba2.ci.i) 
  }
})  # 1 min
names(predy_draws) <- colnames(sdata$c)
save(predy_draws, file = "~/uomShare/wergStaff/ChrisW/git-data/urban_riffle_experiment/model_fits/fit_no_sigma_riffle_baci_predy_draws.rda")

####

predx_cf <- expand.grid(intercept = 1, ba = factor(c(0,1,2)), ci = c(0,1), 
                        i = seq(min(u[,5]),max(u[,5]),length=10), spring = 0)
predx_cf <- data.frame(model.matrix(~ ba + ci + ba:ci + i + ba:ci:i + spring, 
                                    data = predx_cf))
predx_cf$ai <- 10^(predx_cf$i * attr(i_scaled, 'scaled:scale') + 
                     attr(i_scaled, 'scaled:center')) - 0.1
predy_cf <- list()
system.time({
  for(i in 1:sdata$n_taxa){
    drawsi <- mod_draws[grep(paste0(",",i,"]"), names(mod_draws))]
    predy_cf[[i]] <-
      drawsi[,paste0("gamma[1,",i,"]")] %*% t(predx_cf$X.Intercept.) +
      drawsi[,paste0("gamma[2,",i,"]")] %*% t(predx_cf$ba1) +
      drawsi[,paste0("gamma[3,",i,"]")] %*% t(predx_cf$ba2) +
      drawsi[,paste0("gamma[4,",i,"]")] %*% t(predx_cf$ci) +
      drawsi[,paste0("gamma[5,",i,"]")] %*% t(predx_cf$i) + 
      drawsi[,paste0("gamma[6,",i,"]")] %*% t(predx_cf$spring) +
      drawsi[,paste0("gamma[7,",i,"]")] %*% t(predx_cf$ba1.ci) +
      drawsi[,paste0("gamma[8,",i,"]")] %*% t(predx_cf$ba2.ci) +
      drawsi[,paste0("gamma[9,",i,"]")] %*% t(predx_cf$ba0.ci.i) + 
      drawsi[,paste0("gamma[10,",i,"]")] %*% t(predx_cf$ba1.ci.i) +
      drawsi[,paste0("gamma[11,",i,"]")] %*% t(predx_cf$ba2.ci.i) 
  }
})  # 5 s
names(predy_cf) <- names(biota_ct)
save(predy_cf, 
     file = "~/uomShare/wergStaff/ChrisW/git-data/urban_riffle_experiment/model_fits/fit_no_sigma_riffle_baci_predy_cf_draws.rda")

####

phi_draws <- mod_draws[,grep("phi",names(mod_draws))]
p_cf_mat <- list()
n_cf_mat <- list()
# system.time({
for(t in 1:ncol(biota_ct)){
  p_cf_mat[[t]] <- predy_cf[[t]]
  n_cf_mat[[t]] <- exp(predy_cf[[t]])
  for(i in 1:ncol(p_cf_mat[[t]]))
    # Probability of occurrence in 1 sample unit
    p_cf_mat[[t]][,i] <- 1 - pnbinom(0,mu = exp(predy_cf[[t]][,i]),
                                     size = phi_draws[t][,1])
  # Probability of occurrence in 4 sample units
  p_cf_mat[[t]] <- 1 - (1 - p_cf_mat[[t]])^4
  # back-transformed abundance per 4 samples
}
# }) # ~15 s
names(p_cf_mat) <- colnames(biota_ct)
tot_rich_cf_mat <- Reduce('+', p_cf_mat)
# see https://stackoverflow.com/questions/11641701/sum-a-list-of-matrices
ept_rich_cf_mat <- Reduce('+', p_cf_mat[substr(names(p_cf_mat),1,2) %in% c("QE","QP","QT")])
do_rich_cf_mat <- Reduce('+', p_cf_mat[substr(names(p_cf_mat),1,2) %in% c("QD","LO")])

tot_n_cf_mat <- Reduce('+', n_cf_mat)
ept_n_cf_mat <- Reduce('+', n_cf_mat[substr(names(p_cf_mat),1,2) %in% c("QE","QP","QT")])
do_n_cf_mat <- Reduce('+', n_cf_mat[substr(names(p_cf_mat),1,2) %in% c("QD","LO")])

library(vegan)

# Shannon's H' requires proportional abundances
relab_cf_mat <- predy_cf

for(i in 1:length(p_cf_mat)){
  relab_cf_mat[[i]] <- exp(p_cf_mat[[i]])
}
tot_abund_cf_mat <- Reduce('+', relab_cf_mat)
for(i in 1:length(p_cf_mat)){
  relab_cf_mat[[i]] <- relab_cf_mat[[i]] / tot_abund_cf_mat
}


H_mat <- relab_cf_mat[[1]]
for(i in 1:nrow(H_mat)){
  H_mat[i,] <- diversity(sapply(relab_cf_mat,function(x, n) x[n,],i))
}

# Pielou's Evenness
evenness_cf_mat <- H_mat/log(tot_rich_cf_mat)

# Reduce set to family-presence absence for SIGNAL calculation
biotic_indices<- 
  read.csv("~/uomShare/wergStaff/ChrisW/git-data/urban_riffle_experiment/biotic_indices.csv")
signal <- biotic_indices[!is.na(biotic_indices$SIGNALWoV2003),
                         c("shortcode","taxon","SIGNALWoV2003")]
signal <- signal[signal$shortcode %in%
                   c("LO",substr(names(biota_ct),1,4)),]
signal <- signal[order(signal$shortcode),]
signal_mat <- list()
for(i in 1:nrow(signal)){
  ss <- ifelse(signal$shortcode[i] == "LO",
               which(substr(names(p_cf_mat),1,2) == "LO"),
               which(substr(names(p_cf_mat),1,4) %in% signal$shortcode[i]))
  signal_mat[[i]] <-
    Reduce(function(x) 1-(1-x)^length(ss), p_cf_mat[ss])
  signal_mat[[i]] <- matrix(unlist(Map(function(x) rbinom(1,1,x),
                                       signal_mat[[i]])),ncol = ncol(p_cf_mat[[1]]))
  signal_mat[[i]][signal_mat[[i]] > 0] <-
    signal_mat[[i]][signal_mat[[i]] > 0]* signal$SIGNALWoV2003[i]
}
signal_sum_cf_mat <- Reduce('+', signal_mat)
signal_len_cf_mat <- Reduce('+', lapply(signal_mat, function(x)
  replace(x, x > 0, 1)))
signal_cf_mat <- signal_sum_cf_mat/signal_len_cf_mat

assemb_stats <- list(tot_rich = 
                       data.frame(mean =apply(tot_rich_cf_mat,2, 
                                              FUN = mean),
                                  t(apply(tot_rich_cf_mat,2, 
                                          FUN = quantile, 
                                          probs = c(0.025, 0.975)))),
                     evenness = 
                       data.frame(mean =apply(evenness_cf_mat,2, 
                                              FUN = mean),
                                  t(apply(evenness_cf_mat,2, 
                                          FUN = quantile, 
                                          probs = c(0.025, 0.975)))),
                     ept_rich = 
                       data.frame(mean =apply(ept_rich_cf_mat,2, 
                                              FUN = mean),
                                  t(apply(ept_rich_cf_mat,2, 
                                          FUN = quantile, 
                                          probs = c(0.025, 0.975)))),
                     do_rich = 
                       data.frame(mean = apply(do_rich_cf_mat,2, 
                                               FUN = mean),
                                  t(apply(do_rich_cf_mat,2, 
                                          FUN = quantile, 
                                          probs = c(0.025, 0.975)))),
                     signal = 
                       data.frame(mean =apply(signal_cf_mat,2, 
                                              FUN = mean, na.rm = TRUE),
                                  t(apply(signal_cf_mat,2, 
                                          FUN = quantile, 
                                          probs = c(0.025, 0.975), 
                                          na.rm = TRUE))),
                     tot_n = 
                       data.frame(mean =apply(tot_n_cf_mat,2, 
                                              FUN = mean),
                                  t(apply(tot_n_cf_mat,2, 
                                          FUN = quantile, 
                                          probs = c(0.025, 0.975)))),
                     ept_n = 
                       data.frame(mean =apply(ept_n_cf_mat,2, 
                                              FUN = mean),
                                  t(apply(ept_n_cf_mat,2, 
                                          FUN = quantile, 
                                          probs = c(0.025, 0.975)))),
                     do_n = 
                       data.frame(mean =apply(do_n_cf_mat,2, 
                                              FUN = mean),
                                  t(apply(do_n_cf_mat,2, 
                                          FUN = quantile, 
                                          probs = c(0.025, 0.975))))
)
save(assemb_stats, file = "small_data/no_sigma_assemb_stats.rda")

####

prevalence <- apply(biota_ct, 2, FUN = function(x){sum(x > 0)})
prevalence_by_sample <- data.frame(shortcode = colnames(biota_ct), prev = NA)
for(i in 1:ncol(biota_ct)){
  prevalence_by_sample$prev[i] <- 
    sum(aggregate(biota_ct[,i], by = list(sample_no = samples$sample_no), 
                  FUN = function(x) as.numeric(sum(x) > 0))$x)
}
n_obs_samps <- 1000

####

load("~/uomShare/wergStaff/ChrisW/git-data/urban_riffle_experiment/model_fits/riffle_baci_obs_totcount_estimates.rda") #list object called obs

# prepare output matrices for estimates of T (obs), and correlations and slopes betwen O and P
cors <- matrix(nrow = n_obs_samps, ncol = sdata$n_taxa)
oe_slopes <- matrix(nrow = n_obs_samps, ncol = sdata$n_taxa)
system.time({
  for(i in 1:sdata$n_taxa){
    y <- log(exp(apply(predy_draws[[colnames(sdata$c)[i]]],2,FUN = median)) + 1)
    for(j in 1:n_obs_samps){
      x <- aggregate(log(obs[[i]][j,] + 1),
                     by = list(site_no = samples$site_no, 
                               sample_no = samples$sample_no,
                               t = samples$t), FUN = mean)$x
      cors[j,i] <- cor(x,y)
      oe_slopes[j,i] <- coefficients(lm(y ~ x))[2]
      # me[j,i] <- exp(x) - exp(y)
      # 
      #                               meanError = c(((exp(mean(obsAMsByTC$totAbun) + 
      #                                                   mean(predAMsByTC$totAbun_050 - obsAMsByTC$totAbun)) - 1)/
      #                                            (exp(mean(obsAMsByTC$totAbun)) -1) - 1)*100,
      #                                           ((exp(mean(obsAMsByTC$eptAbun) + 
      #                                                   mean(predAMsByTC$eptAbun_050 - obsAMsByTC$eptAbun)) - 1)/
      #                                            (exp(mean(obsAMsByTC$eptAbun)) -1) - 1)*100,
      #                                           mean(obsAMsByTC$totRich - predAMsByTC$totRich_050),
      #                                           mean(obsAMsByTC$eptRich - predAMsByTC$eptRich_050)),
      # 
    }
    if(i %% 5 == 0)
      cat(i,"\n")
  }
}) # 3.8 min
save(cors, oe_slopes, file = "~/uomShare/wergStaff/ChrisW/git-data/urban_riffle_experiment/model_fits/no_sigma_riffle_baci_oe_stats.rda")

####

#### Duplicated in methods/results...Not necessary, if managed differently?
mod_draws <- get(load("~/uomShare/wergStaff/ChrisW/git-data/urban_riffle_experiment/model_fits/fit_no_sigma_riffle_baci_draws.rda"))
params <- c("b_ba1", "b_ba2", "b_ci","b_i","b_spring","b_ba1ci","b_ba2ci","b_ba0cii","b_ba1cii", "b_ba2cii")
df_empty <- data.frame(shortcode = NA, mean = NA, lo95 = NA, lo90 = NA, median = NA, hi90 = NA, hi95 = NA)[0,]
param_summs <- list()
for(i in 1:length(params)){
  drawsi <- mod_draws[substr(names(mod_draws),1,ifelse(i < 9, 7, 8)) == paste0("gamma[", i + 1) & 
                        grepl("gamma", names(mod_draws))]
  qlsi <- as.data.frame(t(apply(drawsi,2, FUN = quantile, 
                                probs = c(0.025,0.05,0.125,0.5,0.875,0.95,0.975))))
  names(qlsi) = c("lo95","lo90","lo75","median","hi75","hi90","hi95")
  meani <- apply(drawsi,2, FUN = mean)
  param_summs[[i]] <- data.frame(shortcode = names(biota_ct), 
                                 mean = meani,
                                 qlsi)
}
names(param_summs) <- params

ba1ci_pos_taxa <- param_summs$b_ba1ci$shortcode[param_summs$b_ba1ci$lo90 > 0] #0
ba1ci_neg_taxa <- param_summs$b_ba1ci$shortcode[param_summs$b_ba1ci$hi90 < 0] #0
ba2ci_pos_taxa <- param_summs$b_ba2ci$shortcode[param_summs$b_ba2ci$lo90 > 0] #0
ba2ci_neg_taxa <- param_summs$b_ba2ci$shortcode[param_summs$b_ba2ci$hi90 < 0] #0
ba0cii_pos_taxa <- param_summs$b_ba0cii$shortcode[param_summs$b_ba0cii$lo90 > 0] #19
ba0cii_neg_taxa <- param_summs$b_ba0cii$shortcode[param_summs$b_ba0cii$hi90 < 0] #7
ba1cii_pos_taxa <- param_summs$b_ba1cii$shortcode[param_summs$b_ba1cii$lo90 > 0] #19
ba1cii_neg_taxa <- param_summs$b_ba1cii$shortcode[param_summs$b_ba1cii$hi90 < 0] #0
ba2cii_pos_taxa <- param_summs$b_ba2cii$shortcode[param_summs$b_ba2cii$lo90 > 0] #4
ba2cii_neg_taxa <- param_summs$b_ba2cii$shortcode[param_summs$b_ba2cii$hi90 < 0] #0
taxa_to_plot <- unique(c(ba0cii_pos_taxa,ba0cii_neg_taxa,
                         ba1ci_pos_taxa,ba1ci_neg_taxa,ba1cii_pos_taxa,ba1cii_neg_taxa,
                         ba2ci_pos_taxa,ba2ci_neg_taxa,ba2cii_pos_taxa,ba2cii_neg_taxa,
                         param_summs$b_ba$shortcode[param_summs$b_ba$lo90 > 0],
                         param_summs$b_ba$shortcode[param_summs$b_ba$hi90 < 0],
                         param_summs$b_ci$shortcode[param_summs$b_ci$lo90 > 0],
                         param_summs$b_ci$shortcode[param_summs$b_ci$hi90 < 0],
                         param_summs$b_i$shortcode[param_summs$b_i$lo90 > 0],
                         param_summs$b_i$shortcode[param_summs$b_i$hi90 < 0],
                         param_summs$b_spring$shortcode[param_summs$b_spring$lo90 > 0],
                         param_summs$b_spring$shortcode[param_summs$b_spring$hi90 < 0]))
higher_taxa_in_plot <- higher_taxa[higher_taxa$higher_taxon %in% 
                                     unique(taxa$higher_taxon[taxa$taxoncode %in% taxa_to_plot]),]

# R results for inclusion in figure caption:
ht_string <- paste(paste0(higher_taxa_in_plot$taxon, " (", higher_taxa_in_plot$higher_taxon, ")"), collapse = "; ")
# Annelida (A); Pelecypoda (B); Coleoptera (C); Diptera (D); Ephemeroptera (E); Turbellaria (F); Gastropoda (G); Lepidoptera (L); Acarina (M); Odonata (O); Plecoptera (P); Trichoptera (T); Crustacea (Z); Cnidaria (I); Nematoda (J)
ntoplot <- length(taxa_to_plot) #59
####

bgs <- rep(NA, ncol(cors))
bgs[names(biota_ct) %in% ba1cii_pos_taxa] <- "red"
  bgs[names(biota_ct) %in% ba0cii_neg_taxa] <- "green"
    cols <- bgs; cols[is.na(cols)] <- "black"
      layout(matrix(c(1,2,3,3),2,2,byrow=TRUE), widths = c(12,12), heights = c(12,1))
      par(mar = c(2,4,1,1))
      xpos <- jitter(prevalence_by_sample$prev)
      plot(xpos, apply(cors,2,mean),
           xlim = c(0,53), ylim = c(-0.3,1), pch = 21, bg = bgs,
           axes = FALSE, ylab = "R (O:P)", xlab = "")
      abline(h = 1, lty = 3)
      for(i in 1:ncol(cors)){
        lines(rep(xpos[i],2), quantile(cors[,i],probs = c(0.025,0.975)),
              col = scales::alpha(cols[i],0.5))
      }
      axis(1); axis(2, las = 1); box(bty = "l")
      title("A.", adj = 0)
      plot(xpos, apply(oe_slopes,2,mean),
           ylim = c(-0.1,1.1),pch = 21, bg = bgs,
           axes = FALSE, ylab = "Regression slope (O ~ P)", xlab = "")
      abline(h = 1, lty = 3)
      for(i in 1:ncol(cors)){
        lines(rep(xpos[i],2), quantile(oe_slopes[,i],probs = c(0.025,0.975)),
              col = scales::alpha(cols[i],0.5))
      }
      axis(1); axis(2, las = 1); box(bty = "l")
      title("B.", adj = 0)
      par(mar = c(0,0,0,0))
      plot.new()
      title(xlab = "Prevalence (number of samples with abundance > 0)", line = -1.5)

#### 
      
      phi_draws <- mod_draws[,grep("phi",names(mod_draws))]
      p_mat <- list()
      # system.time({
      for(t in 1:sdata$n_taxa){
        p_mat[[t]] <- predy_draws[[colnames(biota_ct)[t]]]
        for(i in 1:ncol(p_mat[[t]]))
          p_mat[[t]][,i] <- 1 - pnbinom(0,mu = exp(predy_draws[[colnames(biota_ct)[t]]][,i]),size = phi_draws[t][,1])
      }
      # }) # ~15 s
      names(p_mat) <- colnames(biota_ct)
      p_mean <- matrix(nrow = sdata$n_sample, ncol= sdata$n_taxa)
      for(t in 1:sdata$n_taxa){
        p_mean[,t] <- apply(p_mat[[t]], 2, FUN = mean)
      }
      
      # Function for calculating probability of occurrence in n sample units, given a vector of probabilities, x, of length n
      p_per_n <- function(x) {  # x a vector of probabilities, length n
        #  aggregate(x, by = list(pair = rep(1:(length(x)/n), each = n)), 
        #            FUN = function(y) {
        1 - prod(1 - x)  # })
      }
      
      p_mean_per_sample <- matrix(ncol = sdata$n_sample, nrow = sdata$n_taxa)
      for(t in 1:sdata$n_taxa){
        for(j in 1:sdata$n_sample){
          p_mean_per_sample[t,j] <- p_per_n(rep(p_mean[j,t],sum(samples$sample_no == j)))
        }}
      
      # Probability of presence in observed data
      obs_per_sample <- matrix(ncol = sdata$n_sample, nrow = sdata$n_taxa)
      for(t in 1:sdata$n_taxa){
        obs_per_sample[t,] <- aggregate(biota_ct[,t], 
                                        list(sample_no = samples$sample_no), 
                                        FUN = function(y) {sum(y > 0) > 0})$x
      }
      
      library("PresenceAbsence")
      pa_stats <- data.frame(shortcode = colnames(biota_ct),
                             taxon = taxa$taxon[match(colnames(biota_ct),taxa$shortcode)], 
                             prevalence = prevalence_by_sample$prev[match(prevalence_by_sample$shortcode, colnames(biota_ct))],
                             PCC = NA,
                             sensitivity = NA, specificity = NA, Kappa = NA, AUC = NA)
      for(t in 1:sdata$n_taxa){
        pa_data <- data.frame(taxon = colnames(biota_ct)[t],
                              obs_per_sample[t,],p_mean_per_sample[t,])
        pa_stats[t,-(1:3)] <- presence.absence.accuracy(pa_data)[c("PCC","sensitivity","specificity","Kappa","AUC")] 
      }
      par(mar = c(4,4,1,1),mfrow = c(1,1))
      plot(pa_stats$prevalence, pa_stats$AUC, ylim = c(0.5,1),xlab = "", ylab = "AUC")
      title(xlab = "Prevalence (number of samples with abundance > 0)",)
      
####  

      tot_rich_mat <- Reduce('+', p_mat)
      # see https://stackoverflow.com/questions/11641701/sum-a-list-of-matrices
      ept_rich_mat <- Reduce('+', p_mat[substr(names(p_mat),1,2) %in% c("QE","QP","QT")])
      
      tot_rich_pred_mean <- apply(tot_rich_mat,2, FUN = mean)
      tot_rich_pred_cls <- apply(tot_rich_mat,2, FUN = quantile, probs = c(0.025, 0.975))
      
      ept_rich_pred_mean <- apply(ept_rich_mat,2, FUN = mean)
      ept_rich_pred_cls <- apply(ept_rich_mat,2, FUN = quantile, probs = c(0.025, 0.975))
      
      tot_rich_obs <- apply(obs_per_sample, 2, FUN = function(y){sum(y > 0)})
      ept_rich_obs <- apply(obs_per_sample[substr(colnames(biota_ct),1,2) %in% c("QE","QP","QT"),], 2, FUN = function(y){sum(y > 0)})
      
      par(mfrow = c(1,2))
      tot_lm <- lm(tot_rich_pred_mean ~ tot_rich_obs)
      plot(tot_rich_obs, tot_rich_pred_mean, xlim = c(0,40), ylim = c(0,40),axes = FALSE,
           xlab  = "Observed total richness", ylab = "Predicted total richness")
      for(i in 1:sdata$n_sample){
        lines(rep(tot_rich_obs[i],2), tot_rich_pred_cls[,i])
      }
      axis(1); axis(2, las = 1); box(bty = 'l')
      title("A.", adj = 0)
      abline(0,1)
      abline(tot_lm, lty = 3)
      legend("topleft", legend = NA, box.lty = 0,
             title = paste0("R-sq = ",round(summary(tot_lm)$r.squared,2)))
      
      ept_lm <- lm(ept_rich_pred_mean ~ ept_rich_obs)
      plot(ept_rich_obs, ept_rich_pred_mean, xlim = c(0,10), ylim = c(0,10), axes = FALSE,
           xlab  = "Observed EPT richness", ylab = "Predicted EPT richness")
      for(i in 1:sdata$n_sample){
        lines(rep(ept_rich_obs[i],2), ept_rich_pred_cls[,i])
      }
      axis(1); axis(2, las = 1); box(bty = 'l')
      title("B.", adj = 0)
      abline(0,1)
      abline(ept_lm, lty = 3)
      legend("topleft", legend = NA, , box.lty = 0,
             title = paste0("R-sq = ",round(summary(ept_lm)$r.squared,2)))
