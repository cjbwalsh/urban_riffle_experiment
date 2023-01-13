source("misc_functions.R")
# load and check cmdstanr package
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
# Set random seed to ensure reproducible results
rand_seed = 9430570 #1 #2

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

mod <- cmdstan_model("nbinom_me_rand_sa_si_t_fixedmatrix_noeps_baci.stan", pedantic = TRUE) 

ni <- 3000; nt <- 4; nb <- 2000; nc <- 4
stanfit_i <- mod$sample(data = sdata,
                        seed = rand_seed, chains = nc,
                        parallel_chains = nc, iter_warmup = nb,
                        iter_sampling = ni - nb, refresh = 30)
# 6000 iterations 3.6 h
# #  save csv files rather than the model object to use less RAM
stanfit_i$save_output_files(
  dir = "~/uomShare/wergStaff/ChrisW/git-data/urban_riffle_experiment/model_fits/",
  basename = "fit_riffle_baci_c", timestamp = FALSE, random = FALSE)
stanfit_i$sampler_diagnostics()
saveRDS(stanfit_i, file = "~/uomShare/wergStaff/ChrisW/git-data/urban_riffle_experiment/model_fits/fit_riffle_baci_c.rds")
stanfit_i$diagnostic_summary()
# EBFMI 0.273, 0.289, 0.290, 0.309, zero divergences, zero max treedepth reached.
# # The above three steps required < 500 Mb RAM
summ <- stanfit_i$summary() # This took ~2h and needed >40 Gb RAM
min(summ$ess_bulk,na.rm=TRUE) # 657 from 6500...rerun this with 5500 tonight!
min(summ$ess_tail,na.rm=TRUE) # 1333
max(summ$rhat,na.rm=TRUE)  # 1.005