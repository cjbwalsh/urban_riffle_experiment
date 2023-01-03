library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
rand_seed = 9430572

source("https://tools.thewerg.unimelb.edu.au/documents/mwstr/mwstr_functions.R")
source("https://tools.thewerg.unimelb.edu.au/data/mwbugs/bugDatabaseFunctions.R")

## Load data: ultimately from OSF
# library(osfr); library(dplyr)
# if(!"data" %in% dir()){system("mkdir data")}
# if(!"wq_data_compiled.xlsx" %in% dir("data")){
# wq_files <- osf_retrieve_node("4ywvq") %>% osf_ls_files()
# osf_download(wq_files[wq_files$name == "wq_data_compiled.xlsx",], path = "data")
# }
# compile data 
sites <- readxl::read_excel("~/uomShare/wergStaff/ChrisW/git-data/urban_riffle_experiment/data_for_model.xlsx", sheet = "sites")
samples <- readxl::read_excel("~/uomShare/wergStaff/ChrisW/git-data/urban_riffle_experiment/data_for_model.xlsx", sheet = "samples")
biota <- readxl::read_excel("~/uomShare/wergStaff/ChrisW/git-data/urban_riffle_experiment/data_for_model.xlsx", sheet = "biota")
subcs <- sqlQuery(paste0("SELECT site, reach FROM subcs WHERE reach IN ('",
                         paste(substr(sites$sitecode,1,nchar(sites$sitecode) -1),
                               collapse = "', '"), "');"), "mwstr")
load("~/uomShare/wergStaff/ChrisW/git-data/mwstr/mwstr_v12_corrections/imp_subcs.rda")
sites$site <- subcs$site[match(substr(sites$sitecode,1,nchar(sites$sitecode) -1),
                               subcs$reach)]
sites$ai <- imp_subcs$c_ai[match(sites$site, imp_subcs$site)]

samples$seg <- substr(samples$old_samplecode,nchar(samples$old_samplecode)-1,
                      nchar(samples$old_samplecode)-1)
# aggregate(samples$smpcode, by = list(site = samples$sitecode, 
#                                      t = samples$t, seg = samples$seg), FUN = length)
samples <- samples[samples$seg == "M",]
biota <- biota[biota$smpcode %in% samples$smpcode,]
biota_ct <- as.data.frame(ct(biota$smpcode, biota$shortcode, biota$count))
biota_ct <- biota_ct[match(samples$smpcode,row.names(biota_ct)),]
taxa <- readxl::read_excel("~/uomShare/wergStaff/ChrisW/git-data/urban_riffle_experiment/biota_ct.xlsx", sheet = "taxa")
# # test with subset of taxa 
# biota_ct <- biota_ct[apply(biota_ct, 2, FUN = function(x) sum(x > 0)) > 100] #9 taxa
ss_ct <- biota_ct
for(i in 1:nrow(samples)){
  ss_ct[i,] <- samples$subsample_perc[i]/100
}
for(i in which(biota$coarsepick == 1)){
  ss_ct[row.names(ss_ct) == biota$smpcode[i], biota$shortcode[i]] <- 1
}
samples$ba <- as.numeric(as.numeric(substr(samples$old_samplecode,1,1))  > 2) 
samples$ba[as.numeric(substr(samples$old_samplecode,1,1))  > 4] <- 2
# 0 = before, 1 = 1 year after , 2 = 5 year after
samples$ci <- as.numeric(sites$exp_treatment[match(samples$sitecode,sites$sitecode)] == "riffle") #0 = control, 1 = riffle
samples$baci <- samples$ba*samples$ci
samples$ai <- sites$ai[match(samples$sitecode, sites$sitecode)]
i_scaled <- scale(log10(samples$ai*100 + 0.1))
samples$i <- as.vector(i_scaled)
sites <- sites[order(sites$exp_treatment,sites$ai),]
sites$site_no <- 1:nrow(sites)
samples$site_no <- sites$site_no[match(samples$sitecode,sites$sitecode)]
samples$sample <- substr(samples$old_samplecode,1,nchar(samples$old_samplecode)-1)
sample_nos <- data.frame(sample = unique(samples$sample))
sample_nos$sample_no <- 1:nrow(sample_nos)
samples$sample_no <- sample_nos$sample_no[match(samples$sample, sample_nos$sample)]
samples$t <- as.numeric(substr(samples$old_samplecode, 1,1))
samples$bai <- samples$ba * samples$i
samples$bacii <- samples$baci * samples$i
samples$spring <- as.integer(substr(samples$old_samplecode,1,1) %in% c(1,3,5))
u <- model.matrix(~ ba + ci + baci + i + bacii + spring, data = samples)
sdata <- list(n_obs = nrow(biota_ct),
              n_taxa = ncol(biota_ct),
              n_site = nrow(sites),
              n_sample = nrow(sample_nos),
              n_pred = ncol(u),
              n_t = max(samples$t),
              site_no = samples$site_no,
              samp_no = samples$sample_no,
              t = samples$t,
              u = u,
              c = as.matrix(biota_ct),
              s = as.matrix(ss_ct)
)
mod <- cmdstan_model("nbinom_me_randomsite_fixedmatrix_baci.stan", pedantic = TRUE) 
ni <- 4600; nt <- 4; nb <- 1000; nc <- 4
stanfit_i <- mod$sample(data = sdata,
                        seed = rand_seed, chains = nc,
                        parallel_chains = nc, iter_warmup = nb,
                        iter_sampling = ni - nb, refresh = 200)
# 4400 iters 2.5 h
#  save csv files rather than the model object to use less RAM
stanfit_i$save_output_files(
  dir = "~/uomShare/wergStaff/ChrisW/git-data/urban_riffle_experiment/model_fits/",
  basename = "fit_nbinom_ba012_cii_2000iter_ll", timestamp = FALSE, random = FALSE)
stanfit_i$sampler_diagnostics()
saveRDS(stanfit_i, file = "~/uomShare/wergStaff/ChrisW/git-data/urban_riffle_experiment/model_fits/fit_nbinom_ba012_cii_2000iter_ll.rds")
stanfit_i$diagnostic_summary()

# The above three steps required < 500 Mb RAM
summ <- stanfit_i$summary() # This took > 2h and needed ~80 Gb Ram
min(summ$ess_bulk,na.rm=TRUE) # 475 from 4600 in 3.8 h
min(summ$ess_tail,na.rm=TRUE) # 748
# all diagnostic statistics fine.