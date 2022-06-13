require(pscl)
library(data.table)
library(dplyr)
require(MASS)
library(R.utils)

effect_simulated = 'syn'

####################################### Functions ####################################### 
plot_histogram <- function(model, path, reps){
  predicted_lambdas <- predict(model, type="response")
  predicted <- rpois(length(predicted_lambdas), predicted_lambdas)
  for (j in 1:(reps-1)){
    predicted <- c(predicted, rpois(length(predicted_lambdas), predicted_lambdas))
  }
  m = max(max(data$obs), max(predicted))
  hist_pred <- hist(predicted, breaks = 0:m)
  hist_obs <- hist(data$obs, breaks = 0:m)
  png(file=path,width=1200, height=900)
  log_hist <- rbind(log(hist_pred$counts/reps+1), log(hist_obs$counts+1))
  rownames(log_hist) <- c("Predicted", "Observed")
  barplot(log_hist,beside=TRUE,legend=TRUE)
  dev.off()
}

matrix_to_list_mi <- function(mutation_matrix){
  contexts <- rownames(mutation_matrix)
  colnames(mutation_matrix) <- c(0,1,2,3)
  contexts_muts <- c()
  for(item in contexts){
    for(mut in c('A', 'C', 'G', 'T')){
      if( mut != substr(item,2,2)){
        contexts_muts <- c(contexts_muts, paste(item, mut, sep=''))
      }
    }
  }
  number_muts <- 4**3*3
  values <- c()
  for(i in 0:(4^3-1)){
    for(k in 1:3){
      values <- c(values, as.character(i))
    }
  }
  muts <- c(rep(c('1','2','3'),16),rep(c('0','2','3'),16),rep(c('0','1','3'),16),rep(c('0','1','2'),16))
  values_new <- c()
  for(i in 1:length(muts)){
    values_new <- c(values_new, paste(values[i], muts[i], sep=''))
  }
  values <- values_new
  mut_mat_vec = rep(0, number_muts)
  for(i in 1:number_muts){
    x <- as.integer(substr(values[i], 0, nchar(values[i])-1))
    y <- as.integer(substr(values[i], nchar(values[i]), nchar(values[i])))
    mut_mat_vec[i] = mutation_matrix[x+1,y+1]
  }
  
  #list_mis <- as.list(setNames(as.character(as.integer(values)), as.numeric(mut_mat_vec)))
  return(cbind(as.integer(values), as.numeric(mut_mat_vec)))
}

merge_cbase_info <- function(data, cbase_path, bin_type, cpgi, suffix_cpgi, vars_to_remove=c("tri_tx","mui", "CpGI_status" ,"Nucleosomes","CRG36", "CpGme_Blood", "Recombination_bins", "DHS",  "mi", "ls", "lambdas", "l_s", "lambda_s")){
  data_cbase <- read.csv(paste(cbase_path,'output_data_preparation_', name,suffix_cpgi, '.txt', sep = ''), sep = '\t', header = FALSE)

  pos_vars <- c("chr", "start", "end", "gene", "exon")
  
  if(grepl( 'genes', bin_type, fixed = TRUE)){
    data_cbase <- data_cbase[ ,c( 'V1', 'V10')]
    colnames(data_cbase) <- c('gene', 'expected_s')
    data_combined <- merge(data, data_cbase, by = 'gene')
    selected_cols <- !names(data_combined) %in% vars_to_remove
    data_combined <- data_combined[, ..selected_cols]
    data_combined <- na.omit(data_combined)
    data_position <- data_combined[,c('gene', 'obs')]
    selected_cols <- !names(data_combined) %in% pos_vars
    data_combined <- data_combined[, ..selected_cols]
    
  }else if(grepl( 'exons', bin_type, fixed = TRUE)){
    data_cbase <- data_cbase[ ,c( 'V1', 'V4', 'V10')]
    colnames(data_cbase) <- c('gene', 'l_s', 'lambda_s')
    data_combined <- merge(data, data_cbase, by = 'gene')
    if(grepl( 'align', bin_type, fixed = TRUE)){
      data_sites <- fread(paste('/nfs/users/dweghorn/projects/Selection_sites/data/data_all_covariates_', effect_simulated, '_align_sites/data_all_covariates_', effect_simulated, '_align_sites_', name, '.bed.gz', sep = ''), sep = '\t')
    }else{
      data_sites <- fread(paste('/nfs/users/dweghorn/projects/Selection_sites/data/data_all_covariates_', effect_simulated, '_sites/data_all_covariates_', effect_simulated, '_sites_', name, '.bed.gz', sep = ''), sep = '\t')
    }
    if(cpgi == 0 || cpgi == 1){
      data_sites <- data_sites[data_sites$CpGI == cpgi, ]
    }
    mutation_matrix <- read.csv(paste(cbase_path,'mutation_matrix_', name,suffix_cpgi, '.txt', sep = ''), sep = '\t', header = FALSE, row.names = 1)
    list_mis <- matrix_to_list_mi(mutation_matrix)
    colnames(list_mis) <- c('tri_tx', 'mi')
    data_sites <- merge(data_sites, list_mis, by = 'tri_tx')
    
    data_exons <- aggregate(data_sites$mi, by=list(Category=data_sites$exon), FUN=sum)
    colnames(data_exons) <- c("exon", 'l_s_exon')
    data_combined <- merge(data_combined, data_exons, by = 'exon')
    data_combined$expected_s <- (data_combined$l_s_exon)*data_combined$lambda_s/(data_combined$l_s)
    data_combined <- data_combined[data_combined$expected_s != 0,]
    
    selected_cols <-  !names(data_combined) %in% c(vars_to_remove, 'l_s_exon')
    data_combined <- data_combined[, ..selected_cols]
    data_combined <- na.omit(data_combined)
    data_position <- data_combined[,c('gene', 'exon', 'obs')]
    selected_cols <- !names(data_combined) %in% pos_vars
    data_combined <- data_combined[, ..selected_cols]
    
  }else if(grepl( 'sites', bin_type, fixed = TRUE)){
    data_cbase <- data_cbase[ ,c( 'V1', 'V4', 'V10')]
    colnames(data_cbase) <- c('gene', 'l_s', 'lambda_s')
    data_combined <- merge(data, data_cbase, by = 'gene')
    
    mutation_matrix <- read.csv(paste(cbase_path,'mutation_matrix_', name,suffix_cpgi, '.txt', sep = ''), sep = '\t', header = FALSE, row.names = 1)
    list_mis <- matrix_to_list_mi(mutation_matrix)
    colnames(list_mis) <- c('tri_tx', 'mi')
    data_combined <- merge(data_combined, list_mis, by = 'tri_tx')
    data_combined$expected_s <- data_combined$lambda_s*(data_combined$mi+1e-12)/(data_combined$l_s)
    selected_cols <- !names(data_combined) %in% vars_to_remove
    data_combined <- data_combined[, ..selected_cols]
    data_combined <- na.omit(data_combined)
    data_position <- data_combined[,c('chr', 'start', 'end', 'exon', 'gene', 'obs')]
    selected_cols <- !names(data_combined) %in% pos_vars
    data_combined <- data_combined[, ..selected_cols]
  }
  
  return(list(data_combined, data_position))
}

scale_data <- function(data){
  #data[ , 1:ncol(data)] <- apply(data[ , 1:ncol(data)], 2,function(x) as.numeric(as.character(x)))
  data <- as.data.frame(apply(data, 2,function(x) as.numeric(as.character(x))))
  data$panc_exp <- log(1 + data$panc_exp)
  data$testis_germ_exp <- log(1 +data$testis_germ_exp)
  tryCatch({ data$exp_cancer <- log(1 +data$exp_cancer)}, error = function(e){print('no exp cancer')})
  tryCatch({ data$exp_tissue <-log(1 +data$exp_tissue)}, error = function(e){print('no exp tissue')})
  
  selected_cols <- !names(data) %in% c("obs", "expected_s")
  scaled_data <- as.data.frame(scale(data[, selected_cols]))
  return(scaled_data)
  
}

run_regressions <- function(data, scaled_data, dirname){
  model_pois_o_sd <- NA
  tryCatch({
    model_pois_o_sd = glm( data$obs ~ offset(log(data$expected_s)) + ., scaled_data,  family=poisson(link=log))
    print(summary(model_pois_o_sd))
    data$pred_pois_o_sd = predict(model_pois_o_sd, type = "response")
    png(file=paste(dirname, "/obs_vs_pred/poisson_o_sd_",name,".png",sep=""),width=1200, height=900)
    plot(data$obs, data$pred_pois_o_sd)
    dev.off()
    save( model_pois_o_sd, file=paste(dirname, "/Rmodels/poisson_o_sd_",name,".Rdata" , sep = ''))
    plot_histogram(model_pois_o_sd, paste(dirname, "/histograms/poisson_o_sd_",name,".png",sep="") , 30)
  }, error = function(e){
    print('Pois did not converge')
  })
  
  model_nb_o_sd <- NA
  tryCatch({
    model_nb_o_sd = glm.nb( data$obs ~ offset(log(data$expected_s)) + ., scaled_data)
    print(summary(model_nb_o_sd))
    data$pred_nb_o_sd = predict(model_nb_o_sd, type = "response")
    png(file=paste(dirname, "/obs_vs_pred/negbin_o_sd_",name,".png",sep=""),width=1200, height=900)
    plot(data$obs, data$pred_nb_o_sd)
    dev.off()
    save( model_nb_o_sd, file=paste(dirname, "/Rmodels/negbin_o_sd_",name,".Rdata" , sep = ''))
    plot_histogram(model_nb_o_sd, paste(dirname, "/histograms/negbin_o_sd_",name,".png",sep="") , 30)
  }, error = function(e){
    print('NB did not converge')
    })
  
  scaled_data <- cbind(scaled_data, data$obs, log(data$expected_s))
  colnames(scaled_data)[which(names(scaled_data) == 'data$obs')] <- 'obs'
  colnames(scaled_data)[which(names(scaled_data) == 'log(data$expected_s)')] <- 'log_expected_s'
  
  model_zip_o_sd <- NA
  tryCatch({
    model_zip_o_sd = zeroinfl(obs ~ offset(log_expected_s) + . | ., data=scaled_data,
                              dist ="poisson", link = "logit")
    print(summary(model_zip_o_sd))
    data$pred_zip_o_sd = predict(model_zip_o_sd, type = "response")
    png(file=paste(dirname, "/obs_vs_pred/zip_o_sd_",name,".png",sep=""),width=1200, height=900)
    plot(data$obs, data$pred_zip_o_sd)
    dev.off()
    save( model_zip_o_sd, file=paste(dirname, "/Rmodels/zip_o_sd_",name,".Rdata" , sep = ''))
    plot_histogram(model_zip_o_sd, paste(dirname, "/histograms/zip_o_sd_",name,".png",sep="") , 30)
  }, error = function(e){
    print('ZIP did not converge')})
  
  model_zinb_o_sd <- NA
  tryCatch({
    model_zinb_o_sd = zeroinfl( obs ~ offset(log_expected_s) + . | ., data=scaled_data,
                                dist ="negbin", link = "logit")
    print(summary(model_zinb_o_sd))
    data$pred_zinb_o_sd = predict(model_zinb_o_sd, type = "response")
    png(file=paste(dirname, "/obs_vs_pred/zinb_o_sd_",name,".png",sep=""),width=1200, height=900)
    plot(data$obs, data$pred_zinb_o_sd)
    dev.off()
    save( model_zinb_o_sd, file=paste(dirname, "/Rmodels/zinb_o_sd_",name,".Rdata" , sep = ''))
    plot_histogram(model_zinb_o_sd, paste(dirname, "/histograms/zinb_o_sd_",name,".png",sep="") , 30)
  }, error = function(e){
    print('ZINB did not converge')})
  list_of_models_and_data <- list(mget(c("model_pois_o_sd", "model_nb_o_sd", "model_zip_o_sd", "model_zinb_o_sd")), data)
  return(list_of_models_and_data)
}

create_dirs <- function(dirname){
  dir.create(dirname)
  dir.create(paste(dirname, "/Rmodels", sep = ''))
  dir.create(paste(dirname, "/obs_vs_pred", sep = ''))
  dir.create(paste(dirname, "/histograms", sep = ''))
  dir.create(paste(dirname, "/fcovariates", sep = ''))
}

select_best_model <- function(list_of_models, best_models_df, name, model_names, data){
  colnames <- c('pred_pois_o_sd', 'pred_nb_o_sd', 'pred_zip_o_sd', 'pred_zinb_o_sd')
  mse_list <- c()
  for(i in 1:4){
    if(sum(is.na(list_of_models[[i]]))!=1){
      mse_list <- c(mse_list, sqrt(mean((data$obs - data[[colnames[i]]])^2)))
    }else{
      mse_list <- c(mse_list, 1000000)
    }
  }
  best_model <- which.min(mse_list)
  # Store best model
  best_models_df[best_models_df$cancer_type == name,2] <- model_names[best_model]
  
  return(list(best_models_df, best_model))
}
####################################### Parameters ####################################### 
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

name <- toString(args[1]) 
bin_type <- toString(args[2])      #exons or genes or 100 or align_exons etc....
effect_considered <- toString(args[3])      #synonymous, nonsynonymous or all

print(name)
print(bin_type)
print(effect_considered)

if(effect_considered == 'synonymous'){
  cbase_path <- '/nfs/users/dweghorn/projects/Selection_sites/New_CBaSE_covariates/Output_syn/'
}else if(effect_considered == 'nonsynonymous'){
  cbase_path <- '/nfs/users/dweghorn/projects/Selection_sites/New_CBaSE_covariates/Output_nonsyn/'
}else{
  cbase_path <- '/nfs/users/dweghorn/projects/Selection_sites/New_CBaSE_covariates/Output/'
}

# Dataframe to store best model for each cancer_type
best_models_df <- data.frame('cancer_type'=name,
                             'best_model'='')
best_models_df$best_model <- as.character(best_models_df$best_model)
best_models_df_0 <- data.frame('cancer_type'=name,
                               'best_model'='')
best_models_df_0$best_model <- as.character(best_models_df_0$best_model)
best_models_df_1 <- data.frame('cancer_type'=name,
                             'best_model'='')
best_models_df_1$best_model <- as.character(best_models_df_1$best_model)
best_models_df_01 <- data.frame('cancer_type'=name,
                             'best_model'='')
best_models_df_01$best_model <- as.character(best_models_df_01$best_model)

model_names <- c("model_pois_o_sd", "model_nb_o_sd", "model_zip_o_sd", "model_zinb_o_sd")

if(grepl( 'align', bin_type, fixed = TRUE)){
  suffix <- '_align'
}else{
  suffix <- ''
}

####################################### Main Code for all cancer types #######################################

# Create dirs
dir.create( paste('/nfs/users/dweghorn/projects/Selection_sites/scripts/nextflow/results/', effect_simulated, '/',  sep = ''))
dir.create( paste('/nfs/users/dweghorn/projects/Selection_sites/scripts/nextflow/results/', effect_simulated, '/', effect_considered,  sep = ''))
dir.create( paste('/nfs/users/dweghorn/projects/Selection_sites/scripts/nextflow/results/', effect_simulated, '/', effect_considered,'/', name,   sep = ''))
dirname <- paste('/nfs/users/dweghorn/projects/Selection_sites/scripts/nextflow/results/', effect_simulated, '/', effect_considered,'/', name,'/glm', sep = '')
create_dirs(dirname)

# IF NO CpGI


if (!grepl( 'cpgi', bin_type, fixed = TRUE)){

# Read data covariates
data_cov <- fread(paste('/nfs/users/dweghorn/projects/Selection_sites/data/data_all_covariates_', effect_simulated, '_', bin_type, '/data_all_covariates_', effect_simulated, '_', bin_type, '_', name, '.bed.gz', sep = ''), sep = '\t')

# Merge info with CBaSE values
data_list <- merge_cbase_info(data_cov, cbase_path, bin_type, '', suffix)
data <- data_list[[1]]
data_position <- data_list[[2]]

# Scale data
scaled_data <- scale_data(data)

# Run regressions
list_of_models_and_data <- run_regressions(data, scaled_data, dirname)
list_of_models <- list_of_models_and_data[[1]]
data <- list_of_models_and_data[[2]]

# Select best regression
best_models_list <- select_best_model(list_of_models, best_models_df, name, model_names, data)
best_models_df <- best_models_list[[1]]
best_model <- best_models_list[[2]]

# Compute predictions
tryCatch({data_position$prediction <- predict(list_of_models[[best_model]], type = "response")}, error = function(e){data_position$prediction <<- rep(0, nrow(data_position)); print('Predict did not work :(')})
fwrite(data_position,paste(dirname, '/fcovariates/fcovariates_', bin_type, '.csv.gz', sep=""),compress="gzip", row.names = FALSE)

write.csv(best_models_df,paste(dirname, '/Rmodels/best_models',  bin_type,'.csv', sep=""), row.names = FALSE)

}

# IF WITH CpGI

if (grepl( 'cpgi', bin_type, fixed = TRUE)){
  
# First, do everything for _0
print('do everything for 0')
# Read data covariates
if(grepl( 'sites', bin_type, fixed = TRUE) == FALSE){
  data_cov_0 <- fread(paste('/nfs/users/dweghorn/projects/Selection_sites/data/data_all_covariates_', effect_simulated, '_', bin_type, '/data_all_covariates_', effect_simulated, '_', bin_type, '_', name, '_0.bed.gz', sep = ''), sep = '\t')
  }else{
  data_cov <- fread(paste('/nfs/users/dweghorn/projects/Selection_sites/data/data_all_covariates_', effect_simulated, '_sites/data_all_covariates_', effect_simulated, '_sites_', name, '.bed.gz', sep = ''), sep = '\t')
  data_cov_0 <- data_cov[data_cov$CpGI == 0, ]
  data_cov_1 <- data_cov[data_cov$CpGI == 1, ]
}
# Merge info with CBaSE values
suffix_cpgi <- paste(suffix, '_0', sep='')
data_list <- merge_cbase_info(data_cov_0, cbase_path, bin_type, 0, suffix_cpgi, c("tri_tx","mui", "CpGI", "CpGI_status" ,"Nucleosomes","CRG36", "CpGme_Blood", "Recombination_bins", "DHS",  "mi", "ls", "lambdas", "l_s", "lambda_s"))
data <- data_list[[1]]
data_position <- data_list[[2]]

data_list <- merge_cbase_info(data_cov_0, cbase_path, bin_type, 0, suffix_cpgi)
data_0 <- data_list[[1]]
data_position_0 <- data_list[[2]]

# Scale data
scaled_data <- scale_data(data)

# Run regressions
list_of_models_and_data <- run_regressions(data, scaled_data, dirname)
list_of_models <- list_of_models_and_data[[1]]
data <- list_of_models_and_data[[2]]

# Select best regression
best_models_list <- select_best_model(list_of_models, best_models_df_0, name, model_names,data)
best_models_df_0 <- best_models_list[[1]]
best_model <- best_models_list[[2]]

# Compute predictions
tryCatch({data_position$prediction <- predict(list_of_models[[best_model]], type = "response")}, error = function(e){data_position$prediction <<- rep(0, nrow(data_position)); print('Predict did not work :(')})
fwrite(data_position,paste(dirname, '/fcovariates/fcovariates_', bin_type , '_0.csv.gz', sep=""),compress="gzip", row.names = FALSE)

# Second, do everything for _1
print('do everything for 1')

# Read data covariates
if(grepl( 'sites', bin_type, fixed = TRUE) == FALSE){
  data_cov_1 <- fread(paste('/nfs/users/dweghorn/projects/Selection_sites/data/data_all_covariates_', effect_simulated, '_', bin_type, '/data_all_covariates_', effect_simulated, '_', bin_type, '_', name, '_1.bed.gz', sep = ''), sep = '\t')
}

# Merge info with CBaSE values
suffix_cpgi <- paste(suffix, '_1', sep='')
data_list <- merge_cbase_info(data_cov_1, cbase_path, bin_type, 1, suffix_cpgi,  c("tri_tx","mui", "CpGI", "CpGI_status" ,"Nucleosomes","CRG36", "CpGme_Blood", "Recombination_bins", "DHS",  "mi", "ls", "lambdas", "l_s", "lambda_s"))
data <- data_list[[1]]
data_position <- data_list[[2]]

data_list <- merge_cbase_info(data_cov_1, cbase_path, bin_type, 1, suffix_cpgi)
data_1 <- data_list[[1]]
data_position_1 <- data_list[[2]]

# Scale data
scaled_data <- scale_data(data)

# Run regressions
list_of_models_and_data <- run_regressions(data, scaled_data, dirname)
list_of_models <- list_of_models_and_data[[1]]
data <- list_of_models_and_data[[2]]

# Select best regression
best_models_list <- select_best_model(list_of_models, best_models_df_1, name, model_names, data)
best_models_df_1 <- best_models_list[[1]]
best_model <- best_models_list[[2]]

# Compute predictions
tryCatch({data_position$prediction <- predict(list_of_models[[best_model]], type = "response")}, error = function(e){data_position$prediction <<- rep(0, nrow(data_position)); print('Predict did not work :(')})
fwrite(data_position,paste(dirname, '/fcovariates/fcovariates_', bin_type , '_1.csv.gz', sep=""),compress="gzip", row.names = FALSE)

# Finally, do everything for _01
print('do everything for 01')

# Merge data
data <- rbind(data_0, data_1)
data_position <- rbind(data_position_0, data_position_1)

# Scale data
scaled_data <- scale_data(data)

# Run regressions
list_of_models_and_data <- run_regressions(data, scaled_data, dirname)
list_of_models <- list_of_models_and_data[[1]]
data <- list_of_models_and_data[[2]]

# Select best regression
best_models_list <- select_best_model(list_of_models, best_models_df_01, name, model_names, data)
best_models_df_01 <- best_models_list[[1]]
best_model <- best_models_list[[2]]

# Compute predictions
tryCatch({data_position$prediction <- predict(list_of_models[[best_model]], type = "response")}, error = function(e){data_position$prediction <<- rep(0, nrow(data_position)); print('Predict did not work :(')})
data_position$CpGI <- data$CpGI
fwrite(data_position,paste(dirname, '/fcovariates/fcovariates_', bin_type, '_01.csv.gz', sep=""),compress="gzip", row.names = FALSE)
write.csv(best_models_df_0,paste(dirname, '/Rmodels/best_models',  bin_type,'_1.csv', sep=""), row.names = FALSE)
write.csv(best_models_df_1,paste(dirname, '/Rmodels/best_models',  bin_type,'_0.csv', sep=""), row.names = FALSE)
write.csv(best_models_df_01,paste(dirname, '/Rmodels/best_models',  bin_type,'_01.csv', sep=""), row.names = FALSE)

}


