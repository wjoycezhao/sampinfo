#' Format data for the thoughtListing project
#'
#' Add feature matrix
#' (decision congruence, semantic congruence and cluster resmapling matrix)
#' to the thought listing data.
#' @param cluster_num number of thought cluster for one option. If 3 then option number is 2 * 3 + 1 = 7.
#' @param data_all dataframe witth qID, sID, tNo, and cID
#' @param dist_all distance matrix
#'
#' @export
#'
getFormatData = function(data_all, dist_all, cluster_num = NULL) {
  ## thought listing data
  # data_all = read.csv(paste0(data_folder, exp_name, '_data_c', cluster_num, '.csv'))
  option_num = 2 * cluster_num + 1
  ## add cluster resampling feature for predicting the next trial
  data_all = addFeatureMatrix(diag(option_num),
                              data_all, option_num, 'sameC', 'cID')
  ## add decision congruence for predicting the next trial
  data_all = addFeatureMatrix(getSameA(cluster_num),
                              data_all, option_num, 'sameA', 'cID')
  ## add semantic distance for predicting the next trial
  # dist_all = read.csv(paste0(data_folder, exp_name, '_dist_c', cluster_num, '.csv'),
  #                     row.names = 1)
  data_temp = lapply(unique(data_all$qID), function (qq)
    addFeatureMatrix(dist_all[dist_all$qID == qq, -1],
                     data_all[data_all$qID == qq, ],
                     option_num, 'dist', 'cID'))
  data_all = do.call(rbind.data.frame, data_temp)

  ## move the features to the associated current thoughts
  col_temp = unlist(lapply(c('sameC', 'sameA', 'dist'), function(x)
    paste0(x, '_', 1:option_num)))
  data_all[-1, col_temp] = data_all[-nrow(data_all), col_temp]

  ## feature value should be 0 for the first thoughts
  data_all[data_all$tNo == 1, col_temp] = 0

  ## check thought index
  for (i in 2:nrow(data_all)) {
    if ((data_all$sID[i] == data_all$sID[i - 1] &
         data_all$tNo[i] != data_all$tNo[i - 1] + 1) |
        (data_all$sID[i] != data_all$sID[i - 1] &
         data_all$tNo[i] != 1)) {
      stop('wrong thought index')
    }
  }

  return(data_all)
}

#' Format data for the thoughtListing project
#'
#' Add feature matrix
#' (decision congruence, semantic congruence and cluster resmapling matrix)
#' to the thought listing data.
#' @param cluster_num number of thought cluster for one option. If 3 then option number is 2 * 3 + 1 = 7.
#' @param exp_name name of the experiment (used in the file name)
#' @param data_folder the folder in which .csv is saved
#'
#' @export
#'
getData = function(cluster_num = NULL,
                       exp_name,
                       data_folder = '') {
  ## thought listing data
  data_all = read.csv(paste0(data_folder, exp_name, '_data_c', cluster_num, '.csv'))
  option_num = 2 * cluster_num + 1
  ## add cluster resampling feature for predicting the next trial
  data_all = addFeatureMatrix(diag(option_num),
                              data_all, option_num, 'sameC', 'cID')
  ## add decision congruence for predicting the next trial
  data_all = addFeatureMatrix(getSameA(cluster_num),
                              data_all, option_num, 'sameA', 'cID')
  ## add semantic distance for predicting the next trial
  dist_all = read.csv(paste0(data_folder, exp_name, '_dist_c', cluster_num, '.csv'),
                      row.names = 1)
  data_temp = lapply(unique(data_all$qID), function (qq)
    addFeatureMatrix(dist_all[dist_all$qID == qq, -1],
                     data_all[data_all$qID == qq, ],
                     option_num, 'dist', 'cID'))
  data_all = do.call(rbind.data.frame, data_temp)

  ## move the features to the associated current thoughts
  col_temp = unlist(lapply(c('sameC', 'sameA', 'dist'), function(x)
    paste0(x, '_', 1:option_num)))
  data_all[-1, col_temp] = data_all[-nrow(data_all), col_temp]

  ## feature value should be 0 for the first thoughts
  data_all[data_all$tNo == 1, col_temp] = 0

  ## check thought index
  for (i in 2:nrow(data_all)) {
    if ((data_all$sID[i] == data_all$sID[i - 1] &
         data_all$tNo[i] != data_all$tNo[i - 1] + 1) |
        (data_all$sID[i] != data_all$sID[i - 1] &
         data_all$tNo[i] != 1)) {
      stop('wrong thought index')
    }
  }

  ## get score distributions for all the clusters
  cluster_rating_m = data_all %>%
    dplyr::count(qID, cID, rating) %>%
    tidyr::complete(qID, cID, rating,
                    fill = list(n = 0))  %>%
    dplyr::group_by(qID, cID) %>%
    dplyr::mutate(n = n / sum(n)) %>%
    tidyr::spread(rating, n)
  cluster_rating_m[is.na(cluster_rating_m)] = 0

  ## get base rates
  base_rates = data_all %>% dplyr::group_by(cID) %>% dplyr::tally()
  base_rates$p = base_rates$n / sum(base_rates$n)
  return(
    list(
      format_data = data_all,
      cluster_rating_m = cluster_rating_m,
      base_rates = base_rates
    )
  )
}

#' Generate decision congruence matrix
#'
#' Generate the decision congruence matrix.
#'
#' @inheritParams getData
#'
#' @return A matrix. 1 for clusters supporting the same decision; 0 otherwise.
#' @export
#'
#' @examples
getSameA = function(cluster_num = NULL) {
  sameA_m = diag(cluster_num * 2 + 1)
  temp = matrix(1, ncol = cluster_num, nrow = cluster_num)
  sameA_m[1:cluster_num, 1:cluster_num] = temp
  sameA_m[((1 * cluster_num) + 2):((2 * cluster_num) + 1), (1 * cluster_num +
                                                              2):(2 * cluster_num + 1)] = temp
  return(sameA_m)
}

# addFeature = function(feature_mv, data_df, option_num, feature_name, id=NULL){
#   if (is.vector(feature_mv)){
#     feature_df = addFeatureMatrix(feature_mv, data_df, option_num, feature_name, id=NULL)
#   }else{
#     feature_df = addFeatureVector(feature_mv, data_df, option_num, feature_name, id=NULL)
#   }
#   return(feature_df)
# }

#' Add feature matrix
#'
#' @param feature_mv The features to be added
#' @param data_df The raw data.frame
#' @param option_num Number of options, should be the same as the column number of the  feaature matrix
#' @param feature_name Name of feature to appear in formatted data
#' @param id Indicating which rows to extract from the feature matrix for each row of the raw data frame
#'
#' @export
#'

addFeatureMatrix = function(feature_mv,
                            data_df,
                            option_num,
                            feature_name,
                            id = NULL) {
  if (nrow(feature_mv) != option_num) {
    stop('number of columns of feature matrix should be the same as the option_num')
  }
  feature_mv1 = feature_mv[data_df[, id], ]
  colnames(feature_mv1) = paste0(feature_name, '_', 1:option_num)
  return(cbind.data.frame(data_df, feature_mv1))
}

# addFeatureVector = function(feature, data_df, option_num, feature_name, id=NULL){
#   if(length(feature_mv) != option_num){
#     stop('length of feature vector should be the same as the option_num')}
#   feature_mv1 = feature_mv[data_df[,id]]
#   colnames(feature_mv1) = feature_name
#   return(cbind.data.frame(data_df, feature_mv1))
# }
