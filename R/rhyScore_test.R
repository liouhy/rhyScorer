
#' Run gene set enrichment tests based on rhythmicity scores
#'
#' This function enriches highly rhythmic and phase coherent gene sets within one group,
#' or highly differential rhythmic and phase coherent gene sets between two groups.
#'
#' @param data a matrix of normalized count data.
#' @param exp_des a data.frame describing the experimental design with two columns:
#'  group and time.
#' @param t2g a data.frame containing terms of gene sets and the corresponding genes
#' @param group a vector containing the names of two groups. The first one should
#'  be the reference group and the second one should be the testing group (eg.
#'  c("WT", "KO")).
#' @param bootstrap.n the number of time to perform bootstrapping.
#'
#' @return a data.frame containing the testing result of each gene set.
#' @importFrom dplyr %>%
#' @import foreach
#' @importFrom stats p.adjust
#' @export
#'
#' @examples
rhyScore_test = function(data, exp_des, t2g, group = NULL, bootstrap.n = 10000) {

  `%dopar%` <- foreach::`%dopar%`


  f_t2g = unique(t2g)
  f_t2g = dplyr::filter(f_t2g, f_t2g[,2] %in% rownames(data))

  # create a matrix filled with 0 with rownames as gene sets and colnames as genes we have
  gs_m = matrix(0, nrow = length(unique(f_t2g[,1])),
                ncol = nrow(data),
                dimnames = list(unique(f_t2g[,1]),rownames(data)))

  # replace the values with 1 if a gene is in a gene set
  gs_m[as.matrix(f_t2g)] = rep(1,nrow(f_t2g))


  if (length(group) == 1) {

    # harmonic regression for one group
    c = t(log2(data[,exp_des$group==group]+1))
    t = exp_des$time[exp_des$group==group]
    res = HarmonicRegression::harmonic.regression(c,t,normalize=T)

    ssr = res[['ssr']]
    para = res[['pars']]%>%dplyr::mutate(ex = amp*exp(phi*1i))%>%dplyr::mutate(sco = ex/sqrt(ssr))



    # calculate scores
    s_coh = abs(gs_m %*% para$sco) / gs_m %*% abs(para$sco)
    s_rhy = gs_m %*% abs(para$sco)

    # bootstrapping
    n = bootstrap.n


    # create the matrices with random scores
    n_cores = parallel::detectCores()-1
    my_cl = parallel::makeCluster(n_cores, type = 'PSOCK')
    doParallel::registerDoParallel(cl = my_cl)

    r_s_coh_m = foreach::foreach (j = 1:(n-1), .combine = 'cbind') %dopar% {
      r_s = sample(para$sco)
      r_s_coh = abs(gs_m %*% r_s / gs_m %*% abs(r_s))
      return(r_s_coh)
    }

    r_s_rhy_m = foreach::foreach (j = 1:(n-1), .combine = 'cbind') %dopar% {
      r_s = sample(para$sco)
      r_s_rhy = gs_m %*% abs(r_s)
      return(r_s_rhy)
    }

    parallel::stopCluster(cl = my_cl)

    # copy the column
    n_s_coh_m = s_coh[,rep(1,(n-1))]
    n_s_rhy_m = s_rhy[,rep(1,(n-1))]


    # the function to calculate bootstrapping p values
    calculate_p_val = function(n_m,r_m, n) {
      p_val = c()

      for (i in 1:nrow(n_m)) {
        if (mean(n_m[i,]) >= 0 ){
          p = (sum(r_m[i,] >= n_m[i,])+1) / n
          p_val = c(p_val,p)
        } else {
          p = (sum(r_m[i,] <= n_m[i,])+1) / n
          p_val = c(p_val,p)
        }
      }

      return(p_val)
    }

    p_coh = calculate_p_val(n_s_coh_m, r_s_coh_m, n)
    adjp_coh=p.adjust(p_coh, method='BH')

    p_rhy = calculate_p_val(n_s_rhy_m, r_s_rhy_m, n)
    adjp_rhy=p.adjust(p_rhy, method='BH')

    phi = (Im(log(gs_m %*% para$sco))+2*pi) %% (2*pi)

    count = rowSums(gs_m)

    gene = f_t2g%>%
      dplyr::filter(f_t2g[,2] %in% rownames(para))%>%
      dplyr::group_by(f_t2g[,1])%>%
      dplyr::summarise(gene = paste(f_t2g[,2], collapse = ','))

    result = data.frame('term' = rownames(gs_m),
                        'coh_score' = s_coh,
                        'rhy_socre' = s_rhy,
                        'phi' = phi,
                        'p_coh' = p_coh,
                        'adjp_coh' = adjp_coh,
                        'p_rhy' = p_rhy,
                        'adjp_rhy' = adjp_rhy,
                        'count' = count,
                        'gene' = gene$gene)


    return(result)

  } else {

    # harmonic regression for two groups
    grp_a = group[1]
    c = t(log2(data[,exp_des$group==grp_a]+1))
    t = exp_des$time[exp_des$group==grp_a]
    res = HarmonicRegression::harmonic.regression(c,t,normalize=T)

    ssr = res[['ssr']]
    a_para = res[['pars']]%>%dplyr::mutate(ex = amp*exp(phi*1i))%>%dplyr::mutate(sco = ex/sqrt(ssr))

    grp_b = group[2]
    c = t(log2(data[,exp_des$group==grp_b]+1))
    t = exp_des$time[exp_des$group==grp_b]
    res = HarmonicRegression::harmonic.regression(c,t,normalize=T)

    ssr = res[['ssr']]
    b_para = res[['pars']]%>%dplyr::mutate(ex = amp*exp(phi*1i))%>%dplyr::mutate(sco = ex/sqrt(ssr))



    # calculate diff scores
    diff_s_coh = (abs(gs_m %*% b_para$sco) / gs_m %*% abs(b_para$sco)) - (abs(gs_m %*% a_para$sco) / gs_m %*% abs(a_para$sco))
    diff_s_rhy = gs_m %*% abs(b_para$sco) - gs_m %*% abs(a_para$sco)


    # bootstrapping
    n = bootstrap.n


    # create the matrices with random scores
    n_cores = parallel::detectCores()-1
    my_cl = parallel::makeCluster(n_cores, type = 'PSOCK')
    doParallel::registerDoParallel(cl = my_cl)

    r_diff_s_coh_m = foreach::foreach (j = 1:(n-1), .combine = 'cbind') %dopar% {
      a_r_v = sample(a_para$sco)
      b_r_v = sample(b_para$sco)
      diff_r_s_coh = (abs(gs_m %*% b_r_v) / gs_m %*% abs(b_r_v))-(abs(gs_m %*% a_r_v) / gs_m %*% abs(a_r_v))
      return(diff_r_s_coh)
    }

    r_diff_s_rhy_m = foreach::foreach (j = 1:(n-1), .combine = 'cbind') %dopar% {
      a_r_v = sample(a_para$sco)
      b_r_v = sample(b_para$sco)
      r_s_rhy = (gs_m %*% abs(b_r_v)) - (gs_m %*% abs(a_r_v))
      return(r_s_rhy)
    }

    parallel::stopCluster(cl = my_cl)


    # copy the column
    n_diff_s_coh_m = diff_s_coh[,rep(1,(n-1))]
    n_diff_s_rhy_m = diff_s_rhy[,rep(1,(n-1))]


    # the function to calculate bootstrapping p values
    calculate_p_val = function(n_m,r_m, n) {
      p_val = c()

      for (i in 1:nrow(n_m)) {
        if (mean(n_m[i,]) >= 0 ){
          p = (sum(r_m[i,] >= n_m[i,])+1) / n
          p_val = c(p_val,p)
        } else {
          p = (sum(r_m[i,] <= n_m[i,])+1) / n
          p_val = c(p_val,p)
        }
      }

      return(p_val)
    }

    p_coh = calculate_p_val(n_diff_s_coh_m, r_diff_s_coh_m, n)
    adjp_coh=p.adjust(p_coh, method='BH')

    p_rhy = calculate_p_val(n_diff_s_rhy_m, r_diff_s_rhy_m, n)
    adjp_rhy=p.adjust(p_rhy, method='BH')

    a_phi = (Im(log(gs_m %*% a_para$sco))+2*pi) %% (2*pi)
    b_phi = (Im(log(gs_m %*% b_para$sco))+2*pi) %% (2*pi)

    count = rowSums(gs_m)

    gene = f_t2g%>%
      dplyr::filter(f_t2g[,2] %in% rownames(a_para))%>%
      dplyr::group_by(f_t2g[,1])%>%
      dplyr::summarise(gene = paste(f_t2g[,2], collapse = ','))

    result = data.frame('term' = rownames(gs_m),
                        'diff_coh_score' = diff_s_coh,
                        'diff_rhy_socre' = diff_s_rhy,
                        'a_phi'= a_phi,
                        'b_phi'= b_phi,
                        'p_coh' = p_coh,
                        'adjp_coh' = adjp_coh,
                        'p_rhy' = p_rhy,
                        'adjp_rhy' = adjp_rhy,
                        'count' = count,
                        'gene' = gene$gene)

    colnames(result)[colnames(result) %in% c('a_phi','b_phi')] = c(paste(grp_a,'_phi',sep=''),paste(grp_b,'_phi',sep=''))

    return(result)

  }

}
