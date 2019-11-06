#' Higher Criticism (HC) Test Between two Tables
#' @description
#' Compute HC stasitic directly from two one-way contingency tables.
#' \code{stbl} -- normalize using expeted p-value 
#' (stbl==True) or observed (stbl==False)
#' \code{alpha} -- lower fraction of p-values to use
#'
#' @param tb1 A one-way table with integer counts.
#' @param tb2 A one-way table with integer counts.
#' @param alpha A number between 0 and 1.
#' @param stbl A boolean.
#' @return A list containing the following fields:
#'  \code{HC} -- Higher Critcism (HC) score 
#'  \code{HC.star} -- HC score corrected to finite sample
#'  \code{p} -- p-value attaining \code{HC}
#'  \code{p.star} -- p-value attaining \code{HC.star}
#' 
#' @examples
#' text1 = "But our fish said, No ! No ! Make that cat go away ! 
#'            Tell that Cat in the Hat You do NOT want to play . 
#'            He should not be here . He should not be about . 
#'          He should not be here When your mother is out !"
#' text2 = "Now ! Now ! Have no fear . Have no fear! said the cat . 
#'            My tricks are not bad , Said the Cat in the Hat . Why , 
#'             we can have Lots of good fun, if you wish, with a game
#'              that I call UP - UP - UP with a fish !"
#'
#' tb1 = table(strsplit(tolower(text1),' '))
#' tb2 = table(strsplit(tolower(text2),' '))
#' pv = tables.pval(tb1,tb2)
#' HC.vals(pv$pv)
#'
#' @export
tables.HC = function(tb1,tb2, alpha = 0.45, stbl = TRUE) {
    pv = tables.pval(tb1,tb2)
    HC.vals(pv$pv, alpha = alpha, stbl = stbl)
}

#' Binomial p-values of Two-Tables
#' @description
#' Align tables and compute p-values of features using a 
#' binomial allocation model. Use outer mergeing and fill
#' missing values with zeros. 
#' 
#' @param tb1 A one-way table with integer counts.
#' @param tb2 A one-way table with integer counts.
#' @return table of pair of counts per feature and a p-value 
#' associated with each pair.
#'
#' @examples
#' tb1 = table(c(1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,4,4,4,5,6,6,7,7,7))
#' tb2 = table(c(1,1,1,1,1,1,1,1,1,2,3,3,3,3,3,4,4,4,5,5,5,6))
#' PV = tables.pval(tb1, tb2) # compute P-values 
#' HC.vals(PV$pv)  # combine P-values using the HC statistics
#'
#' text1 = "But our fish said, No ! No ! Make that cat go away ! 
#'            Tell that Cat in the Hat You do NOT want to play . 
#'            He should not be here . He should not be about . 
#'          He should not be here When your mother is out !"
#' text2 = "Now ! Now ! Have no fear . Have no fear! said the cat . 
#'            My tricks are not bad , Said the Cat in the Hat . Why , 
#'             we can have Lots of good fun, if you wish, with a game
#'              that I call UP - UP - UP with a fish !"
#'
#' tb1 = table(strsplit(tolower(text1),' '))
#' tb2 = table(strsplit(tolower(text2),' '))
#' pv = tables.pval(tb1,tb2)
#' HC.vals(pv$pv)
#'
#'
#' @export
tables.pval = function(tb1, tb2) {
    mrg = merge(tb1,tb2, all = TRUE, by = 1) #merge tables 
    mrg[is.na(mrg)] <- 0  # replace missing values with 0

    T1 = sum(mrg$Freq.x)
    T2 = sum(mrg$Freq.y) 

    bin.pv = function(n1, n2, MIN_CNT = 2){
    N = as.integer(n1 + n2);
    p.null = (T1 - n1) / (T1 + T2 - n1 - n2)
    if (is.na(p.null)) {
      p.null <- 0.5
    }
    p.value = NA;
    if(N >= MIN_CNT){ #at least MIN_CNT
            test = stats::binom.test(x = n1, n = N,
             p = p.null, alternative = "two.sided");
             p.value = test$p.value;
            }
    
    p.value
    }

    # compute p-value of each entry:
    mrg$pv = mapply(mrg$Freq.x, mrg$Freq.y, FUN = bin.pv)

    mrg
}

#' Higher Criticism Statistics
#' @description
#' Compute HC stasitic and p-value attaining it from a list of P-values.
#' Can be used with function \code{\link{tables.pval}} to 
#' get a list of p-values discriminating each feature
#' between the two tables. 
#' 
#' \code{stbl} -- normalize using expeted p-value 
#' (stbl==True) or observed (stbl==False)
#' \code{alpha} -- lower fraction of p-values to use
#'
#' @param pv A list of numbers betwee 0 and 1.
#' @param alpha A number between 0 and 1.
#' @param stbl A boolean.
#' @return A list containing the following fields:
#'  \code{HC} -- Higher Critcism (HC) score 
#'  \code{HC.star} -- HC score corrected to finite sample
#'  \code{p} -- p-value attaining \code{HC}
#'  \code{p.star} -- p-value attaining \code{HC.star}
#'
#' @examples
#' tb1 = table(c(1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,4,4,4,5,6,6,7,7,7))
#' tb2 = table(c(1,1,1,1,1,1,1,1,1,2,3,3,3,3,3,4,4,4,5,5,5,6))
#' PV = tables.pval(tb1, tb2) # compute P-values 
#' HC.vals(PV$pv)  # combine P-values using the HC statistics
#'
#'
#' text1 = "But our fish said, No ! No ! Make that cat go away ! 
#'            Tell that Cat in the Hat You do NOT want to play . 
#'            He should not be here . He should not be about . 
#'          He should not be here When your mother is out !"
#' text2 = "Now ! Now ! Have no fear . Have no fear! said the cat . 
#'            My tricks are not bad , Said the Cat in the Hat . Why , 
#'             we can have Lots of good fun, if you wish, with a game
#'              that I call UP - UP - UP with a fish !"
#'
#' tb1 = table(strsplit(tolower(text1),' '))
#' tb2 = table(strsplit(tolower(text2),' '))
#' pv = tables.pval(tb1,tb2)
#' HC.vals(pv$pv)
#'
#' @export
HC.vals = function(pv, alpha = 0.45, stbl = TRUE){ 
    pv = pv[!is.na(pv)]
    n = length(pv);
    if (n == 0) 
    {
      return(list(hc = NA, hc.stbl = NA,
         hc.star = NA,  hc.star.stbl = NA,
     p.star = NA, p.star.stbl = NA))
    }

    uu = seq(from = 1/n, to = 1, length = n) #((1:n) - 0.5) / n; 
                            #approximate expectation of p-values 
    srtd = sort(pv, index.return = TRUE);  #sorted p-vals
    ps = srtd$x

    if(stbl == TRUE) { # the flag 'stbl' indicates how to compute
                       # z-scores of p-values. The version 
                       # stbl==FALSE usually leads to more extreme
                       # values of HC.
      z = (uu - ps) / sqrt(uu * (1 - uu) + 1e-10) * sqrt(n) }
    else {
      z = (uu - ps) / sqrt(ps * (1 - ps) + 1e-10 ) * sqrt(n) }
    

    i.lim.up = floor(alpha * n + 0.5)
    
    i.max = which.max(z[1:i.lim.up]);

    if (ps[i.max] <= 2) {
      i.lim.down = min(which(ps > 1 / n))
      if (i.lim.up < i.lim.down)
        { i.lim.down = 1}
    }
    else
    {
      i.lim.down = i.max
    }

    i.max.star = which.max(z[i.lim.down:i.lim.up]) + i.lim.down - 1;
    
    list(HC = z[i.max], HC.star = z[i.max.star], p = ps[i.max],
     p.star = ps[i.max.star]);
}
