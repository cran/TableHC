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
#' text1 = "On the day House Democrats opened an impeachment inquiry of 
#' President Trump last week, Pete Buttigieg was being grilled by Iowa
#' voters on other subjects: how to loosen the grip of the rich on government, 
#' how to restore science to policymaking, how to reduce child poverty. At an
#' event in eastern Iowa, a woman rose to say that her four adult children
#' were `stuck' in life, unable to afford what she had in the 1980s when a 
#' $10-an-hour job paid for rent, utilities and an annual vacation."
#'
#' text2 = "How can the federal government help our young people that want to
#' do whats right and to get to those things that their parents worked so hard
#' for? the voter asked. This is the conversation Mr. Buttigieg wants to have. 
#' Boasting a huge financial war chest but struggling in the polls, Mr.
#' Buttigieg is now staking his presidential candidacy on Iowa, and particularly
#' on connecting with rural white voters who want to talk about personal
#' concerns more than impeachment. In doing so, Mr. Buttigieg is also trying to
#' show how Democrats can win back counties that flipped from Barack Obama to 
#' Donald Trump in 2016 — there are more of them in Iowa than any other state — 
#' by focusing, he said, on “the things that are going to affect folks’ lives in
#' a concrete way."
#'
#' tb1 = table(strsplit(tolower(text1),' '))
#' tb2 = table(strsplit(tolower(text2),' '))
#' pv = two.sample.pvals(tb1,tb2)
#' HC.vals(pv$pv)
#'
#' @export
two.sample.HC = function(tb1,tb2, alpha = 0.45, stbl = TRUE) {
    pv = two.sample.pvals(tb1,tb2)
    HC.vals(pv$pv, alpha = alpha, stbl = stbl)
}

#' Feature-by-feature exact binomial test between two tables
#' @description
#' Align tables and use an exact binomial test (binom.test)
#' on each feature. Alignment is done using "outer mergeing";
#' missing values are filled with zeros. 
#' 
#' @param tb1 A one-way table with integer counts.
#' @param tb2 A one-way table with integer counts.
#' @return table of pair of counts per feature and a p-value 
#' associated with each pair.
#'
#' @examples
#' tb1 = table(c(1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,4,4,4,5,6,6,7,7,7))
#' tb2 = table(c(1,1,1,1,1,1,1,1,1,2,3,3,3,3,3,4,4,4,5,5,5,6))
#' PV = two.sample.pvals(tb1, tb2) # compute P-values 
#' HC.vals(PV$pv)  # use the Higher-Criticism to combine the P-values 
#'                 # for a global test
#'
#' # Can be used to check similarity of word-frequencies in texts:
#'
#' text1 = "On the day House Democrats opened an impeachment inquiry of 
#' President Trump last week, Pete Buttigieg was being grilled by Iowa
#' voters on other subjects: how to loosen the grip of the rich on government, 
#' how to restore science to policymaking, how to reduce child poverty. At an
#' event in eastern Iowa, a woman rose to say that her four adult children
#' were `stuck' in life, unable to afford what she had in the 1980s when a 
#' $10-an-hour job paid for rent, utilities and an annual vacation."
#'
#' text2 = "How can the federal government help our young people that want to
#' do whats right and to get to those things that their parents worked so hard
#' for? the voter asked. This is the conversation Mr. Buttigieg wants to have. 
#' Boasting a huge financial war chest but struggling in the polls, Mr.
#' Buttigieg is now staking his presidential candidacy on Iowa, and particularly
#' on connecting with rural white voters who want to talk about personal
#' concerns more than impeachment. In doing so, Mr. Buttigieg is also trying to
#' show how Democrats can win back counties that flipped from Barack Obama to 
#' Donald Trump in 2016 — there are more of them in Iowa than any other state — 
#' by focusing, he said, on “the things that are going to affect folks’ lives in
#' a concrete way."
#'
#' tb1 = table(strsplit(tolower(text1),' '))
#' tb2 = table(strsplit(tolower(text2),' '))
#' pv = two.sample.pvals(tb1,tb2)
#' HC.vals(pv$pv)
#'
#'
#' @export
two.sample.pvals = function(tb1, tb2) {
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

#' Higher Criticism (HC) test 
#' @description
#' Compute the HC stasitic and the HC threshold given a list of P-values.
#' Can be used with function \code{\link{two.sample.pvals}} to 
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
#' PV = two.sample.pvals(tb1, tb2) # compute P-values 
#' HC.vals(PV$pv)  # combine P-values using the HC statistics
#'
#' # Can be used to check similarity of word-frequencies in texts:
#' text1 = "On the day House Democrats opened an impeachment inquiry of 
#' President Trump last week, Pete Buttigieg was being grilled by Iowa
#' voters on other subjects: how to loosen the grip of the rich on government, 
#' how to restore science to policymaking, how to reduce child poverty. At an
#' event in eastern Iowa, a woman rose to say that her four adult children
#' were `stuck' in life, unable to afford what she had in the 1980s when a 
#' $10-an-hour job paid for rent, utilities and an annual vacation."
#'
#' text2 = "How can the federal government help our young people that want to
#' do whats right and to get to those things that their parents worked so hard
#' for? the voter asked. This is the conversation Mr. Buttigieg wants to have. 
#' Boasting a huge financial war chest but struggling in the polls, Mr.
#' Buttigieg is now staking his presidential candidacy on Iowa, and particularly
#' on connecting with rural white voters who want to talk about personal
#' concerns more than impeachment. In doing so, Mr. Buttigieg is also trying to
#' show how Democrats can win back counties that flipped from Barack Obama to 
#' Donald Trump in 2016 — there are more of them in Iowa than any other state — 
#' by focusing, he said, on “the things that are going to affect folks’ lives in
#' a concrete way."
#'
#' tb1 = table(strsplit(tolower(text1),' '))
#' tb2 = table(strsplit(tolower(text2),' '))
#' pv = two.sample.pvals(tb1,tb2)
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
