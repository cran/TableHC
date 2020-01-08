library(testthat)

test_that("p-value from two texts", {
text1 = "On the day House Democrats opened an impeachment inquiry of 
 President Trump last week, Pete Buttigieg was being grilled by Iowa
 voters on other subjects: how to loosen the grip of the rich on government, 
 how to restore science to policymaking, how to reduce child poverty. At an
 event in eastern Iowa, a woman rose to say that her four adult children
 were stuck in life, unable to afford what she had in the 1980s when a 
 $10-an-hour job paid for rent, utilities and an annual vacation."

text2 = "How can the federal government help our young people that want to
 do what’s right and to get to those things that their parents worked so hard
 for?” the voter asked. This is the conversation Mr. Buttigieg wants to have. 
 Boasting a huge financial war chest but struggling in the polls, Mr.
 Buttigieg is now staking his presidential candidacy on Iowa, and particularly
 on connecting with rural white voters who want to talk about personal
 concerns more than impeachment. In doing so, Mr. Buttigieg is also trying to
 show how Democrats can win back counties that flipped from Barack Obama to 
 Donald Trump in 2016 — there are more of them in Iowa than any other state — 
 by focusing, he said, on “the things that are going to affect folks’ lives in
 a concrete way."

tb1 = table(strsplit(tolower(text1),' '))
tb2 = table(strsplit(tolower(text2),' '))

expect_equal(abs(HC.vals(two.sample.pvals(tb1,tb2)$pv)$p-0.555)<0.01,TRUE)
})

