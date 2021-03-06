% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TableHC.R
\name{two.sample.HC}
\alias{two.sample.HC}
\title{Higher Criticism (HC) Test Between two Tables}
\usage{
two.sample.HC(tb1, tb2, alpha = 0.45, stbl = TRUE)
}
\arguments{
\item{tb1}{A one-way table with integer counts.}

\item{tb2}{A one-way table with integer counts.}

\item{alpha}{A number between 0 and 1.}

\item{stbl}{A boolean.}
}
\value{
A list containing the following fields:
 \code{HC} -- Higher Critcism (HC) score 
 \code{HC.star} -- HC score corrected to finite sample
 \code{p} -- p-value attaining \code{HC}
 \code{p.star} -- p-value attaining \code{HC.star}
}
\description{
Compute HC stasitic directly from two one-way contingency tables.
\code{stbl} -- normalize using expeted p-value 
(stbl==True) or observed (stbl==False)
\code{alpha} -- lower fraction of p-values to use
}
\examples{
text1 = "On the day House Democrats opened an impeachment inquiry of 
President Trump last week, Pete Buttigieg was being grilled by Iowa
voters on other subjects: how to loosen the grip of the rich on government, 
how to restore science to policymaking, how to reduce child poverty. At an
event in eastern Iowa, a woman rose to say that her four adult children
were `stuck' in life, unable to afford what she had in the 1980s when a 
$10-an-hour job paid for rent, utilities and an annual vacation."

text2 = "How can the federal government help our young people that want to
do whats right and to get to those things that their parents worked so hard
for? the voter asked. This is the conversation Mr. Buttigieg wants to have. 
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
pv = two.sample.pvals(tb1,tb2)
HC.vals(pv$pv)

}
