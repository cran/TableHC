# TableHC -- Higher Criticism Test between Two Frequency Tables

An adaptation of the Donoho-Jin-Tukey Higher-Critisim (HC) test to frequency tables. This adapatation uses a binomial allocation model for the number of occurances of each feature in two samples, each of which is associated with a frequency table. An exact binomial test on each feature yields a p-value. The HC statistic is used to combine these P-values to a global hypothesis test against the null hypothesis that the two tables are sampled from the same discrete distribution.

Use this test to check if two frequency tables are sampled from the same parent distribution over each entry in the table. This test is particularly useful in identifying non-null effects under weak and sparse alternatives, i.e., when the difference between the tables is due to few features, and the evidence each such feature provide is realtively weak. 

## Example:
```
n = 1000  #number of features
N = 10 * n  #number of observations

seq = seq(1,n)
P = 1 / seq  #sample from Zipf distribution
P = P / sum(P)
tb1 = data.frame(Feature = seq(1,n),  # sample 1
      Freq = rmultinom(n = 1, prob = P, size = 10*n))

k = 0.05*n  # change few features:
seq[sample(seq,k)] <- seq[sample(seq,k)] # change k features
Q = 1 / seq 
Q = Q / sum(Q) 

tb2 = data.frame(Feature = seq(1,n), # sample 2
Freq = rmultinom(n = 1, prob = Q, size = 10*n))

PV = tables.pval(tb1, tb2) # compute P-values 

HC.vals(PV$pv)  # combine P-values using the HC statistic

# The same result can be obtained using a single function call:
tables.HC(tb1,tb2) 
```
