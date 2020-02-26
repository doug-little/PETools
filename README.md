# PETools
Tools to analyse and visualise ordinal data using Permutation Entropy

## Introduction
Permutation entropy (PE) is a statistic applicable to ordered (indexed) data sets [1]. PE is similar to the Shannon entropy [2], except symbols are encoded based on the *relative ordering* of data values, rather than the data values themselves. Symbols are therefore interpreted as representations of *dynamical states* - or more colloquially, patterns made by neighbouring points. The PE therefore is a measure of how evenly represented these dynamical states are in a data set, or the level of predictability in the dynamical state of the probabilistic machine/system that produced the data.

Permutation entropy is a popular statistical measure for complex-structured data, as it is
* A high-fidelity map of the continuum from fully stochastic (e.g. Gaussian white noise) to fully predictable behaviour.
* Excellent at spotting patterns, even in noisy data.
* Readily customisable to analyse behaviour over different sampling intervals.
* Relatively straightforward to compute versus other complexity measures.
* Based on evaluating inequalities only, and so is robust to noise and unaffacted by monotonic transforms in the data.
* An entropy, meaning it is readily interpretable.

#### Simple Example
To highlight the basic idea of this process, it is instructive to look at a simple example of how this process works. Let us examine the following data string {0, 1, 0.1, 0.5, 1.2, 0.4, 1.5, 2.0, 1.7, 1.2, 1.5, 100}. For the purpose of this demonstration, we will encode data based on the relative ordering of D = 2 data points (D is referred to in literature as the *embedding dimension*). First, we group the data string into lots of two as follows;
* {0, 1}
* {0.1, 0.5}
* {1.2, 0.4}
* {1.5, 2.0}
* {1.7, 1.2}
* {1.5, 100}

Next, we classify each group (or *substring*) in the following way. If x(0) > x(1) we assign the symbol "1" (or Heads), or if x(1) < x(0) we instead assign the symbol "2" (or Tails). This gives us the encoded symbol string {2, 2, 1, 2, 1, 2}. The PE is just the Shannon entropy of this string, i.e. -4/6(log(4/6)) - 2/6(log(2/6)) = 0.6365 Nats. Note the final value of 100: the fact it is an outlier of sorts has absolutely zero bearing on the computation of PE. This is an important intuitive hurdle that needs jumping - for PE magnitudes don't matter, just relative values. Whether the final data element was 1.500001 or 500 billion, it still would have returned the same symbol, and thus the same PE.

#### Some things to note
There are a number of caveats to address in this simplex example. Firstly, what if we were to group the data into groups of 3 or 4? The PE will generally change depending on our choice of this parameter. In other words the PE will be a function of some user-defined inputs. This is analogous to quantising continuous variables into discrete values before taking the Shannon Entropy - the SE we measure will depend on the quantisation scheme. It is often the case that this parameter D (called the *embedding dimension*) is held fixed for a specific application, however it is common practise to normalise the permutation entropy by dividing it by log(D!), to ensure that it falls on the interval [0, 1]. This enables PE's taken with different D to be better compared. **In fact, this normalisation is so common, that many authors implicitly include it in the definition of the PE**, a practice that we will follow here.

Secondly, some readers might question why we don't re-use data elements, i.e. why we don't group the data as {0, 1}, {1, 0.1}, {0.5, 1.2}... It is certainly possible to compute the PE grouping substrings this way (and it is included as a user-option in the `OrdEncode` function. In fact, this "overlapping" approach is arguably more common in the literature as on the surface it yields more symbols, and thus better statistical certainty. But there is a hidden cost - overlapping subsets introduce spurious *correlations* between symbols, which can greatly complicate further statistical analysis [3]. This issue is fairly nuanced and will be covered further in depth later on, but **it is recommended that, unless there is a compelling reason to do so, substrings should not overlap when grouping**. 

Next, readers may have noted that the case where x(0) == x(1) was not explicitly handled. It is generally assumed for continuous data that equalities of this nature are not possible (justified on the basis of measure theory), and that equalities if they do arise are due to inadvertent quantisation or rounding. Some authors will map equalities to a specific symbol (i.e. they will include an "or equals to" in one of the above conditions), however this method is inherently biased. A better (but not perfect) approach is to assign equalities randomly among eligible candidates. In this release this is implemeneted by **adding a vanishingly small amount of Gaussian white noise to each data element**.

Finally, what if the length of the original data set does not neatly fit into groups of D? What if our original data had 11 elements instead of 12? Quite simply, the excess symbols are dropped entirely and are not used in the calculation. Only the first D * floor(length(y)/D) data elements are used (not factoring in embedding delays > 1).

#### Embedding Delay
A further useful user parameter is the *embedding delay* which is defined as the number of indices separating the data elements in each subgrouping. Varying the embedding delay can be colloquially thought of as "probing" for patterns at different scales (time scales or otherwise). For example with an embedding delay (or "tau") or 2, the grouping in our above example would be {0, 0.1}, {1.2, 1.5}, {1.7, 1.5} and {1, 0.5}, {0.4, 2.0}, {1.2, 100}, giving {2, 2, 1} and {1, 2, 2} as the symbol strings. Note here the output string has been subdivided into 2, one for the odd-indexed elements and one for the even-indexed elements (generally the mod(tau) indexed elements). Semantically, this is expressed in the `OrdEncode` function as {2, 2, 1, 0, 1, 2, 2} with the zero indicating a break between sets. This is done to expedite the analysis of correlations between symbols, and does not affect the computation of the PE overall.

#### References
1. C. Bandt and B. Pompe, Phys. Rev. Lett. 88, 174102 (2002).
2. C.E. Shannon, W. Weaver, The Mathematical Theory of Communication (1949), Univ of Illinois Press.
3  D.J. Little, D. M. Kane, Phys. Rev. E 95, 052126 (2017).
