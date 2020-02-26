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

Next, we classify each group (or *substring*) in the following way. If x(0) > x(1) we assign the symbol "0" (or Heads), or if x(1) < x(0) we instead assign the symbol "1" (or Tails). This gives us the encoded symbol string {1, 1, 0, 1, 0, 1}. The PE is just the Shannon entropy of this string, i.e. -4/6(log(4/6)) - 2/6(log(2/6)) = 0.6365 Nats.

#### Some things to note
There are a number of caveats to address in this simplex example. Firstly, what if we were to group the data into groups of 3 or 4? The PE will generally change depending on our choice of this parameter. In other words the PE will be a function of some user-defined inputs. This is analogous to quantising continuous variables into discrete values before taking the Shannon Entropy - the SE we measure will depend on the quantisation scheme. It is often the case that this parameter D (called the *embedding dimension*) is held fixed for a specific application, however it is common practise to normalise the permutation entropy by dividing it by log(D!), to ensure that it falls on the interval [0, 1]. This enables PE's taken with different D to be better compared. **In fact, this normalisation is so common, that many authors implicitly include it in the definition of the PE**, a practice that we will follow here.

Next, some readers might question why we don't re-use data elements, i.e. why we don't group the data as {0, 1}, {1, 0.1}, {0.5, 1.2}... It is certainly possible to compute the PE grouping substrings this way (and it is included as a user-option in the 'OrdEncode' function.

Secondly, readers may have noted that the case where x(0) == x(1) was not explicitly handled. It is generally assumed for continuous data that equalities of this nature are not possible (justified on the basis of measure theory), and that equalities if they do arise are due to inadvertent quantisation or rounding. Some authors will map equalities to a specific symbol (i.e. they will include an "or equals to" in one of the above conditions), however this method is inherently biased. A better (but not perfect) approach is to assign equalities randomly among eligible candidates. In this release this is implemeneted by **adding a vanishingly small amount of Gaussian white noise to each data element**.

Finally, what if the length of the original data set does not neatly fit into groups of D? What if our original data had 11 elements instead of 12?

#### References
1. C. Bandt and B. Pompe, Phys. Rev. Lett. 88, 174102 (2002).
2. C.E. Shannon, W. Weaver, The Mathematical Theory of Communication (1949), Univ of Illinois Press.
