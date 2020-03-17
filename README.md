# PETools
Tools to analyse and visualise ordinal data using Permutation Entropy

## 1. Introduction
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

Next, we classify each subgroup (or *substring*) in the following way. If x(0) > x(1) we assign the symbol "1" (or Heads), or if x(1) < x(0) we instead assign the symbol "2" (or Tails). This gives us the encoded symbol string {2, 2, 1, 2, 1, 2}. The PE is just the Shannon entropy of this string, i.e. -4/6(log(4/6)) - 2/6(log(2/6)) = 0.6365 Nats. Note the final value of 100: the fact it is an outlier of sorts has absolutely zero bearing on the computation of PE. This is an important intuitive hurdle that needs jumping - for PE magnitudes don't matter, just relative values. Whether the final data element was 1.500001 or 500 billion, it still would have returned the same symbol, and thus the same PE.

#### Some things to note
There are a number of caveats to address in this simple example. Firstly, what if we were to group the data into groups of 3 or 4? The PE will generally change depending on our choice of this parameter. In other words the PE will be a function of some user-defined inputs. This is analogous to quantising continuous variables into discrete values before taking the Shannon Entropy - the SE we measure will depend on the quantisation scheme. It is often the case that this parameter D (called the *embedding dimension*) is held fixed for a specific application, however it is common practise to normalise the permutation entropy by dividing it by log(D!), to ensure that it falls on the interval [0, 1]. This enables PE's taken with different D to be better compared and removes the dependence on the logarithm base. **In fact, this normalisation is so common, that many authors implicitly include it in the definition of the PE**, a practice that we will follow here.

Secondly, some readers might question why we don't re-use data elements, i.e. why we don't group the data as {0, 1}, {1, 0.1}, {0.5, 1.2}... It is certainly possible to compute the PE grouping substrings this way (and it is included as a user-option in the `OrdEncode` function. In fact, this "overlapping" approach is arguably more common in the literature as on the surface it yields more symbols, and thus better statistical certainty. But there is a hidden cost - overlapping subsets introduce spurious *correlations* between symbols, which can greatly complicate further statistical analysis [3]. This issue is fairly nuanced and will be covered further in depth later on, but **it is recommended that, unless there is a compelling reason to do so, substrings should not overlap when grouping**. 

Next, readers may have noted that the case where x(0) == x(1) was not explicitly handled. It is generally assumed for continuous data that equalities of this nature are not possible (justified on the basis of measure theory), and that equalities if they do arise are due to inadvertent quantisation or rounding. Some authors will map equalities to a specific symbol (i.e. they will include an "or equals to" in one of the above conditions), however this method is inherently biased. A better (but not perfect) approach is to assign equalities randomly among eligible candidates. In the `OrdEncode` function this is implemeneted by **adding a vanishingly small amount of Gaussian white noise to each data element**.

Finally, what if the length of the original data set does not neatly fit into groups of D? What if our original data had 11 elements instead of 12? Quite simply, the excess symbols are dropped entirely and are not used in the calculation. Only the first D * floor(length(y)/D) data elements are used (not factoring in embedding delays > 1, more on embedding delays below).

## 2. Definitions
### Permutation Entropy
In this section, we will take the ideas introduced in the previous section and place them on a firmer mathematical footing. The (normalised) permutation entropy (or just PE for short) is computed as

(2.1) *H*<sub>est</sub>(*D*,*&tau;*) = log(*D*!)<sup>-1</sup> &Sigma;<sub>s</sub> (*p*<sub>*i*</sub>/*N*)log(*p*<sub>*i*</sub>/*N*)

where *D* is the embedding dimension, *N* denotes the total number of symbols, *p*<sub>*i*</sub> denote counts of the *i*th symbol, and the summation is over the complete set of symbols, *i* &isin; *s*. The symbols are assigned to each data subgroup based on which permutation of {1:*D*!} matches the condition

*x*<sub>1</sub> < *x*<sub>2</sub> < ... < *x*<sub>*D*!</sub>,

for the given group. Note here the subscripts refer to indices *within the subgroup*, not the indices of the original full data set.

The "est" subscript in equation 2.1 denotes this to be an estimated PE, which will be discussed in the next section below. Note that as a consequence of the limit of *x*log(*x*) being 0 as *x* &rarr; 0, empty (zero) counts can be safely removed from the summation. It is noted explicitly that the PE is a funtion of the embedding dimension as well as *D*, as it is common in PE analysis to compute *H* as a function of *&\tau;*, though it should also be remarked that there is an overall dependence on how the data elements are mapped to the groups before encoding and things like e.g. whether there is overlap in sub-groupings are reliant on context, so it is imperative that these details are noted for the purposes of reproducibility.

### Interpretation of Permutation Entropy
It is possible to define both the PE of the input ordinal data, and the PE of the machine/system that generated the ordinal data. The former is a *deterministic* interpretation, while the latter must be regarded as a *probabilistic* interpretation. Unfortunately, these two interpretations are often conflated in scientific literature.

The connection between the deterministic/probabilistic interpretations of PE is nothing more than the connection between the statsitics of a data set, and a sampled subset of the original set. In general the sample mean will *not* equal the true mean, the sample variance will not equal the true variance, etc. The sampled subset can be used to make inferences about the larger data set, but not with infinite precision. So it is with permutation entropy: we can use the PE of the ordinal data to make inferences of the machine/system that generated it, but not with infinite precision.

The permuation entropy of the machine/system (hereafter referred to as the *true* permutation entropy) is given by

(2.2) *H*(*D*,*&tau;*) = log(*D*!)<sup>-1</sup> &Sigma;<sub>s</sub> *P*<sub>*i*</sub>log(*P*<sub>*i*</sub>),

which depend on *P*<sub>*i*</sub>, the *probability* of obtaining the symbol corresponding to the label *i*. Note that this is *not* equivalent to the formula above, since *p*<sub>*i*</sub>/*N* will not equal *P*<sub>*i*</sub> due to sampling error. Equation 2.1 is therefore an *estimator* for the "true" permutation entropy, equation 2.2. It is this nuance that many users fail to navigate, leading to errors in analysis. As PE is often used in inference testing e.g. to detect changes in the probabilistic output of a machine/system, knowing how to calculate the precision of a PE estimate is critical, especially in instances where the amount of available data is limited.

### Embedding Dimension
The size of the subsets the ordinal data is divided into before the encoding step of the calculation. Incidentally, the name *permutation entropy* derives from the fact that each ordinal pattern (or *dynamical state*, to borrow the language used earlier) corresponds to a permutation of the string {1:D}, thus there are *D*! possible dynamic states to be summed over.

Because the number of dynamic states increases rapidly with increasing *D*, the recommended range of this parameter is quite narrow. In the `OrdEncode` function, this will return an error message unless 2 <= *D* <= 8. Not only do higher *D* natively require more computation time due to the increased number of inequality operations that need to be performed, they also need a longer input data string to establish an equivalent level of statistical precision (which blows out computation times even further). It is generally rare to see meaningful patterns detected in data that cannot be seen using D <= 5. It is therefore strongly recommended to remain within these bounds.

So given all this, what is the most appropriate *D* to use for a given data set (in addition to the above considerations pertaining to computation time)? The benefit of using higher *D* is the ability to spot more complex patterns in the input data. Consider an input string with four data points - using *D* = 2 gives two patterns, each with two possible outcomes, so 4 possible outcomes in total. Using *D* = 4 on the other hand yields 24 possible outcomes. The general rule of thumb is to make *D* large enough to spot any patterns of interest, but no higher since the costs rapidly outweight the benefits of having a higher *D*. Of course, knowing which *D* will spot the patterns of interest will initially require some exploration. A good habit in PE analysis is to derive what insights you can from low *D* values before progressively working your way to higher *D*. Analyses at lower *D* can also help provide context for outputs at higher *D*, which can sometimes be tricky to analyse in the absence of this context.

### Embedding Delay
A further useful user parameter is the *embedding delay* which is defined as the index separation between data elements in each subgrouping. Unless specified, it can usually be accepted that &tau; = 1 as the default case. Varying the embedding delay can be colloquially thought of as "probing" for patterns at different scales (time scales or otherwise). For example with an embedding delay (or "tau") or 2, the grouping in our above example would be {0, 0.1}, {1.2, 1.5}, {1.7, 1.5} and {1, 0.5}, {0.4, 2.0}, {1.2, 100}, giving {2, 2, 1} and {1, 2, 2} as the symbol strings. Note here the output string has been subdivided into 2, one for the odd-indexed elements and one for the even-indexed elements (generally the mod(tau) indexed elements). Semantically, this is expressed in the `OrdEncode` function as {2, 2, 1, 0, 1, 2, 2} with the zero indicating a break between sets. This is done to expedite the analysis of correlations between symbols, and does not affect the computation of the PE overall.

Aside from probing patterns at different scales, using a larger embedding delay can be beneficial for data sets exhibiting strong short-range correlations. Oversampled data (where the sample rate exceeds the bandwidth of the system producing the output) commonly exhibits these correlations, and so it is useful in such circumstances to set &tau; to a higher value, so that the effective sampling rate matches the system bandwidth; in effect, removing these strong short-range correlations, which can lead to correlations in the resultant symbol string that can mask other significant behaviours.

## 3. PE statistics for iid-distributed symbols
The data analysis procedure outlined thus far involved taking an ordinal data set, breaking it up into subsets and assigning each subset a label (which we will hereafter refer to as a *symbol*) based on the positional order (or *index*) of each element within the subset going from lowest to highest or vice versa. As a side point: the actual symbol labels are completely interchangable. In fact, a property of permutation entropy (or any statistic derived from the categorisation of subsets in this fashion) is invariance under interchanges of labels. It is useful to include some kind of iterator in the labels for notational convenience, labels of the form &pi;<sub>*i*</sub> are popular in the literature. 

In statistics parlance, "iid" is a standard abbreviation that stands for *independent and identically-distributed*. *Independent* here means that the probabilities (labelled *P*<sub>*i*</sub>) of observing different symbols remain invariant over all subsets, and do not depend on previous observations. It **does not** mean that the original ordinal data set contains iid elements. Taking symbols to iid is a strong approximation that greatly simplifies statistical analysis, however a note of caution for the user - real data sets are *rarely* iid. Nonetheless, looking at the iid case is a useful starting point, and in some circumstances, a meaningful limit case.

Note too that allowing data subsets to contain overlapping element renders any iid assumptions invalid for obvious reasons. If two data subsets share elements, they cannot be regarded as independent, hence why it is generally recommended not to have these subsets overlap.

### PE of White noise
White noise (GWN) is defined as a set of uncorrelated samples drawn from a distribution. Because PE is distribution-agnostic when it is applied to stochatsic data series, the distribution the noise samples are drawn from does not matter. Further, we can consider the case where the mean is zero and the standard deviation is unity without loss of generality. It is straightforward to show that the "true" PE as defined by equation 2.2 is exactly 1 via symmetry arguments. Because the samples are uncorrelated, any order of points must be equally likely to any other order, thus *P*<sub>*i*</sub> must all be equal to 1/*D*!.

Numerically however, it can be seen that the estimated PE for any finite set of white noise elements will almost never be 1. In fact, if the number of data subsets is not an even multiple of *D*!, the probability of *H*<sub>*est*</sub> = 1 is exactly zero. It is possible to rapidly simulate many realisations of Gaussian white noise, from which a *distribution* of estimated PEs can be constructed. An example of such a distribution is shown below for 100,000 realisations of GWN (*D* = 3, *N* = 3,333).

<p align="center">
  <img width="560" height="400" src="/Images/GWNPE.jpg">
</p>

The mean value of this distribution can be determined by expressing the logarithm term in equation 2.1 as a Taylor series expansion of *count deviations* [4], defined as

(3.1) &Delta;*q*<sub>*i*</sub> = *Np*<sub>*i*</sub> - E[*Np*<sub>*i*</sub>],

where E denotes the expectation value of the quantity inside the square brackets. Equation 2.1 then takes the form

*H*<sub>est</sub> = 1 - &alpha;<sub>1</sub>&Delta;(*q*<sub>*i*</sub>)<sup>2</sup> - &alpha;<sub>2</sub>(&Delta;*q*<sub>*i*</sub>)<sup>3</sup> - ...

To a first approximation, the expectation of equation [3] is

(3.2) E[*H*<sub>est</sub>] = 1 - (*D*!-1)/(2*N*log(*D*!)).

This is derived from the fact that &Delta;*q*<sub>*i*</sub> (and *p*<sub>*i*</sub>) are *multinomially distributed* variables, and that the second moment can be expressed as E[&Delta;*q*<sub>*i*</sub><sup>2</sup>] = *p*<sub>*i*</sub>(1-*p*<sub>*i*</sub>) = 1/*D*!(1-1/*D*!). The expectation of the distribution of *H*<sub>est</sub> can be viewed as the "true" *H* (i.e. 1) minus a deviation term. The presence of this deviation term is indicative of the fact that equation 2.1 is a biased estimator. For our above example of (*D* = 3, *N* = 10,000), the deviation term calculated in equation 3.2 agrees with the numerically calculted deviation to within a few percent. It is also possible to compute the variance in a similar fashion by truncating equating 1.4 and taking the expectation of the square, however this process is a little more involved, and so the reader is referred to [4] for details.

Instead, let us look at the *distribution* of PE values. The observed distribution is very close to a *chi-squared* (&chi;<sup>2</sup>) distribution, except reflected, scaled and shifted compared to the standard &chi;<sup>2</sup> distribution. Is is not difficult to see however that the deviation term in equation 3.2 follows an orthodox, scaled &chi;<sup>2</sup> distribution. We can see how this emerges from equation 3.2 by considering that &chi;<sup>2</sup> distributions arise from the sum of squared, normally-distributed random variables. With some rearrangement, the first two terms of equation 3.1 can be expressed as

(3.3) H = 1 - 1/(2*N*log(*D*!) * (1 - 1/*D*!) * &sum; *u*<sub>*i*</sub><sup>2</sup>,

where *u* is a random variable of unit variance and will converge to Gaussian random variables in the limit of large *N* (a general property of multinomially distributed variables). This corresponds to exactly *D*!-1 independent variables (the *D*!th variable is not considered independent, as it is constrained by the fact that the &sum;&Delta;*q*<sub>*i*</sub> = 0). We can directly read off equation 3.3 that the &chi;<sup>2</sup> distribution has *D*!-1 degrees of freedom and a *scale parameter* of 2*N*log(*D*!).

### PE of Brownian noise
Brownian noise (often referred to as *random walks*) in its discrete form is defined as noise whereby the difference in consecutive elements follows a GWN distribution. We can think of GWN as the *derivative* of Brownian noise in the continuous limit. Brownian noise as an example of a non-iid distribution, since the distribution of a given point depends on the sum of sample values that preceded it. Alternatively we say that it is *non-stationary*. Brownian noise is also an example of a *Markov* process, defined as a system where the output probability depends only on the previous state. Markov processes play an important role in statistics, as they describe processes with strong short-range correlations.

While Brownian noise is non-iid, the encoded symbol string that is produced via ordinal analysis *is* iid. For *D* = 3, the set of probabilities is {1/4, 1/8, 1/8, 1/8, 1/8, 1/4}, which can be reasoned out by noting that there is always a 50% chance for a point to be greater than (or less than) the previous point. Thus the probability of obtaining patterns 1<2<3 and 3<2<1 is (1/2)<sup>2</sup>. The remaining probabilities must be equal by symmetry.

<p align="center">
  <img width="560" height="400" src="/Images/BNPE.jpg">
</p>

Above is a histogram of estimated PEs from 100,000 realisations of Brownian noise with *D* = 3 and *N* = 3,333 as before. We know that the "true" PE of this system is -1/log(6)(0.5log(1/8) + 0.5log(1/4)) = 0.96713 (to 5 decimal places). The mean of the above distribution is 0.96672 = 0.96713 - 1/2/3,333/log(6) = 0.96713 - 0.00042, and so is consistent with equation 3.2.

A natural question that arises at this point is "can we distinguish between white noise and Brownian noise through the PE"? The answer is unquestionably yes, in fact, we can do this without computing the PE at all - we'd only need to count the instances of each pattern. Distinguishing things via PE is nontheless useful because it reduces a multi-dimensional problem to a 1D one. The next question is "with what certainty?", which taps into the notion of *inference*. It is clear looking at the probability distributions above that they are well separated in that the difference in expectation values is well in excess of their combined standard deviations, however this is only true for noise series of this specific length (10,000 points, or 3,333 ordinal pattern subsets). Shorter data series will yield larger standard deviations, and so the ability to distinguish one from of noise from the other will be reduced. Conversely, more data points is usually favourable for the purpose of inference. We will look at inference testing in a bit more detail later.

### PE analysis using Bayesian methods

To this point, we have looked at probability distributions by running repeated noise realisations over and over to build up a histogram. In other words, we computed the probability of obtaining *H*<sub>est</sub> from some known set of probabilities (with a corresponding "true* *H*). In statistics parlance, we computed *p*(*H*<sub>est</sub>|*H*), where the | symbol denotes a conditional probability, and is read as "the probability of observing *H*<sub>est</sub> given *H*".

The problem with this approach is that *H* is often not known *a priori*, so trying to compute histograms is iffy at best. Thankfully Bayes' rule, given by

(3.4) *p*(*A*|*B*) = *p*(*B*|*A*)*p*(*A*)/*p*(*B*),

allows us to frame the problem in an alternative way. Instead of seeking the distribution of *H*<sub>est</sub> given some *H*, we instead seek to find the distribution of *H* given some *H*<sub>est</sub>, flipping the problem "on its head". The advantage of this approach (apart from being more intuitive, since *H* is the unknown thing we are looking to characterise), it does not require the generation of multiple data realisations. A single data set, with a single estimate of *H* is sufficient to generate the required probability distribution.

Examining each term in equation 3.4 in turn;
* *p*(*A*|*B*) is called the *posterior distribution* and is the thing we want to calculate.
* *p*(*B*|*A*) is called the *likelihood function*, since a) it is not a probability distribution as it is not normalised, and b) it expresses the relative probability (likelihood) of obtaining the observed *B* given some *A*.
* *p*(*A*) is called the *prior distribution*, and is noteworthy in that it does not depend on *B* at all. It is the probability distribution of *A* *before* we have considered any observations. It can seem counter-logical at first to sensibly define a probability distribution of this type. In Bayesian statistics probability distributions are interpreted as a "confidence of belief" rather than an expected frequency of appearance over many realisations, thus the prior distribution essentially codifies all prior information/expectations of what *A* is likely to be before any symbols are observed. The logic of prior distributions is a topic that inspires plenty of debate that we will not rehash here. Note though that it is often possible to set *p*(*A*) in such a way as to delegate all influence over the posterior distribution to the likelihood function (thus allowing the data to "speak for itself").
* *p*(*B*) can be thought of as a normalising term that ensures that the RHS is a proper probability distribution.

To simplify the evaluation of the posterior distribution, it is common practice to assume a *conjugate prior* distribution in order to make the numerator in 3.4 more tractable. For this reason, we initially set *A* to **P** (the vector containing the marginal probabilities {*P*<sub>1</sub>,...,*P*<sub>*D*!</sub>} rather than *H*, because the likelihood function can't be expressed in a closed, analytic form. We set *B* to be **O**, the vector of symbols that have been observed. In this case, the likelihood function becomes

(3.5) *p*(**O**|**P**) = &prod;<sub>*i*=1</sub><sup>*D*!</sup> *P*<sub>*i*</sub><sup>*ni*</sup>,

where *ni* are the number of occurances of the *i*th symbol appearing in the vector **O** (thus *ni* are the *sufficient statistics*). It is important to reaffirm that the *P*<sub>*i*</sub> in equation 3.5 **are now variables**, not fixed/known quantities as in previous equations. Equation 3.5 arises from **P** being a multinomial distribution, and it is here that the iid assumption is important - if this assumption did not hold, then this would not be the correct likelihood function and *ni* would no longer be sufficient statistics for calculating the likelihood.

For a likelihood function of this form, it is common to choose a prior in the form of a *Dirichlet* distribution [5], defined as

(3.6) *p*(**P**) = (1/**B**)&prod;<sub>*i*=1</sub><sup>*D*!</sup> *P*<sub>*i*</sub><sup>&alpha;*i*-1</sup>,

where &alpha;<sub>*i*</sub> are the *hyperparameters* of the distribution and **B** is the Beta function, a normalising factor which is commonly defined in terms of gamma functions as

(3.7) **B** = &Gamma;(&alpha;<sub>1</sub>)...&Gamma;(&alpha;<sub>*D*!</sub>)/&Gamma;(&alpha;<sub>1</sub> + ... + &alpha;<sub>*D*!</sub>).

The convenience of this prior is evident when examining its product with the likelihood function. This product is proportional to the posterior distribution, hence we can see that

(3.8) p(**P**|**O**) &alpha; &prod;<sub>*i*=1</sub><sup>*D*!</sup> *P*<sub>*i*</sub><sup>&alpha;*i* + *ni* -1</sup>,

which is another Dirichelt distributions, where the observations *ni* have been added to the hyperparameters. For this reason, the hyperparameters are interpreted as pseudo/equivalent observations of their respective symbols. Equation 3.8 is easily normalised by identifying the correct Beta function from the modified hyperparameters. 

#### Calculating Moments of p(H)
Although expressing distributions in terms of **P** is analytically convenient for the purposes of constructing the likelihood function and the prior, it is inconvenient for numerical computation due to the high dimensionality of **P**. It is thus desirable to perform numerical computations in *H* rather than **P**. One way to determine *p*(*H*) without directly invoking *p*(**P**) is through the calculation of *moments*.

In statistics, the *n*th moment of a distribution of some variable is defined as the expectation of the *n*th power of the variable, i.e.

(3.9) E(*x*<sup>n</sup>) = 1/*N* &sum; *x*<sup>*n*</sup>*P*<sub>*x*</sub>,

converging to an integral on the RHS in the case of a continuous variable where *P*<sub>*x*</sub> becomes *p*(*x*). For the sake of completeness, it is worth remarking that for moments of second order or higher, it is common to also define a *centralised moment* where

(3.10) E((*x*-E(*x*))<sup>n</sup>) = 1/*N* &sum; (*x*-E(*x*))<sup>*n*</sup>*P*<sub>*x*</sub>.

The variance is the second *centralised* moment and **not** the ordinary (or *raw*) second moment. Higher order centralised moments are the skewness (3rd centralised moment) and the kurtosis (4th centralised moment).

Here, we seek to exploit the fact that for variables on closed, finite intervals, the moments uniquely define the functional form of the probability distribution. This is known as the Hausdorff moment problem [6]. Although the number of moments is countably infinite, we can often compute a good approximation with only the lower order moments; analgous to truncating a Taylor series. The *n*th moment of *p*(*H*) is defined as

(3.11) E(*H*<sup>*n*</sup>) = &int;<sub>0</sub><sup>1</sup> *H*<sup>*n*</sup> *p*(*H*) d*H*.

We can compute the moments of *p*(*H*) through a change of variable from *H* to **P**. Since probability distributions must remain normalised, expectations must be invariant under changes of variable, hence

(3.12) E(*H*<sup>*n*</sup>) = &int;<sub>*S*</sub> *H*<sup>*n*</sup> *p*(**P**) d**P**,

where *S* is the simplex defined by the condition &sum; *P*<sub>*i*</sub> = 1. Substituting the form of the prior in for *p*(**P**), we are now in a position to begin computing moments of the prior *p*(*H*). Similarly, we can find *p*(*H*|**O**). Because our judicious choice of prior resulted in the prior and posterior sharing the same form, they are tantamount to the same problem. Furthermore, the posterior distribution from one set of observations can be subsequently utlised as the prior distribution of the next set of observations with little difficulty.

##### First Moment
It can be seen from direct substitution of equation 2.2 into 3.12 for *H* we arrive at a series of integrals of the form

(3.13) &int; 1/**B** log(*P*<sub>*i*</sub>) &prod;<sub>*i*=1</sub><sup>*D*!</sup> *P*<sub>*i*</sub><sup>&alpha;*i*+*kj*-1</sup> d**P**,

where we have substituted the prior distribution (equation 3.6) for *p*(**P**). Here, *kj* is an index term that is 1 for *i* = *j* and 0 otherwise. The trick to solving these integrals is realising that the log terms can be eliminated by expressing the integrand as a derivative with respect to the &alpha; hyperparameters

(3.14) 1/**B** &int; &part;/&part;&alpha;<sub>*j*</sub> &prod;<sub>*i*=1</sub><sup>*D*!</sup> *P*<sub>*i*</sub><sup>&alpha;*i*+*kj*-1</sup> d**P**.

Exchanging the order of the derivative and integral allows these integrals to be solved purely through the normalisation condition of the Dirichlet distribution. With some algebra, the first moment of the prior *p*(*H*) can be expressed as

(3.15) E(*H*) = -1/log(*D*!) &sum; &Gamma;(&alpha;<sub>0</sub>)/&Gamma;(&alpha;<sub>*i*</sub>) &part;/&part;&alpha;<sub>*i*</sub> &Gamma;(&alpha;<sub>0</sub> + 1)/&Gamma;(&alpha;<sub>*i*</sub> + 1),

which is equivalent to

(3.16) E(*H*) = -1/log(*D*!) &sum; &alpha;<sub>i</sub>/&alpha;<sub>0</sub> (&psi;(&alpha;<sub>i</sub> + 1) - &psi;(&alpha;<sub>*0*</sub> + 1)),

where the summations are over the indices *i* = 1 to *D*!. Here, &psi; denotes the *digamma* function, which is available in most numerical computing libraries, meaning that the expectation of *p*(*H*) and *p*(*H*|**O**) can be rapidly computed knowing only the hyperparameters (which are typically user-defined) and the observation counts *ni*. With this approach, there is no need to numerically model anything directly in the **P** space.

For completeness, we can utilise the approximations

(3.17) &psi;(*x*) = log(*x*) - 1/2*x*

(3.18) log(*x*+1) - log(*x*) = 2/(2*x*+1)

to show this expression for E(*H*) for the posterior distribution converges to equation 3.2 in the limit of large *N*. With these approximations, equation 3.16 can be expressed as

(3.19) E(*H*) = -1/log(*D*!) &sum; &alpha;<sub>i</sub>/&alpha;<sub>0</sub> (log(&alpha;<sub>i</sub>/&alpha;<sub>0</sub>) + 2/(2&alpha;<sub>i</sub>+1) - 1/(2&alpha;<sub>i</sub>+2) - 2/(2&alpha;<sub>0</sub>+1) + 1/(2&alpha;<sub>0</sub>+2)).

Equation 3.2 tacitly assumes a uniform prior distribution, so &alpha;<sub>*i*</sub> = *ni* and &alpha;<sub>0</sub> = *N*. Thus the first term in E(*H*) is just *H*<sub>est</sub>, and the remaining terms make up the deviation term. With some algebra, it is not too difficult to show that equation 3.19 in this limit becomes

(3.20) E(*H*) = *H*<sub>est</sub> - 1/2*N*log(*D*!) &sum; 1 - *ni*/*D*!,

which is equivalent to equation 3.2.

#### Second Moment
For the second moment, we must square the expression for *H* to determine the integrals that need to be solved. A handy way to visualise this product is to view each term of *H* as the element of a vector, **V**. The square of *H* can then be viewed as the sum of all the terms in the outer-product **V**<sup>T</sup>**V** (represented as *V*<sub>*ij*</sub><sup>*j*</sup> in Einstein notation). The diagonal terms in **V**<sup>T</sup>**V** are

(3.21) 1/**B** &int; &part;<sup>2</sup>/&part;&alpha;<sub>*j*</sub><sup>2</sup> &prod;<sub>*i*=1</sub><sup>*D*!</sup> *P*<sub>*i*</sub><sup>&alpha;*i*+2*kj*-1</sup> d**P**,

while the cross-terms are represented by

(3.22) 1/**B** &int; &part;<sup>2</sup>/&part;&alpha;<sub>*j*</sub>&part;&alpha;<sub>*m*</sub> &prod;<sub>*i*=1</sub><sup>*D*!</sup> *P*<sub>*i*</sub><sup>&alpha;*i*+*kj*+*km*-1</sup> d**P**,

where *i* &ne; *m*. Following the same process as before, namely exchanging the order of integral and differntiation, using the normalisation condition of a Dirichlet distribution, and then finding the derivatives, E(*H*<sup>2</sup>) can be expressed as

(3.23) E(*H*<sup>2</sup>) = 1/log(*D*!)<sup>2</sup> &sum;<sub>*i*</sub> &alpha;<sub>*i*</sub>(&alpha;<sub>*i*</sub>+1)/&alpha;<sub>*0*</sub>(&alpha;<sub>*0*</sub> + 1) ((&psi;(&alpha;<sub>*i*</sub>+2) - &psi;(&alpha;<sub>*i*</sub>+2))<sup>2</sup> + &psi;<sub>1</sub>(&alpha;<sub>*i*</sub>+2) - &psi;<sub>1</sub>(&alpha;<sub>*0*</sub>+2)) + 1/log(*D*!)<sup>2</sup> &sum;<sub>*i*&ne;*j*</sub> &alpha;<sub>*i*</sub>&alpha;<sub>*j*</sub>/&alpha;<sub>*0*</sub>(&alpha;<sub>*0*</sub> + 1) ((&psi;(&alpha;<sub>*i*</sub>+1) - &psi;(&alpha;<sub>*i*</sub>+2)) (&psi;(&alpha;<sub>*j*</sub>+1) - &psi;(&alpha;<sub>*i*</sub>+2)) - &psi;<sub>1</sub>(&alpha;<sub>*0*</sub>+2)),

where &psi;<sub>1</sub> is the *trigamma* function. The gamma, digamma and trigamma functions are the first in a series of functions generated by taking repeated derivatives over the gamma function. With the second moment, the variance can be computed as E(*H*<sup>2</sup>) - E(*H*)<sup>2</sup>. As with the first moment, it can be shown to converge to the variance obtained in the limit of large *N* from the multinomial statistics.

#### Higher-order Moments
This procedure can be generalised to any higher order moment with some notational tweaks.

(3.24) E(*H*<sup>*n*</sup>) = (-1/log(*D*!))<sup>*n*</sup> &sum; 1/**B** &int; &part;<sup>*n*</sup>/&part;**&alpha;**<sup>**K**</sup> &prod;<sub>*i*=1</sub><sup>*D*!</sup> *P*<sub>*i*</sub><sup>&alpha;*i*+*Kj*-1</sup> d**P**,

where **K** is a *multi-index* where the sum of the index components (denoted |**K**|) is *n* for the *n*th moment. The summation in equation 3.24 is over all the possible **K** satisfying |**K**| = *n*. Here, &part;<sup>*n*</sup>/&part;**&alpha;**<sup>**K**</sup> is interpreted as the combined *Ki*th derivative with respect to &alpha;<sub>*i*</sub>. For example, a multinomial index of {3,0,0,0,0,0} (for *D* = 3) would correspond to taking the third derivative w.r.t. &alpha;<sub>*1*</sub>, {1,0,1,1,0,0} would correspond to &part;<sup>*3*</sup>/&part;&alpha;<sub>1</sub>&part;&alpha;<sub>3</sub>&part;&alpha;<sub>4</sub> and so on. Analytic expressions up to the fourth moment have been computed as part of this package.

## 4. PE statistics for non iid-distributed symbols
It is important to re-iterate that everything in section 3 is applicable only when the observed symbols are iid. In this section, we will relax this requirement and look at cases where observed symbols may be non-iid. An important concept in this context is *correlation*. By definintion, iid data is uncorrelated (independent), so this is not something that needed considering in the previous section. In essence, correlation is the tendency for something to following another thing. You've probably heard of correlation in the context of determining relationships between quantities/variables; if changes in one variable tends to instigate change in another variable, they are said to be correlated.

In the context of PE analysis, correlation and the tendency for one event to follow another can be summarised in terms of conditional probabilities. We encountered conditional probabilities in the previous section in the context of Bayesian analysis, where *p*(*X*|*Y*) is the probability of observing *X* given some observation *Y*. This can be applied to the observation of symbols as well. Most generally, we can express the probability of observing a symbol *S* in a non-iid system as *P*(*S*|**O**), where **O** denotes all the observations made up to the observation of *S*. For iid systems *P*(*S*|**O**) = *P*(*S*), that is, prior observations have no influence over future probabilities.

Clearly *P*(*S*|**O**) defines an enormous-dimensional probability space, so it is normal to constrain **O** to a limited subset of prior observations. Obviously, the simplest approach is to constrain **O** to a single value. This defines a class of system known as a *Markovian* (or memory-less) system, where the next state/observation depends only on the current state/observation. We will assume for brevity that **O** is the prior observation, *S*<sub>i</sub> and we seek to characterise the probability of observing *S*<sub>i+1</sub>.

Under the Markov assumption, the probabilities associated with the different state transitions can be characterised as a *vector* of probability distributions, {**P**<sub>1</sub>,...,**P**<sub>*D*!</sub>}, where **P**<sub>*i*</sub> represents the probability distribution when *S*<sub>*i*</sub> is the *i*th symbol. Alternatively, they can be represented as a *matrix* where *P*<sub>*i*|*j*</sub> represents the probability of observing the *i*th symbol given the previous observation was the *j*th symbol. The vertical bar is included to distinguish it from the *joint* distribution *P*<sub>*ij*</sub>, which is the probability of observing an ith and jth symbol together in a pair of observations. In short, non-iid distributions require a higher-dimensional characterisation of the probabilities. 

### PE analysis of overlapping subsets of white noise
As an instructive exercise, let us now look at the effect of using overlapping subsets when dividing up our original ordinal data set. Recall that it was recommended not to allow these subsets to overlap as it introduces correlations in resultant symbol string. We will now analyse this effect for white noise data sets.

Let us first consider the simple case where *D* = 2. There are two possible ordinal patterns, which we shall label subsets where *x*<sub>*i*</sub> < *x*<sub>*i+1*</sub> with the symbol "0" and *x*<sub>*i*</sub> > *x*<sub>*i+1*</sub> with the symbol "1" (remember, we rely on tiebreaking to avoid equalities). A clever way to determine the matrix of conditional probabilities is to decompose *D* = 3 inequalities into consecutive *D* = 2 inequalities like so;
* {*x*<sub>1</sub> < *x*<sub>2</sub> < *x*<sub>3</sub>} ⇒ {*x*<sub>1</sub> < *x*<sub>2</sub>}, {*x*<sub>2</sub> < *x*<sub>3</sub>} = {0, 0}
* {*x*<sub>1</sub> < *x*<sub>3</sub> < *x*<sub>2</sub>} ⇒ {*x*<sub>1</sub> < *x*<sub>2</sub>}, {*x*<sub>2</sub> > *x*<sub>3</sub>} = {0, 1}
* {*x*<sub>2</sub> < *x*<sub>1</sub> < *x*<sub>3</sub>} ⇒ {*x*<sub>1</sub> > *x*<sub>2</sub>}, {*x*<sub>2</sub> < *x*<sub>3</sub>} = {1, 0}
* {*x*<sub>2</sub> < *x*<sub>3</sub> < *x*<sub>1</sub>} ⇒ {*x*<sub>1</sub> > *x*<sub>2</sub>}, {*x*<sub>2</sub> < *x*<sub>3</sub>} = {1, 0}
* {*x*<sub>3</sub> < *x*<sub>1</sub> < *x*<sub>2</sub>} ⇒ {*x*<sub>1</sub> < *x*<sub>2</sub>}, {*x*<sub>2</sub> > *x*<sub>3</sub>} = {0, 1}
* {*x*<sub>3</sub> < *x*<sub>2</sub> < *x*<sub>1</sub>} ⇒ {*x*<sub>1</sub> > *x*<sub>2</sub>}, {*x*<sub>2</sub> > *x*<sub>3</sub>} = {1, 1}

Since each of these six possibilities must be equal by symmetry, we can work out the conditional probabilities to be
 **P**|0|1|
---|---|---
0|1/3|2/3
1|2/3|1/3

In words, we deduce that the likelihood of observing a "0" is twice as likely when the previous observation was a "1" than when the previous observation was a "0" and vice versa, a non-negligible effect! for *D* = 2 the conditional probabilities are nice and symmetric, and in fact, we can reduce this to an iid system through a second XOR encoding step, where "1" corresponds to a change in observed pattern and "0" if the pattern is repeated, which would mean that *P*<sub>0</sub> = 1/3 and *P*<sub>1</sub> = 2/3. Unfortunately this is generally not the case for higher *D*, as the following conditional probabilities for *D* = 3 shows
 **P**|0|1|2|3|4|5
---|---|---|---|---|---|---
0|6/40|6/40|3/40|6/40|6/40|13/40
1|6/40|6/40|13/40|6/40|6/40|3/40
2|6/40|6/40|3/40|6/40|11/40|8/40
3|13/40|3/40|4/40|3/40|8/40|9/40
4|3/40|13/40|9/40|8/40|3/40|4/40
5|6/40|6/40|8/40|11/40|6/40|3/40

where the labels have been arranged such that the pairs {0,1}, {2,3} and {4,5} correspond to ordinal patterns with the same central index, i.e. *x*<sub>1</sub> < *x*<sub>2</sub> < *x*<sub>3</sub> and *x*<sub>3</sub> < *x*<sub>2</sub> < *x*<sub>1</sub> are grouped in a pair. Perhaps surprisingly, the matrix of conditional probabilities is not symmetric and with a fair deal on non-trivial structure to it, so the simple tricks used to simplify the system for *D* = 2 don't work.







#### References
1. C. Bandt and B. Pompe, Phys. Rev. Lett. 88, 174102 (2002).
2. C.E. Shannon, W. Weaver, The Mathematical Theory of Communication (1949), Univ of Illinois Press.
3. D.J. Little, D.M. Kane, Phys. Rev. E 95, 052126 (2017).
4. D.J. Little, D.M. Kane, Phys. Rev. E 94, 022118 (2016).
5. C. Bishop, Pattern Recognition and Machine Learning *Chapter 2*, Springer (2006).
6. F. Hausdorff, Math. Z 9, 74 (1921).
