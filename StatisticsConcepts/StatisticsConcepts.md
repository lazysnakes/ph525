Statistical Concepts in PH525
========================================================

# Central Limit Theorem

This is introduced in Week 3's Inference module, videos 3 and 4.

## Background

The problem is: we have two populations that we want to compare. We can't measure the entire population, so we must measure a sample. However, there is sampling bias. How do we quantify sampling bias?

* Population Averages: $\mu_{X} = \frac{1}{m} \sum_{i=1}^m x_{i}$
* Observations: $X_{1},\cdots,X_{M}$
* Sample Averages: $\overline{X} = \frac{1}{M} \sum_{i=1}^M X_{i}$ -- the samples averages are
'Random Varialbles'
* Population Standard Deviation: $\sigma_{X}^2 = \frac{1}{m} \sum_{i=1}^m (x_{i} - \mu_{x})^{2}$

## Definition

The sample average follows a normal distribution centered at the population average with a standard deviation equal to th population standard deviation divided by the sample size, $M$.

*Central Limit Theorem*

$\overline{X} \sim N(\mu_x,\frac{\sigma_{x}}{\sqrt{M}})$

## T-Test

Is there a difference between two populations? How can we predict based on measurements of a sample of each group? The CLT tells us that:

$\frac{\overline{X} - \overline{Y}}{\sqrt{\frac{\sigma_{X}^2}{M}+\frac{\sigma_{Y}^2}{N}}} \sim N(0,1)$

The difference between the sample averages also follows a normal distribution and the width of the distribution is inversely proportional to the sample sizes.

*Note:* There is a problem with this formula: we don't know the population standard deviations.

We can use the sample variance to estimate the population standard deviation:

$s_{X}^2 = \frac{1}{M-1} \sum_{i=1}^M(X_{i} - \overline{X})^{2}$

This gives the *T-Statistic*:

$\frac{\overline{X} - \overline{Y}}{\sqrt{\frac{s_{X}^2}{M}+\frac{s_{Y}^2}{N}}} \sim N(0,1)$

For the T-Statistic to work, the sample sizes must be large. If they are not large, then the T distribution can be substituted for the normal distribution and you say that you have $M + N - 2$ degrees of freedom.

## Confidence Intervals

Confidence intervals gives 

$\overline{X} \pm 2 \frac{s_{X}}{\sqrt{M}}$


