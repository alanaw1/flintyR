# flintyR 

This is the homepage of **flintyR**, the R version of the software **flinty** (*Fl*exible and *I*nterpretable *N*on-parametric *T*ests of Exchangeabilit*y*). The Python version is under development. 


![Bruno appears](bruno_intro.png)
<p align="center">
*Bruno is named after a famous statistican who studied exchangeability. Who might that [be](http://www.brunodefinetti.it/)?*
</p>

## What does this package offer?

**flinty** provides exact tests of exchangeability in multivariate datasets. 

- It is *non-parametric* (i.e., makes no distributional assumptions of the features), which makes it suitable for settings where the user might prefer not to make distributional assumptions about their data.   
- It is *flexible*, meaning that the practitioner can specify feature dependencies based on their knowledge of the problem. Our tests handle dependent features, should the dependencies satisfy partitionability. See tutorials for details.   
- It is *scalable*, so the user does not have to worry about the sample size $N$ or the number of features $P$ of the data. 
- It is *robust*, meaning that it controls for false positive rate (FPR) and remains powerful in realistic settings including uneven representation of subpopulations, sparsity of discriminative features, and small sample sizes.   

## Installation

To install our package, run

```
devtools::install.packages("flintyR")
```

## Example Usage

The code below demonstrates running our test on a binary matrix. 

```
# library(flintyR)
X <- matrix(nrow = 5, ncol = 10, rbinom(50, 1, 0.5))
getPValue(X) # perform exact test with 5000 permutations
# Output should be larger than 0.05
```

Examples involving real datasets can be found in the tutorials.

## Tutorials

We offer several tutorials on using our software.

- Conceptual explorations of exchangeability
- Application to avoiding risking double dipping
- Application to single cell genomics 
- Application to World Values Survey 

We love to see our methods and software used across multiple fields, so please reach out to us if you are interested in using them! If there is enough interest, we are happy to include more tutorials.    
