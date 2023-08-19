# Probabilistic-AHP

1. MAIN.m is main code

2. Input of problem is defined in problem.m. Input is matrix consist of number only at upper triangular part.
Each component of matrix could be cosist of several number from different opinion. This matrix is similar to pairwise comparison matrix
commonly used in AHP (Analytic Hierarchy Process) algorithm. The only different is that, we only need to input upper triangular part.

3. SOS (Symbiotic Organism Search) optimization algorithm is utilized to find maximum beta PDF (Probability Density Function) 
while being constrained by Consistency Ratio. The idea of this proposed AHP is to defined each element inside pairwise matrix as random varible 
and chosen number is highest possible PDF that still satisfy Consistency Ratio Limit.
