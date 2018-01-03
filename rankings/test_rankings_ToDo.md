+ Use simulated matrices/a priori known differences.

+ Add ROC curve analysis only for the ranking methods (no comparison with different normalizations). That is, the goal is to see which ranking method produces best ROC curve. May need to convert ranks to a scale used by the ROC curve function.

+ Distance-weighting idea:  
    + Use ranks of $M$, or $A$, or $IF_1 - IF_2$, or a combination of them
    + Ranks are calculated for the whole matrix, the correspoiding distances are remembered
    + (Mean/max summaries of) ranks is up-weighted (multiplied) by the distance weight $1 + (d_i/D)$, where $d_i$ is the discance corresponding to the rank measured at that distance, and $D$ is the full distance (matrix dimension)


- Amygdala vs. DLPFC comparison, so the positive $M$s are higher in DLPFC? 

- When using "top 5% of rnkDiff", the majority of $M$ differences are positive. That implies the raw difference is mostly one-directional, with interaction frequencies in one (DLPFC?) dataset are always higher across the range of distances. Need to investigate why.


+ Try to convert M to Z score and raw diff to Z score, take mean then multiply by distance weight see what we get

+ Try z scores for M with a cut off for raw difference so only z-scores with raw diff bigger than cut off are considered. Maybe also try weighting by distance for raw difference

- look at fit-hic see if we can just use their method and overlap results for 2 matrices

- try mean of z score for M and z score for A

- use smyth RWPE data as gold standard - make his regions into our regions for comparison

- decide whether to remove M values with low A before or after converting to Z scores
  - Z before removing low A: less significant values due to higher std dev
  - Z after removing low A: more significant values due to lower std dev
  - removing low A before Z might be better for multiple testing correction but produces more false positives

- implement distance adjusted multiple testing correction, test on smyth RWPE data


- grep -i hicdiff *.R in R folder, replace old instances in documentation
- Redo Vignette
- go through all functions to make sure unnecessary ones added for testing are removed


- go back to ranks
- do ROC to figure out A and M cutoffs for what to call significant
- use smyth data as gold standard
  - regions with 50% overlap with their significant regions for gold standard
  
  
- do 2 ROCs one at short distance and one at large distance
- spike in differences only at single distance at a time and run our method on it then generate ROC

