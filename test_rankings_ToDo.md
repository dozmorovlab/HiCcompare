- Use simulated matrices/a priori known differences.

- Add ROC curve analysis only for the ranking methods (no comparison with different normalizations). That is, the goal is to see which ranking method produces best ROC curve. May need to convert ranks to a scale used by the ROC curve function.

- Distance-weighting idea:  
    - Requires adapting $M$ and $A$ calculation and $MD\ plot$ visualization on the $log_{10}$ scale
    - Use ranks of $M$, or $A$, or $IF_1 - IF_2$, or a combination of them
    - Ranks are calculated for the whole matrix, the correspoiding distances are remembered
    - (Mean/max summaries of) ranks is downweighted (multiplied) by the distance weight $1 - (d_i/D)$, where $d_i$ is the discance corresponding to the rank measured at that distance, and $D$ is the full distance (matrix dimension)


- Amygdala vs. DLPFC comparison, so the positive $M$s are higher in DLPFC? 

- When using "top 5% of rnkDiff", the majority of $M$ differences are positive. That implies the raw difference is mostly one-directional, with interaction frequencies in one (DLPFC?) dataset are always higher across the range of distances. Need to investigate why.