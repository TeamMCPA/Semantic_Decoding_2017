# Decoding semantic representations from fNIRS signals
Software for performing representational similarity analysis (RSA)-based decoding on fNIRS data. Demonstrated by Zinszer, Bayet, Emberson, Raizada, &amp; Aslin (2017) in *Neurophotonics* (forthcoming).

Analyses for each experiment have their own demo scripts:
(commentary provided within)

Semantic_8class_Neurophotonics_Exp1_30June2017.m 
Semantic_8class_Neurophotonics_Exp2_botharrays_30June2017.m
Semantic_8class_Neurophotonics_Exp2_latlarray_30June2017.m
Semantic_8class_Neurophotonics_Exp2_postarray_30June2017.m

The demo script relies on several support functions also included in this repository. fNIRS data are also included in the repository, but you should be able to run this script over your own data with minor modifications.

Permutation-based significance testing is included at the end of each script. This code can take a long time (several hours) to run, depending on the number of permutations performed. To find the results of each analysis and the saved output of the permutation tests we ran, look in the directory "analysis_results" and load the .mat file.

--------------------------------------------
See also our initial poster presentatation at SfNIRS 2016: http://benjaminz.com/SfNIRS-2016-demo/

Please see Emberson, Zinszer, Raizada, &amp; Aslin (2016) regarding Multichannel Pattern Analysis for fNIRS.

- bioRxiv: http://biorxiv.org/content/early/2016/06/30/061234

- Github: http://teammcpa.github.io/EmbersonZinszerMCPA/

Please see Anderson, Zinszer, &amp; Raizada (2015) regarding representaitonal similarity based decoding and encoding methods for fMRI, the foundation upon which this work was built.

- paywall access via NeuroImage journal: http://www.sciencedirect.com/science/article/pii/S1053811915011489

- free access via Raizada Lab at University of Rochester: http://www.bcs.rochester.edu/people/raizada/papers/AndersonZinszerRaizada_RepSimEncoding_NeuroImage2016.pdf
