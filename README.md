# Predective_Potential_FTemp_Indices
Scripts for generating results and figures for determining predictive potential of F-Temperature for thermodynamic properties.

**File Name	Description**

**Combined_Corr_Curves_FTemp.m**	Generates combined correlation curves plots, grouped by experimental data, for the data, highlighting the curves' peaks and intersection between two curves with highest correlation in some given alpha.

**Good_alpha_intervals_FTemp.m**	Generate correlation curve plots, one for each pair of data, highlighting the alpha intervals which results in good correlation.

**Scatter_plots_Ftemp.m**	Generate scatter plots, one for each pair of data.

**GoldenSectionSearch_Maximum.m**	IMPORTANT: This must be placed in the working directory as it is referenced in the three other scripts. This is the Octave 7.2 implementation of the Golden Section Search algorithm. I manually translated this from the Python code. 
