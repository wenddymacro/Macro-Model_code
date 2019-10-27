This archive contains Matlab code and data to replicate the methods presented in #1527 . The software comes with no warranty, and for any further problems or questions contact Hippolyte.d-albis@univ-paris1.fr.

Replication material:

Age Groups and the Measure of Population Aging 
	
	Hippolyte d’Albis (Paris School of Economics, University Paris 1)
	Fabrice Collard(Department of Economics, University of Bern)

All codes are written in Matlab. They were first developed using Matlab R2008, but should be compatible with the lastest version too.

Figure_2_3.m replicates figures 2 and 3 of the paper.

Figure_5_8.m replicates figures 5 through 8 of the paper. This file makes use of the function optimal_grouping_c.m which computes our optimal groups using the method developed by Aghevli and Mehran (1981).

Figure_9.m replicates figure 9 and generates Table 1.

table_2.m produces the table 2 in the paper in a LaTeX form. It makes us of the econometrics toolbox developed by James Lesage and available from his web page. Note that this file assumes that the function main_optigroup.m was run first. This function produces various indicators and saves them into files that that are used in the table_2.m file. Also note that these files are present in the data directory.

All needed data are in the data directory and should be readable in any text editor. 