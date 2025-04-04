# QUANDO
Quantitative analysis of DNA counterstains



%%   User-fiendly code based on the following article:
%$£   
%$£  "Location of oncogene-induced DNA damage sites revealed by quantitative analysis of a DNA counterstain"
%$£  by Greta Paterno' et al , European Biophysics Journal 2025


% OUTPUT text FILE contain: 
% First 3 columns:
% cell analyzed, DNA density, Colocalization with Heterochromatin
% The last 3 columns: 
% file analyzed, frame analyzed, cell # in the count mask


INPUT file (TIF):

-COUNT MASK: Count mask containing the numbered nuclei (can be obtained from ImageJ analyze particles)

-Intensity (2 channels): channel 1: DNA damage or other marker  - Channel 2: DNA counterstain 

-Binary mask corresponding to the selected marker (e.g. dna damage region or other marker)


