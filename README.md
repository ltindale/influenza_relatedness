# Relatedness between HA1 residues in historic relative to contemporary influenza A(H3N2) viruses

[Link to interactive shiny influenza heatmap.](https://tindale.shinyapps.io/influenza_relatedness/)

## Description of shiny tabs:

### Heatmaps by influenza subtype

Heatmap of historic-contemporary match scores. The historic AA is the the most frequent AA for each combination of year and AA position. The contemporary AA is the current variant consensus sequence at the same position.

Scores are based on a normalized BLOSUM80 matrix where a match = 1.

### Heatmaps using blended scores

Blended Score = (historic vs vaccine score) + (historic vs contemporary clade score) + (contemporary clade vs vaccine score)

Historic-Contemporary-Vaccine Match and Mismatch were assigned discrete 1/0 values and added together to generate a blended 3-way theoretical score and to compare to BLOSUM80 matrix scoring in heatmaps.

0 = historic, contemporary and vaccine all mismatch

1 = only one combination matches (historic-contemporary or historic-vaccine or contemporary-vaccine)

3 = historic, contemporary and vaccine all match

### Historic amino acids

Bargraphs showing the percent of each amino acid as a total of all GISAID sequences at each amino acid residue by year.
