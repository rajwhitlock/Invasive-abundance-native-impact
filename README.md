# Invasive-abundance-native-impact


https://zenodo.org/badge/176828931.svg


R-scripts for meta-analysis of the impacts of increasing invasive species' abundance on native animal and plant populations and communities

These scripts support meta-analysis of the relationship between increasing invasive species abundance and native species' responses at the population level and the community level. The scripts are designed for analysis of a dataset synthesising invasive species impacts in a wide range of study taxa and habitat types globally. Datasets referred to in the script are published separately, and the location of these files will be given here on publication. The meta-analysis code presented in this repository is, together with the seperately published raw data, sufficient to replicate analyses presented in the following manuscript:

Bradley, B. A., Laginhas, B. B., Whitlock, R., Allen, J. M., Bates, A. E., Bernatchez, G., Diez, J. M., Early, R., Lenoir, J., Vilà, M. & Sorte, C. J. B. 2019. Disentangling the abundance-impact relationship for invasive species. In review.

A key aim of the meta-analysis presented here is to summarise the shape of the invasive abundance – native response relationship, capturing both linear effects and non-linear effects in the form of second order polynomials. Two meta analysis scripts are presented. In the first, effect sizes are computed from partial correlation coefficients, separately for the linear and polynomial regression coefficients (see AvI_metaanalysis_partial-r_24_03_19.R). In the second, linear and polynomial regression coefficients are used directly as effect sizes, following data rescaling (see AvI_metaanalysis_slopes__24_03_19.R). Code is provided to conduct meta-analysis using Baysian mixed effects models, and to plot the results as bivariate forest plots, and as invasive abundance–native response curves, with 95% credible zone. Code is also provided to create funnel plots.

The script Bar_plot_error_bars_24_03_19.R provides utility functions for plotting error bars/whiskers about plotted points in an R plot.
The script AvI_sensitivity_analyses_24_03_19.R provides additional sensitivity analyses, described in the manuscript above, and that support the meta-analyses.
