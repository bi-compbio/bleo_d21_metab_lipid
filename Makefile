RSCRIPT=Rscript

# tasks
all:: wrangling descriptive diff plots

wrangling: 01_metab_se.html 02_lipid_se.html 03_mae.html 
descriptive: 10_descriptive.html 
diff: 20_differential_analysis.html 21_pathway_analysis.html 25_network_analysis.html 
plots: 30_differential_plots.html 31_differential_age_plots.html

# rmd files
%.html: %.Rmd; $(RSCRIPT) -e "require(rmarkdown); render('$<');"
