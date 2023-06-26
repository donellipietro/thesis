all: visualize_results compare_results latex

generate_data:
	Rscript scripts/functions/*
	Rscript scripts/generate_data.R
	Rscript scripts/generate_data_comparison.R

visualize_results:
	Rscript scripts/visualize_results.R
	
compare_results:
	Rscript scripts/compare_results.R
	
latex: visualize_results compare_results
	pdflatex export_results.tex
	evince export_results.pdf
	
clean:
	rm -f -r *.log *.aux *.fdb_latexmk *.fls *.synctex.gz
	
distclean: clean
	rm -f *.pdf
	rm -f scripts/functions/*.RData