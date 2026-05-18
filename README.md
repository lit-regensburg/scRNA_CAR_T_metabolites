# Prerequisites

* direnv must be installed and allowed
* Snakemake must be present on the system for running as a pipeline (otherwise, running `rmarkdown::render()` on the individual scripts is also possible)
* A Singularity container containing the dependencies must be present for running the Snakemake workflow
* *container_setup.yaml* must be adapted to the actual container location and location of reference files
* *config.yaml* must be adapted to the location of reference files and location of Cellranger files (path specified relative to the .projectroot file)

# Running the pipeline

* Navigate to *experiment/scRNA/analysis/hep62703_2025-08-13_pub/scripts*
* Then run the Snakemake pipeline with the command `snakemake --cores 6 -s Snakefile_scRNA.py` (you can 
specify the `cores` parameter as you like, use the `n` parameter for a dry run)
* The pipeline graph in *experiment/scRNA/analysis/hep62703_2025-08-13_pub/pipeline_graph.pdf* displays the workflow with the individual rules/scripts that are run
* To run an individual script, use the command `rmarkdown::render("script_name", params = list(config = "config.yaml"))` from within an R session (Singularity container) in the directory *experiment/scRNA/analysis/hep62703_2025-08-13_pub/scripts*


