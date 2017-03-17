# ariba-publication
Supplementary data and scripts for the ARIBA publication.


## How to rerun the analysis

Build the Docker container:

    docker build -t ariba-publication .

Run the analysis:

    docker run --rm -it -v $PWD:/data ariba-publication

Then the main figures can be found in `Figures/`, the
supplementary tables file is `supplementary_tables.xlsx`,
and the supplementary PDF file is `Supplementary_LaTeX/main.pdf`.

