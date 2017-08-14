# ariba-publication
Supplementary data and scripts for the ARIBA publication.


## How to rerun the analysis

Build the Docker container:

    docker build -t ariba-publication .

Run the analysis to reproduce the publication results and figures:

    docker pull martinghunt/ariba-publication:publication
    git clone https://github.com/martinghunt/ariba-publication.git
    cd ariba-publication
    docker run --rm -it -v $PWD:/data martinghunt/ariba-publication:publication

Then the main figures can be found in `Figures/`, the
supplementary tables file is `supplementary_tables.xlsx`,
and the supplementary PDF file is `Supplementary_LaTeX/main.pdf`.

