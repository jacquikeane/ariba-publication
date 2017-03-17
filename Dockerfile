FROM ubuntu:16.04

RUN apt-get update \
  && apt-get install --no-install-recommends -y git build-essential libffi-dev mummer cd-hit zlib1g-dev wget gdebi-core python python3-dev python3-setuptools python3-pip texlive-full \
  && wget -q https://mirrors.ebi.ac.uk/CRAN/bin/linux/ubuntu/xenial/r-base-core_3.3.2-1xenial0_amd64.deb \
  && gdebi -n r-base-core_3.3.2-1xenial0_amd64.deb \
  && rm r-base-core_3.3.2-1xenial0_amd64.deb \
  && R -e 'install.packages(c("ggplot2"), repos="http://cran.r-project.org")' \
  && wget -q http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip \
  && unzip bowtie2-2.2.9-linux-x86_64.zip \
  && rm bowtie2-2.2.9-linux-x86_64.zip \
  && pip3 install ariba==2.8.1 cairosvg==2.0.0 openpyxl==2.4.1

RUN apt-get install -y python3-tk
RUN pip3 install scipy

ENV ARIBA_BOWTIE2=$PWD/bowtie2-2.2.9/bowtie2 ARIBA_CDHIT=cdhit-est

# this makes matplotlib work without X11, otherwise get the error
# _tkinter.TclError: no display name and no $DISPLAY environment variable
ENV MPLBACKEND="agg"

CMD /data/Scripts/run_all_analysis.sh
