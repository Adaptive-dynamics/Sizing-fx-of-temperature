FROM docker.dragonfly.co.nz/squid-butterfish
MAINTAINER philipp@dragonfly.co.nz

RUN Rscript -e 'install.packages(c("cowplot","purrr"),repo="http://cloud.r-project.org/")'