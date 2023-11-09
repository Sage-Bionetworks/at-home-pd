FROM rocker/tidyverse:4.0.3

RUN apt-get install libffi-dev
ENV CRYPTOGRAPHY_DONT_BUILD_RUST=true
RUN Rscript -e 'install.packages("devtools")'
RUN Rscript -e 'devtools::install_github("Sage-Bionetworks/synapser", ref="0.8.71")'
RUN Rscript -e 'devtools::install_github("philerooski/foxden", ref="0827600")'
RUN Rscript -e 'devtools::install_github("philerooski/bridgeclient", ref="8cf6c22")'
RUN git clone https://github.com/Sage-Bionetworks/at-home-pd.git /root/at-home-pd
