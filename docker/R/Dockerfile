FROM rocker/tidyverse:3.6.1

RUN apt install libffi-dev
ENV CRYPTOGRAPHY_DONT_BUILD_RUST=true
RUN Rscript -e 'install.packages("devtools")'
RUN Rscript -e 'devtools::install_github("Sage-Bionetworks/synapser", ref="0.8.71")'
RUN Rscript -e 'devtools::install_github("philerooski/foxden")'
RUN Rscript -e 'devtools::install_github("philerooski/bridgeclient")'
RUN git clone https://github.com/Sage-Bionetworks/at-home-pd.git /root/at-home-pd
