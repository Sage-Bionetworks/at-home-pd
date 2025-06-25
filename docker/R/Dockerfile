FROM rocker/tidyverse:4.1.0

RUN apt install libffi-dev
RUN Rscript -e 'devtools::install_github("philerooski/foxden", ref="0827600")'
RUN Rscript -e 'devtools::install_github("philerooski/bridgeclient", ref="8cf6c22")'
RUN Rscript -e 'devtools::install_github("Sage-Bionetworks/PythonEmbedInR", ref="0.9")'
RUN Rscript -e 'devtools::install_github("Sage-Bionetworks/synapser", ref="0.10.89")'

RUN git clone https://github.com/Sage-Bionetworks/at-home-pd.git /root/at-home-pd
