FROM amancevice/pandas:1.5.3

RUN pip install synapseclient==4.8.0 \
                boto3~=1.38.0 \
                git+https://github.com/larssono/bridgeclient.git@b90eb4c \
                git+https://github.com/Sage-Bionetworks/synapsebridgehelpers@d22afb2 \
                git+https://github.com/redcap-tools/PyCap.git@3d4964b
RUN git clone https://github.com/Sage-Bionetworks/at-home-pd.git /root/at-home-pd

CMD ["/bin/bash", "-l"]