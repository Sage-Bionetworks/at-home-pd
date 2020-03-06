FROM amancevice/pandas:0.23.4-python3

RUN pip install synapseclient boto3 git+https://github.com/larssono/bridgeclient.git git+https://github.com/Sage-Bionetworks/synapsebridgehelpers.git
RUN git clone https://github.com/Sage-Bionetworks/at-home-pd.git --branch superusers --single-branch /root/at-home-pd

CMD python /root/at-home-pd/user_add/user_add.py
