FROM kbase/kbase:sdkbase.latest
MAINTAINER KBase Developer
# -----------------------------------------

# Make sure SSL certs are properly installed
RUN apt-get install python-dev libffi-dev libssl-dev \
    && pip install pyopenssl ndg-httpsclient pyasn1 \
    && pip install requests --upgrade \
    && pip install 'requests[security]' --upgrade

RUN mkdir -p /kb/module/work && mkdir -p /kb/module/lib
RUN chmod 777 /kb/module
WORKDIR /kb/module

# Genbank uploader uses the Transform script_utils.py, so install that with
# dependencies to /kb/deployment/lib 
RUN git clone https://github.com/kbase/transform && \
    cd transform && \
    git checkout 3d0f47b && \
    cd /kb/module && \
    cp -a transform/lib/biokbase /kb/deployment/lib/ && \
    . /kb/module/transform/deps/pylib.sh

# Genbank uploader uses the data_api, so install that to our lib directory
RUN git clone https://github.com/kbase/data_api && \
    cd data_api && \
    git checkout 0.4.0-dev && \
    cd /kb/module && \
    cp -a data_api/lib/doekbase lib/ && \
    pip install -r /kb/module/data_api/requirements.txt


COPY ./ /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
