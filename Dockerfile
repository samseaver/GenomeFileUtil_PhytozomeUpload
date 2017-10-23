FROM kbase/kbase:sdkbase.latest
MAINTAINER KBase Developer
# -----------------------------------------

# Make sure SSL certs are properly installed
RUN apt-get install python-dev libffi-dev libssl-dev \
    && pip install pyopenssl ndg-httpsclient pyasn1 \
    && pip install requests --upgrade \
    && pip install 'requests[security]' --upgrade

RUN mkdir -p /kb/module/work && mkdir -p /kb/module/lib
WORKDIR /kb/module

# Genbank uploader uses the Transform script_utils.py, so install that with
# dependencies to /kb/deployment/lib 
RUN git clone https://github.com/kbase/transform && \
    cd transform && \
    git checkout 3d0f47b && \
    cd /kb/module && \
    cp -a transform/lib/biokbase /kb/deployment/lib/ && \
    cp transform/plugins/scripts/upload/trns_transform_FASTA_DNA_Assembly_to_KBaseGenomeAnnotations_Assembly.py lib/. && \
    . /kb/module/transform/deps/pylib.sh

# Install the data_api dependencies.  The code is directly copied into this repo
# right now so we can make hotfixes
RUN git clone https://github.com/kbase/data_api && \
    cd /kb/module/data_api && \
    git checkout 0.4.1-dev && \
    cd /kb/module && \
    pip install thrift && \
    pip install -r /kb/module/data_api/requirements.txt && \
    rm -rf /kb/module/data_api


COPY ./ /kb/module
RUN chmod -R 777 /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
