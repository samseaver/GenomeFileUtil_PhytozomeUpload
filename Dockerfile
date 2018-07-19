FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------

RUN pip install biopython==1.70
RUN pip install mock

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
