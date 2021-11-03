FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# RUN apt-get update

RUN conda install -c conda-forge phantomjs

RUN pip install --upgrade pip && \
    pip install bokeh && \
    pip install selenium && \
    pip install pandas && \
    pip install scipy

# Loading PlantSEED data
RUN echo "Current PlantSEED commit: 5520f62ed57e95589bc9e5e30871a0adfa27bff3"
RUN git clone -b kbase_release https://github.com/ModelSEED/PlantSEED /kb/module/PlantSEED

# Loading Biochemistry data
RUN git clone -b plant_kbase_release https://github.com/ModelSEED/ModelSEEDDatabase /kb/module/ModelSEEDDatabase

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

COPY ./data/compartments.txt /kb/module/data/
WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
