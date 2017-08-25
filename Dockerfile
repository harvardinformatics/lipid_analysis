FROM conda/miniconda3-centos7
RUN conda install -y \
        flask \
        pillow \
        numpy \
        scipy \
        pandas \
        bokeh
RUN pip install \
        flask_bootstrap \
        flask_wtf
RUN conda install -c \
        conda-forge phantomjs
RUN conda install -c \
        conda-forge selenium
ADD . /code
WORKDIR /code
RUN conda create --name lipidx
ENV PYTHONPATH /code:/code/lipidx
RUN source activate lipidx
ENTRYPOINT ["python", "run.py"]
