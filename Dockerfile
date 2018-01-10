FROM continuumio/miniconda3
EXPOSE 5000
RUN conda install -y \
        flask \
        pillow \
        numpy \
        scipy \
        pandas \
        bokeh=0.12.6 \
        scikit-learn
RUN pip install \
        flask_bootstrap \
        flask_wtf
RUN conda install -c \
        conda-forge phantomjs
RUN conda install -c \
        conda-forge selenium
ADD . /code
WORKDIR /code
ENV PYTHONPATH /code:/code/lipidx
ENTRYPOINT ["python", "run.py"]
