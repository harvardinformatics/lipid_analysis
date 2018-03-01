# Lipidx - A lipid analysis tool for small molecule mass spec?

## Development
Lipidx is a Flask application with a number of graphics dependencies including the Anaconda bokeh library

    $ conda create -n lipidx bokeh=0.12.6 scikit-learn flask pillow pandas
    $ conda install -c conda-forge phantomjs selenium
    $ pip install flask_bootstrap flask_wtf
    $ git clone https://github.com/harvardinformatics/lipidx.git
    $ cd lipidx
    $ export PYTHONPATH=`pwd`
    $ export FLASK_APP=lipidx
    $ export FLASK_DEBUG=1
    $ flask run

You should be able to see the application at http://localhost:5000/lipidx/lipid_analysis/

## Installation
Lipidx can be deployed as a Docker container.

    $ git clone https://github.com/harvardinformatics/lipidx.git
    $ cd lipidx
    $ docker build -t lipidx .
    $ docker run -d -p 8000:80 --name lipidx lipidx

Once running, navigate to http://localhost:8000/lipid_analysis

For production, you should set the following environment variables:

    $ export LIPIDX_ADMIN_EMAILS='one@harvard.edu,two@harvard.edu'
    $ export LIPIDX_KEY='someunintelligiblestring'
    $ export WTF_CSRF_KEY='someotherunintelligiblestring'

