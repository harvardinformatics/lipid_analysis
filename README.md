# Lipidx - A lipid analysis tool for small molecule mass spec?

## Installation
Lipidx can be deployed as a Docker container.

    $ git clone https://github.com/harvardinformatics/lipidx.git
    $ cd lipidx
    $ docker build -t lipidx .
    $ docker run -d -p 5000:5000 --name lipidx lipidx

Once running, navigate to http://localhost:5000/lipid_analysis 

For production, you should set the following environment variables:

    $ export LIPIDX_ADMIN_EMAILS='one@harvard.edu,two@harvard.edu'
    $ export LIPIDX_KEY='someunintelligiblestring'
    $ export WTF_CSRF_KEY='someotherunintelligiblestring'

