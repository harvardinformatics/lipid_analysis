from flask import Flask
from flask_bootstrap import Bootstrap
from lipidx.config import Config
import logging
import sys, os

class ReverseProxied(object):
    def __init__(self, app, script_name = None, scheme = None, server = None):
        self.app = app
        self.script_name = script_name
        self.server = server

    def __call__(self, environ, start_response):
        script_name = environ.get('HTTP_X_SCRIPT_NAME', '') or self.script_name
        if script_name:
            environ['SCRIPT_NAME'] = script_name
            path_info = environ['PATH_INFO']
            if path_info.startswith(script_name):
                environ['PATH_INFO'] = path_info[len(script_name):]
        scheme = environ.get('HTTP_X_SCHEME', '') or self.server
        if scheme:
            environ['wsgi.url_scheme'] = scheme
        server = environ.get('HTTP_X_FORWARDED_SERVER', '') or self.server
        if server:
            environ['HTTP_HOST'] = server
        return self.app(environ, start_response)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'lipid_analysis.log'))
logging.basicConfig(filename=os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
                                                            'lipid_analysis.log'), level=logging.DEBUG)
app = Flask(__name__)
app.wsgi_app = ReverseProxied(app.wsgi_app, script_name='/lipidx')
conf = Config()
app.config.from_object(conf)
bootstrap = Bootstrap()
bootstrap.init_app(app)

import lipidx.views
