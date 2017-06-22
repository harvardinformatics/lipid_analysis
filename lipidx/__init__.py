from flask import Flask
from flask.ext.bootstrap import Bootstrap
from lipidx.config import Config
import logging
import sys, os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'lipid_analysis.log'))
logging.basicConfig(filename=os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
                                                            'lipid_analysis.log'), level=logging.DEBUG)
app = Flask(__name__)
conf = Config()
app.config.from_object(conf)
bootstrap = Bootstrap()
bootstrap.init_app(app)

import lipidx.views
