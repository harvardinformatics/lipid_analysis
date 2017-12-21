from flask import Flask
from flask_bootstrap import Bootstrap
import config
import logging
import sys, os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'lipid_analysis.log'))
logging.basicConfig(filename=os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
                                                            'lipid_analysis.log'), level=logging.DEBUG)
app = Flask(__name__)
app.config.from_object(config)
bootstrap = Bootstrap()
bootstrap.init_app(app)

import lipidx.views
