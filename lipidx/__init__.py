from flask import Flask
from flask_bootstrap import Bootstrap
from lipidx.views import lipidx_bp
from lipidx import config
import logging
import sys, os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'lipid_analysis.log'))
logging.basicConfig(filename=os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
                                                            'lipid_analysis.log'), level=logging.DEBUG)
app = Flask(__name__)
app.register_blueprint(lipidx_bp, url_prefix='/lipidx')
app.config.from_object(config)
bootstrap = Bootstrap()
bootstrap.init_app(app)

import lipidx.views
