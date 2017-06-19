from flask import (Flask, request, current_app, render_template,
    send_from_directory)
from flask.ext.bootstrap import Bootstrap
from lipid_analysis import LipidAnalysis
from forms import LipidAnalysisForm
from config import Config
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

if __name__ == '__main__':
    app.run()

@app.route('/')
def hello():
    return 'hello'

@app.route('/lipid_analysis/', methods=['GET', 'POST'])
def lipid_analysis():
    config = Config()
    form_data = request.form
    form = LipidAnalysisForm()
    zip_path = None
    debug = 'debug' in request.args
    if form.validate_on_submit():
        root_path = Config.UPLOAD_FOLDER
        file1 = request.files[form.file1.name]
        file1.save(root_path + 'file1.txt')
        file2 = request.files[form.file2.name]
        file2.save(root_path + 'file2.txt')
        file1_path = root_path + 'file1.txt'
        file2_path = root_path + 'file2.txt'

        la = LipidAnalysis([file1_path, file2_path])
        la.filter_rows(form.data['retention_time_filter'],
                form.data['group_pq_filter'],
                form.data['group_sn_filter'],
                form.data['group_area_filter'],
                form.data['group_height_filter']
        )
        la.subtract_blank(form.data['blank'], form.data['mult_factor'])
        la.remove_columns(form.data['remove_cols'])
        la.normalize(form.data)
        subclass_stats, class_stats = la.calc_class_stats(form.data['class_stats'])
        zip_path = la.write_results()

    context = {'params': {}}
    if debug:
        context['params'] = {'debug': True}
    return render_template('lipid_analysis.html', form=form, zip_path=zip_path, **context)

@app.route('/file/<filename>')
def file(filename):
    file_dir = Config.UPLOAD_FOLDER
    return send_from_directory(file_dir, filename)


