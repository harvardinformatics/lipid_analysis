from flask import (request, current_app, render_template,
    send_from_directory)
from lipidx.lipid_analysis import LipidAnalysis
from lipidx.forms import LipidAnalysisForm
from lipidx import app
import logging
import sys, os
#import plotly
#from plotly.graph.objs import Scatter, Layout
from bokeh.plotting import figure, output_file, show
from bokeh.models import HoverTool
from bokeh.embed import components

@app.route('/')
def hello():
    x = [1,2,3,4,5]
    y = [6,7,2,4,5]
    output_file('test.html')
    TOOLS = "hover"
    p = figure(title='test', tools=TOOLS, x_axis_label = 'x', y_axis_label = 'y')
    p.circle(x, y, size=10, color="red", legend='Temp.', alpha=0.5)
    hover = p.select_one(HoverTool)
    hover.point_policy = "follow_mouse"
    hover.tooltips = [
            ("Name", "test")
    ]
    script, div = components(p)
    '''plotly.offline.plot({
        'data': [Scatter(x=[1,2,3,4], y=[4,3,2,1])],
        'layout': Layout(title='test')
    })'''
    return render_template('scatter.html', script = script, div = div)

@app.route('/lipid_analysis/', methods=['GET', 'POST'])
def lipid_analysis():
    form_data = request.form
    form = LipidAnalysisForm()
    zip_path = None
    debug = 'debug' in request.args
    if form.validate_on_submit():
        root_path = app.config['UPLOAD_FOLDER']
        file1 = request.files[form.file1.name]
        file1.save(root_path + 'file1.txt')
        file2 = request.files[form.file2.name]
        file2.save(root_path + 'file2.txt')
        file1_path = root_path + 'file1.txt'
        file2_path = root_path + 'file2.txt'

        la = LipidAnalysis([file1_path, file2_path], debug)
        la.remove_rejects()
        la.group_ions(form.data['group_ions_within'])
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
    file_dir = app.config['UPLOAD_FOLDER']
    return send_from_directory(file_dir, filename)


