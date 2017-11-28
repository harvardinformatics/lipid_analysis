from flask import (request, current_app, render_template,
    send_from_directory)
from lipidx.lipid_analysis import LipidAnalysis
from lipidx import forms
from lipidx import app
import logging
import sys, os
from bokeh.plotting import figure, output_file, show
from bokeh.models import HoverTool, Whisker, ColumnDataSource, Span, Range1d
from bokeh.embed import components
from bokeh.palettes import d3
from sklearn.decomposition import PCA
import numpy

@app.route('/lipid_analysis/', methods=['GET', 'POST'])
def lipid_analysis():
    form_data = request.form
    form = forms.LipidAnalysisForm()
    zip_path = None
    script = None
    div = None
    debug = 'debug' in request.args
    context = {'params': {}}
    if debug:
        context['params'] = {'debug': True}
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
        if form.data['class_stats']:
            subclass_stats, class_stats = la.calc_class_stats()
            context['class_script'], context['class_div'] = la.class_plot()
        context['volcano_script'], context['volcano_div'] = la.volcano_plot(form.data)
        zip_path = la.write_results()
    return render_template('lipid_analysis.html', form=form, zip_path=zip_path, **context)

@app.route('/volcano/', methods=['GET', 'POST'])
def volcano():
    form_data = request.form
    form = forms.VolcanoForm()
    zip_path = None
    script = None
    div = None
    context = {}
    if form.validate_on_submit():
        root_path = app.config['UPLOAD_FOLDER']
        file1 = request.files[form.file1.name]
        file1.save(root_path + 'file1.txt')
        file1_path = root_path + 'file1.txt'

        la = LipidAnalysis([file1_path])
        context['volcano_script'], context['volcano_div'] = la.volcano_plot(form.data)
        zip_path = la.write_results()
    return render_template('volcano.html', form=form, zip_path=zip_path, **context)

@app.route('/pca/', methods=['GET', 'POST'])
def pca():
    form_data = request.form
    form = forms.PCAForm()
    zip_path = None
    context = {}
    if form.validate_on_submit():
        root_path = app.config['UPLOAD_FOLDER']
        file1 = request.files[form.file1.name]
        file1.save(root_path + 'pca_file.txt')
        #path = '/Users/portermahoney/sites/lipidx_dev/lipidx/lipidx/files/compounds_mc.csv'
        path = root_path + 'pca_file.txt'
        samples = []
        names = []
        data = []
        with open(path,'r') as f:
            for ln in f:
                row = ln.split(',')
                if 'name' in ln:
                    samples = [x.replace('/n', '') for x in row[1:]]
                else:
                    names.append(row[0])
                    data.append([float(x.replace('/n', '')) for x in row[1:]])
        cnt = len(data[0])
        pca = PCA(n_components=2)
        pca.fit(data)
        print(pca)
        print(pca.components_)
        print(pca.explained_variance_ratio_)
        print(samples)
        p = figure(
                title = 'pca ploti',
                x_axis_label = ('pc1: %.2f%%' % (pca.explained_variance_ratio_[0] *
                    100)),
                y_axis_label = ('pc2: %.2f%%' % (pca.explained_variance_ratio_[1] *
                    100)),
                width = 1000, height = 800, toolbar_location = "above"
        )
        class_points = p.circle(pca.components_[0][0:4], pca.components_[1][0:4], size=10, color=d3['Category20'][20][1])
        #legend_items.append((class_name, [class_points]))
        class_points = p.circle(pca.components_[0][4:8], pca.components_[1][4:8], size=10, color=d3['Category20'][20][2])
        class_points = p.circle(pca.components_[0][8:], pca.components_[1][8:], size=10, color=d3['Category20'][20][3])
        context = {}
        context['pca_script'], context['pca_div'] = components(p)
    return render_template('pca.html', form=form, zip_path=zip_path, **context)

@app.route('/file/<filename>')
def file(filename):
    file_dir = app.config['UPLOAD_FOLDER']
    return send_from_directory(file_dir, filename)


