from flask import (request, current_app, render_template,
    send_from_directory)
from lipidx.lipid_analysis import LipidAnalysis
from lipidx import forms
from lipidx import app
import logging
import sys, os
from bokeh.plotting import figure, output_file, show
from bokeh.models import (HoverTool, Whisker, ColumnDataSource, Span, Range1d,
        Legend)
from bokeh.embed import components
from bokeh.palettes import d3
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
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

@app.route('/pca_test/', methods=['GET', 'POST'])
def pca_test():
    context = {}
    form_data = request.form
    form = forms.PCAForm()
    zip_path = None
    rng = numpy.random.RandomState(1)
    data = numpy.dot(rng.rand(2,2), rng.randn(2,20)).T
    p = figure(
            title = 'input test',
            x_axis_label = ('y'),
            y_axis_label = ('x'),
            width = 800, height = 800, toolbar_location = "above"
    )
    p.circle(data[:, 0], data[:, 1], size=10)
    context['input_script'], context['input_div'] = components(p)
    pca = PCA(n_components=2)
    data_pca = pca.fit_transform(data)
    p = figure(
            title = 'pca test',
            x_axis_label = ('pc1: %.2f%%' % (pca.explained_variance_ratio_[0] *
                100)),
            y_axis_label = ('pc2: %.2f%%' % (pca.explained_variance_ratio_[1] *
                100)),
            width = 800, height = 800, toolbar_location = "above"
    )
    p.circle(data_pca[:, 0], data_pca[:, 1], size=10)
    context['pca_script'], context['pca_div'] = components(p)
    return render_template('pca.html', form=form, zip_path=zip_path, **context)

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
        sample_grps = {}
        names = []
        data = []
        with open(path,'r') as f:
            for ln in f:
                row = ln.split(',')
                if 'name' in ln:
                    samples = [x.replace('/n', '') for x in row[1:]]
                    for i, s in enumerate(samples):
                        prefix = s.split('_')[0]
                        if prefix not in sample_grps:
                            sample_grps[prefix] = []
                        sample_grps[prefix].append(i)
                else:
                    names.append(row[0])
                    data.append([float(x.replace('/n', '')) for x in row[1:]])
        cnt = len(data[0])
        print(data)
        #data = numpy.matrix.transpose(numpy.array(data))
        p = figure(
                title = 'pca input',
                x_axis_label = ('y'),
                y_axis_label = ('x'),
                width = 800, height = 800, toolbar_location = "above"
        )
        #p.circle(data[:, 0], data[:, 1], size=10)
        #context['input_script'], context['input_div'] = components(p)
        pca = PCA(n_components=2)
        #pca.fit(data)
        # what is scale - changes numbers to be close to one or 2 not extremes,
        # with std
        # what is centering -
        # try with R for comparison
        std_data = StandardScaler().fit_transform(data)
        data_pca = pca.fit_transform(data)
        #data_pca = pca.inverse_transform(data_pca)
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
                width = 800, height = 800, toolbar_location = "above"
        )
        legend_items = []
        palette_key = 1
        #p.circle(pca.components_[0], pca.components_[1], size=10)
        for grp, indexes in sample_grps.items():
            pc_1 = [x for i, x in enumerate(pca.components_[0]) if i in indexes]
            pc_2 = [y for i, y in enumerate(pca.components_[1]) if i in indexes]
            class_points = p.circle(pc_1, pc_2, size=10, color=d3['Category20'][20][palette_key])
            palette_key += 1
            legend_items.append((grp, [class_points]))
        legend = Legend(
                items = legend_items,
                click_policy = 'hide',
                location = (0, -30)
        )
        p.add_layout(legend, 'right')
        context = {}
        context['pca_script'], context['pca_div'] = components(p)
    return render_template('pca.html', form=form, zip_path=zip_path, **context)

@app.route('/file/<filename>')
def file(filename):
    file_dir = app.config['UPLOAD_FOLDER']
    return send_from_directory(file_dir, filename)


