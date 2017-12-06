import csv
import os
import numpy
import re
from math import pi, isnan
from scipy.stats import ttest_ind
import zipfile
import re
from flask import request
from flask import current_app as app
from collections import OrderedDict
from lipidx.forms import LipidAnalysisForm
import logging
from bokeh.layouts import gridplot
from bokeh.plotting import figure, show
from bokeh.models import (HoverTool, ColumnDataSource, Whisker, DataRange1d,
    Range1d, BoxAnnotation, Legend)
from bokeh.embed import components
from bokeh.sampledata.autompg import autompg as df
from bokeh.palettes import d3
from bokeh.io import export_svgs, export_png

class LipidAnalysis:
    MAX_GROUPS = 10
    MAX_CLASSES = 20
    ROUND_TO = 2
    POST_NORMAL_ROUND = 8
    NEGATIVE_IONS_WITH_PLUS = ['HCOO', 'CH3COO', 'CL']
    # limited cols as a string because remove cols takes a string (for now)
    LIMITED_COLS = '''name, ret_time, lipidion, class, fattyacid, fa1, fa2, fa3,
        fa4, calcmz, ionformula, area, ratio, p_value'''
    MAX_VOLCANO_PLOTS = 3

    def __init__ (self, paths, debug = False):
        self.paths = paths
        # debug adds cols to results to show pre normalized values
        self.debug = debug
        self.area_start = 'area['
        self.group_area_start = 'grouparea['
        self.groups = {}

        # cols should be the same in all files
        # col names will be taken from first file
        self.rows = self.get_rows_from_files(self.paths)
        # several functions will need to know the group names
        self.groups = self.get_groups()

        # filled in calc_class_stats
        self.class_stats = {}
        self.subclass_stats = {}
        self.class_dict = {}
        self.subclass_dict = {}

        # file paths, eventualy these may not be hardcoded
        self.root_path = app.config['UPLOAD_FOLDER']
        lipid_class_file = 'lipidKey.csv'
        self.lipid_class_path = app.config['BASE_DIR'] + '/' + lipid_class_file
        self.lipid_results_file = 'lipid_analysis.csv'
        self.lipid_results_path = self.root_path + self.lipid_results_file
        self.lipid_results_limited_file = 'lipid_analysis_summary.csv'
        self.lipid_results_limited_path = self.root_path + self.lipid_results_limited_file
        self.paths_to_zip = {}
        zip_file = 'lipid_results.zip'
        self.zip_path = self.root_path + zip_file

        # set lipid classes
        self.class_keys = self.load_lipid_classes()

    def get_rows_from_files(self, paths):
        rows = OrderedDict()
        cols = [] # cols common to all files
        set_name = True
        for path in paths:
            if path:
                with open(path,'r') as f:
                    for i,ln in enumerate(f):
                        if i == 0:
                            row_cols = []
                        if (ln.startswith('#') or ln.startswith('\t') or
                                ln.startswith('\n')):
                            continue
                        if not row_cols:
                            ln = ln.replace('\n', '')
                            file_cols = re.split('\t|,', ln)
                            # Only add name and ret_time if not there
                            # for volcano only it will already be there
                            if 'name' in file_cols:
                                set_name = False
                                row_cols = file_cols
                            else:
                                row_cols = ['name', 'ret_time'] # these two columns added first
                                row_cols.extend(file_cols)
                        else: # data lines
                            ln = ln.replace('\n', '')
                            file_row = re.split('\t|,', ln)
                            if set_name:
                                row = ['', ''] # filler vals for name and ret_time
                                row.extend(file_row)
                            else:
                                row = file_row
                            # remove trailing newline
                            row[(len(row) - 1)] = row[(len(row) - 1)].strip('\n')
                            # lowercase the cols before creating dictionary
                            row_cols = [x.lower() for x in row_cols]
                            row_d = OrderedDict(zip(row_cols, row))
                            if set_name: # not necessary if volcano only
                                # calc retention time: average of GroupTopPos
                                ret_time = round(numpy.mean(self.list_col_type(row_d, 'grouptoppos')), self.ROUND_TO)
                                row_d['ret_time'] = ret_time # add to row
                                # unique name for row LipidIon + ret_time
                                name = row_d['lipidion'] + '_' + str(ret_time)
                                row_d['name'] = name
                            else:
                                name = row_d['name']
                            if name in rows: # rare case
                                # if lipid has same name then keep the one with
                                # greater area
                                avg_areas = self.list_col_type(row_d, self.area_start)
                                avg_prev_areas = self.list_col_type(rows[name], self.area_start)
                                if avg_areas > avg_prev_areas:
                                    rows[name] = row_d
                            else: # add the new row
                                rows[name] = row_d

                # ensure cols are the same for all rows
                if not cols: # set cols to first file
                    cols = row_cols
                else: # after second file ensure cols in each row are the same
                    cols = set(cols).intersection(row_cols)
                    for key, row in rows.items():
                        rows[key] = self.limit_row_cols(cols, row)
        return rows

    def limit_row_cols(self, cols, row):
        if cols:
            extra = set(row.keys()) - set(cols)
            for k in extra:
                del row[k]
        return row

    def get_cols(self, start = None):
        first = list(self.rows.keys())[0]
        keys = list(self.rows[first].keys())
        if start:
            keys = [i for i in keys if i.startswith(start)]
        return keys

    def write_results(self):
        # get a list of results sorted by key
        res = [x for y, x in sorted(self.rows.items(), key=lambda t: t[0].lower())]
        # get cols from first row
        cols = list(res[0].keys())
        self.write_csv(self.lipid_results_path, cols, res)
        # save results with limited cols
        self.remove_columns(self.LIMITED_COLS, whitelist = True)
        res = [x for y, x in sorted(self.rows.items(), key=lambda t: t[0].lower())]
        cols = list(res[0].keys())
        self.write_csv(self.lipid_results_limited_path, cols, res)
        # create a zip file for lipids and stats
        z = zipfile.ZipFile(self.zip_path, "w")
        z.write(self.lipid_results_path, self.lipid_results_file)
        z.write(self.lipid_results_limited_path, self.lipid_results_limited_file)
        # save any other files (volcano, class_summary)
        for filename, path in self.paths_to_zip.items():
            if os.path.exists(path):
                z.write(path, filename)
        z.close
        return self.zip_path

    def write_csv(self, path, cols, rows):
        success = False
        if self.rows:
            if not os.path.exists(self.root_path):
                os.makedirs(self.root_path)
            with open(path,'w') as c:
                w = csv.DictWriter(c, cols)
                w.writeheader()
                w.writerows(rows)
                success = True
        return success

    def subtract_blank(self, blank, mult_factor):
        if blank and self.rows:
            subtracted = {}
            area_cols = self.get_cols(self.area_start)
            blank_start = self.area_start + blank
            blank_cols = self.get_cols(blank_start)
            for name, row in self.rows.items():
                # calculate avg blank
                avg_blank = self.calculate_avg_blank(blank_cols, row)
                include_row = False
                # subtract blank * mult_factor from all area cols
                for col in area_cols:
                    sub = round((float(row[col]) - (avg_blank * mult_factor)),
                            self.ROUND_TO)
                    # negative areas are not permited, neg becomes 0
                    if sub < 0:
                        sub = 0
                    row[col] = sub
                    # only include row if atleast one non blank area is non-zero
                    if sub > 0 and col not in blank_cols:
                        include_row = True
                if include_row:
                    row['avg_blank'] = avg_blank
                    subtracted[name] = row
                    for col in blank_cols:
                        del(row[col])
            self.rows = subtracted
            self.groups = self.get_groups()

    def calculate_avg_blank(self, blank_cols, row):
            blank_vals = []
            for col in blank_cols:
                if col in row:
                    blank_vals.append(float(row[col]))
            return round(numpy.mean(blank_vals), self.ROUND_TO)

    def remove_columns(self, remove_cols, whitelist = False):
        if self.rows:
            # ensure user entered cols are lowercase and without spaces
            remove_cols = [x.strip().lower() for x in remove_cols.split(',')]
            clean_selected = {}
            removed_cols = []
            for col in self.get_cols():
                prefix = col.split('[')[0]
                # remove cols contains the start of the column name before [
                if prefix.lower() in remove_cols:
                    removed_cols.append(col)
            for name, row in self.rows.items():
                new_row = OrderedDict()
                if whitelist: # include only cols from input
                    for col, val in row.items():
                        if col in removed_cols:
                            new_row[col] = val
                else: # remove cols from input
                    for col, val in row.items():
                        if col not in removed_cols:
                            new_row[col] = val
                clean_selected[name] = new_row
            self.rows = clean_selected

    def list_col_type(self, row, col_type):
        # get a list of all cols that startwith col_type
        lst = []
        for name, val in row.items():
            if name.startswith(col_type):
                lst.append(float(val))
        return lst

    def group_ions(self, within):
        # for lipid charges with different ions but ret time within 0.9 only the
        # lipid ion with the greatest area will be kept
        lc_grps = {}
        ion_dups = {}
        # group rows by lipid charge and ret time
        for name, row in self.rows.items():
            # capture lipid_charge
            grps = re.search('(.*[+,-])(.*)_.*', name)
            lipid_charge = grps.group(1)
            adduct = grps.group(2)
            # some negative ions actually start with a +, change those to - for
            # the purpose of grouping
            if adduct in self.NEGATIVE_IONS_WITH_PLUS:
                lipid_charge = lipid_charge.replace(')+', ')-')
            ret_time = row['ret_time']
            if lipid_charge not in lc_grps:
                # for ret_time store a list of lipid rows that are within 0.9
                lc_grps[lipid_charge] = {ret_time: [name]}
            else:
                found = False
                # see if curr ret time is within 0.9 of any others in this group
                # if so choose the closest
                min_diff = None
                min_prev_ret = None
                for prev_ret in lc_grps[lipid_charge]:
                    diff = abs(float(ret_time) - float(prev_ret))
                    if diff < within:
                        if not min_diff or min_diff > diff:
                            min_diff = diff
                            min_prev_ret = prev_ret
                            found = True
                if found:
                    # add to closest prev_ret group
                    lc_grps[lipid_charge][min_prev_ret].append(name)
                    # keep list of charges and ret times that have dups, for
                    # faster filtering
                    if lipid_charge not in ion_dups:
                        ion_dups[lipid_charge] = []
                    if min_prev_ret not in ion_dups[lipid_charge]:
                        ion_dups[lipid_charge].append(min_prev_ret)
                else:
                    # start new ret_time bucket
                    lc_grps[lipid_charge][ret_time] = [name]

        for lc, r_times in ion_dups.items():
            # loop through dups and keep only the one with largest area
            for r in r_times:
                max_area = 0
                keep = None
                for k in lc_grps[lc][r]:
                    areas = self.list_col_type(self.rows[k], self.area_start)
                    avg_area = numpy.mean(areas)
                    if avg_area > max_area:
                        max_area = avg_area
                        keep = k
                for j in lc_grps[lc][r]:
                    # delete non-max areas of dup ions
                    if j != keep:
                        del(self.rows[j])

    def remove_rejects(self):
        # rejects must be removed before grouping ions
        selected = {}
        for name, row in self.rows.items():
            if row['rej.'] == '0':
                selected[name] = row
        self.rows = selected

    def filter_rows(self, ret_time_fil, group_pq_fil, group_sn_fil, group_area_fil,
            group_height_fil):
        selected = {}
        for name, row in self.rows.items():
            # select only rows that pass all filters
            if self.filter_in(row, ret_time_fil, group_pq_fil, group_sn_fil,
                    group_area_fil, group_height_fil):
                selected[name] = row
        self.rows = selected

    def filter_in(self, row, ret_time_fil, group_pq_fil, group_sn_fil,
            group_area_fil, group_height_fil):
        if float(row['ret_time']) <= ret_time_fil:
            return False
        group_pq_max = max(self.list_col_type(row, 'grouppq'))
        if group_pq_max <= group_pq_fil:
            return False
        group_sn_max = max(self.list_col_type(row, 'groups/n'))
        if group_sn_max <= group_sn_fil:
            return False
        if group_area_fil > 0: # all values are pos ints, so skip if 0
            group_area_max = max(self.list_col_type(row, 'grouparea'))
            if group_area_max <= group_area_fil:
                return False
        if group_height_fil > 0: # all values are pos ints, so skip if 0
            group_height_max = max(self.list_col_type(row, 'groupheight'))
            if group_height_max <= group_height_fil:
                return False
        return True

    def normalize(self, form_data):
        # most common is to not normalize
        # or one can use avg intensity calculated from data
        # or input manual values
        if self.rows:
            normal = self.rows
            #TODO: comma to seperate sample values in input S1-1 S1-2 etc
            if form_data['normalize'] != 'none':
                area_cols = self.get_cols(self.area_start)
                # use manual values
                if form_data['normalize'] == 'values':
                    for name, row in normal.items():
                        for col in area_cols:
                            group, num = self.get_group_from_col(col)
                            form_name = 'normal_' + group
                            # TODO: what if they don't fill it out
                            if form_data[form_name]:
                                # form field is comma separated values for each
                                # numeric sample
                                val = form_data[form_name].split(',')
                                num = int(num)
                                if num < len(val):
                                    val = val[num].strip()
                                    normal[name][col] = round(float(row[col]) / float(val), self.POST_NORMAL_ROUND)
                # use calculated intensity
                elif form_data['normalize'] == 'intensity':
                    intensities = self.calc_intensities(area_cols)
                    for name, row in normal.items():
                        for col in area_cols:
                            sam = self.get_sample_from_col(col)
                            if intensities[sam] > 0:
                                normal[name][col] = round(float(row[col])/intensities[sam],
                                self.POST_NORMAL_ROUND)
                normal = self.recalc_avg(normal)
            self.rows = normal

    def calc_intensities(self, area_cols):
        # calc average from area_cols (intensity)
        intensities = {}
        for name, row in self.rows.items():
            for col in area_cols:
                sam = self.get_sample_from_col(col)
                if sam not in intensities:
                    intensities[sam] = []
                intensities[sam].append(float(row[col]))
        for sam, sum_lst in intensities.items():
            intensities[sam] = round(numpy.mean(intensities[sam]), self.POST_NORMAL_ROUND)
        return intensities

    def get_groups(self):
        groups = OrderedDict()
        area_cols = self.get_cols(self.area_start)
        for a_col in area_cols:
            group, num = self.get_group_from_col(a_col)
            if group not in groups:
                groups[group] = []
            groups[group].append(num)
        return groups

    def get_group_from_col(self, col):
        gr = col.split('[')[1]
        gr = gr.split('-')
        num = gr[1].split(']')[0]
        return gr[0], num

    def get_sample_from_col(self, col):
        sam = col.split('[')[1]
        sam = sam.split(']')
        return sam[0]

    def recalc_avg(self, normal):
        for name, row in normal.items():
            stats = OrderedDict()
            # for each group recalc the avg and std from areas
            for group, nums in self.groups.items():
                if group not in stats: # group like c, s1, s2
                    stats[group] = []
                for num in nums: # num replicates per group
                    num_col = self.area_start + group + '-' + num + ']'
                    stats[group].append(float(row[num_col]))
            for group, val_lst in stats.items():
                normal[name]['grouparea[' + group + ']'] = round(numpy.mean(val_lst), self.POST_NORMAL_ROUND)
                normal[name]['arearsd[' + group + ']'] = round(numpy.std(val_lst), self.POST_NORMAL_ROUND)
        return normal

    def calc_class_stats(self):
        # set to false if file not saved
        sub_success = True
        class_success = True
        class_stats = {}
        subclass_stats = {}
        for name, row in self.rows.items():
            # take subclass key from row
            subclass_key = row['class']
            if subclass_key in self.class_keys:
                # get corresponding names from class_keys
                subclass_name = self.class_keys[subclass_key]['subclass']
                class_name = self.class_keys[subclass_key]['class']
                # populate new subclasses and classes
                if subclass_name not in subclass_stats:
                    subclass_stats[subclass_name] = OrderedDict()
                    if class_name not in class_stats:
                        class_stats[class_name] = OrderedDict()
                # add row to group areas for the class and subclass that
                # correspond to the row
                subclass_stats[subclass_name] = self.group_areas(row,
                        subclass_stats[subclass_name])
                class_stats[class_name] = self.group_areas(row,
                        class_stats[class_name])
        # write files
        class_file = 'class_stats.csv'
        class_path = self.root_path + class_file
        subclass_file = 'subclass_stats.csv'
        subclass_path = self.root_path + subclass_file
        subclass_cols = self.stats_cols('subclass')
        self.subclass_stats, self.subclass_dict = self.compute_stats('subclass', subclass_stats)
        sub_success = self.write_csv(subclass_path, subclass_cols, self.subclass_dict.values())
        self.paths_to_zip[subclass_file] = subclass_path
        class_cols = self.stats_cols('class')
        self.class_stats, self.class_dict = self.compute_stats('class', class_stats)
        class_success = self.write_csv(class_path, class_cols, self.class_dict.values())
        self.paths_to_zip[class_file] = class_path
        return sub_success, class_success

    def stats_cols(self, cat):
        cols = [cat]
        for g in self.groups:
            cols.append(g + ' cnt')
            cols.append(g + ' avg')
            cols.append(g + ' std')
        return cols

    def compute_stats(self, cat, stats):
        rows = {}
        for name, groups in stats.items():
            row = {}
            # fill the class type with the name
            row[cat] = name
            # for each group add stats based on grp area lists
            for group, info in groups.items():
                row[group + ' cnt'] = info['cnt']
                avg = numpy.mean(info['grp_areas'])
                gr_sum = numpy.sum(info['grp_areas'])
                row[group + ' avg'] = avg
                stats[name][group]['sum'] = gr_sum
                stats[name][group]['avg'] = avg
                std = numpy.std(info['grp_areas'])
                #log_std = numpy.std(info['log_grp_areas'])
                row[group + ' std'] = std
                stats[name][group]['std'] = std
                #stats[name][group]['log_std'] = log_std
            rows[name] = row
        return stats, rows

    def group_areas(self, row, grp_info):
        # add stats for each group from the row
        for key in self.groups.keys():
            if key not in grp_info:
                # keep cnt and a list of group areas for group
                grp_info[key] = {'cnt': 0, 'grp_areas': [],
                        #'log_grp_areas': []
                }
            areas = self.list_col_type(row, self.area_start + key)
            '''for a in areas:
                # use 0.0 as the log if the area is 0.0 this will show there's
                # nothing in the subclass/group
                log = 0.0
                if a > 0.0:
                    log = numpy.log10(a)
                log_areas.append(log)'''
            if max(areas) > 0.0:
                grp_info[key]['cnt'] += 1
            # append one area and log area per row group because there will be
            # multiple rows in a subclass
            # later another average is taken across all rows in the subclass
            grp_info[key]['grp_areas'].append(numpy.mean(areas))
            #grp_info[key]['log_grp_areas'].append(numpy.mean(log_areas))
        return grp_info

    def load_lipid_classes(self):
        classes = {}
        with open(self.lipid_class_path,'r') as f:
            for i,ln in enumerate(f):
                ln = ln.replace('\n', '')
                ln = ln.rstrip(',')
                row = ln.split(',')
                if i == 0:
                    cols = [x.lower() for x in row]
                else:
                    row_dict = dict(zip(cols, row))
                    classes[row_dict['key']] = row_dict
        return classes

    def class_plot(self):
        # reorganize class_stats data for plotting
        # TODO: can we avoid the need for this regoranization?
        # TODO: remove log stats once we're happy with the log access
        data = {
                'lipid': [],
                'cnt': [],
                'sum': [],
                #'log_sum': [],
                'std': [],
                #'log_std': [],
                'x': [],
                'lower': [],
                'upper': [],
                'relative': [],
                #'log_relative': []
        }
        x = 1
        # TODO: FactorRange for groups lipid axis
        gr_data = OrderedDict()
        group_list = []
        for lipid_class, groups in self.class_stats.items():
            for group, stats in groups.items():
                if group not in group_list:
                    group_list.append(group)
                # get all group level data into a list
                if group not in gr_data:
                    gr_data[group] = {
                            'cnt': [],
                            'sum': [],
                            #'log_sum': [],
                            'x': [],
                            'relative':[],
                            #'log_relative':[]
                    }
                gr_data[group]['cnt'].append(stats['cnt'])
                gr_data[group]['sum'].append(stats['sum'])
                # don't try to take log of 0
                '''log = 0.0
                if stats['sum'] > 1.0:
                    log = numpy.log10(stats['sum'])
                gr_data[group]['log_sum'].append(log)'''
                gr_data[group]['x'].append(x)

                # lists of the data not organized by group used in plotting
                data['lipid'].append(lipid_class)
                data['cnt'].append(stats['cnt'])
                data['sum'].append(stats['sum'])
                #data['log_sum'].append(log)
                data['std'].append(stats['std'])
                #data['log_std'].append(stats['log_std'])
                data['lower'].append(stats['sum'] - stats['std'])
                data['upper'].append(stats['sum'] + stats['std'])
                data['x'].append(x)
                x += 1
            # put space between lipid groups
            data['lipid'].append(lipid_class)
            x += 1

        # get relative sums
        for lipid, groups in self.class_stats.items():
            for group, stats in groups.items():
                gr_sum = numpy.sum(gr_data[group]['sum'])
                # replace with a small value rather than divide by zero
                if gr_sum <= 0.0:
                    gr_sum = 0.1
                relative = stats['sum'] / gr_sum
                relative_percent = relative * 100
                gr_data[group]['relative'].append(relative_percent)
                data['relative'].append(relative_percent)
                # prevent taking the log of 0
                '''log_relative = 0.0
                if relative > 0.0:
                    log_relative = numpy.log10(relative)
                gr_data[group]['log_relative'].append(log_relative)
                data['log_relative'].append(log_relative)
                '''
        bar_cnt = self.bar_chart(gr_data, data, 'x', 'cnt', 'Number of Lipids', 'nb of lipids', data['lipid'], data['cnt'])
        bar_sum = self.bar_chart(gr_data, data, 'x', 'sum', 'Area', 'sum of area per group', data['lipid'], data['sum'], 'std')
        bar_log_sum = self.bar_chart(gr_data, data, 'x', 'sum', 'Area, log', 'sum of area per group', data['lipid'], data['sum'], 'std', 'log')
        bar_relative = self.bar_chart(gr_data, data, 'x', 'relative', 'Relative area', 'sum of area per group / total sum of area %', data['lipid'], data['relative'])
        bar_log_relative = self.bar_chart(gr_data, data, 'x', 'relative', 'Relative area, log', 'sum of area per group / total sum of area', data['lipid'], data['relative'], None, 'log')

        bars = gridplot(
                [bar_cnt],
                [bar_sum],
                [bar_log_sum],
                [bar_relative],
                [bar_log_relative]
        )
        script, div = components(bars)
        return script, div

    def bar_chart(self, gr_data, data, x, y, title, y_label, x_range, y_range, std = None, y_axis_type = None):
        if y_axis_type:
            bottom = 0.0000000001
        else:
            bottom = 0
            y_axis_type = 'auto'

        bar = figure(title=title, y_axis_label = y_label,
                x_range = x_range, width = 1500, y_axis_type
                = y_axis_type)
        bar.y_range.start = bottom
        bar.xaxis.major_label_orientation = pi/4

        if std:
            base = data[x]
            upper = [u + data[std][i] for i, u in enumerate(data[y])]
            lower = [u - data[std][i] for i, u in enumerate(data[y])]
            errors = ColumnDataSource(data=dict(base=base, lower=lower, upper=upper))
            bar.add_layout(Whisker(source=errors, base='base', upper='upper',
                lower='lower', level='annotation'))
        palette_key = 0
        for group, d in gr_data.items():
            bar.vbar(x = d[x], width = 0.5, top = d[y], bottom = bottom, legend
                    = group, color = d3['Category10'][self.MAX_GROUPS][palette_key])
            palette_key += 1
        # TODO: clear all somehow and select all
        bar.legend.click_policy = 'hide'
        bar.output_backend = 'svg'
        # save the png for the zip file
        chart_file_png = 'class_chart_' + title.replace(' ', '_').replace(',',
        '').lower() + '.png'
        chart_path_png = self.root_path + chart_file_png
        export_png(bar, filename=chart_path_png)
        # save paths to zip
        self.paths_to_zip[chart_file_png] = chart_path_png
        return bar

    def calc_ratio(self, group1, group2):
        ratio_name = group1 + '-over-' + group2
        ratio_col_name = 'ratio[' + ratio_name + ']'
        for key, row in self.rows.items():
            # only calculate ratio if it's not already there
            # it will already be there for ppl using volcano endpoint
            if ratio_col_name in row:
                return ratio_name
            dividend = float(row[self.group_area_start + group1 + ']'])
            divisor = float(row[self.group_area_start + group2 + ']'])
            # set ratio artificaly to 0.1 or 10 if zero
            if dividend == 0.0:
                ratio = float(0.1)
            elif divisor == 0.0:
                ratio = float(10.0)
            else:
                ratio = dividend/divisor
            self.rows[key][ratio_col_name] = ratio
            self.rows[key]['log_ratio[' + ratio_name + ']'] = numpy.log2(ratio)
            dividend_list = self.list_col_type(row, self.area_start + group1)
            divisor_list = self.list_col_type(row, self.area_start + group2)
            # if both lists are zero there is no significance so set to 1
            if (numpy.sum(divisor_list) + numpy.sum(dividend_list)) == 0.0:
                p = 1
            else: # perform ttest to get p_value
                t, p = ttest_ind(dividend_list, divisor_list, equal_var = False)
            self.rows[key]['p_value[' + ratio_name + ']'] = p
            self.rows[key]['log_p_value[' + ratio_name + ']'] = numpy.log10(p) * -1
        return ratio_name

    def get_plots(self, form_data):
        plots = []
        prefix = 'group'
        for g in range(1, (self.MAX_VOLCANO_PLOTS * 2), 2):
            group1 = prefix + str(g)
            group2 = prefix + str(g + 1)
            if (form_data[group1] and form_data[group2] and
                    form_data[group1].lower() in
            self.groups and form_data[group2].lower() in self.groups):
                plots.append((form_data[group1].lower(),
                    form_data[group2].lower()))
        return plots

    def volcano_plot(self, form_data):
        plots = self.get_plots(form_data)
        # set params for highlight of up and down reg regions
        ratio_highlight = form_data['ratio_highlight']
        pvalue_highlight = form_data['pvalue_highlight']
        if ratio_highlight <= 0:
            ratio_highlight = LipidAnalysisForm.RATIO_HIGHLIGHT_DEFAULT
        if pvalue_highlight <= 0:
            pvalue_highlight = LipidAnalysisForm.PVALUE_HIGHLIGHT_DEFAULT
        ratio_highlight = numpy.log2(ratio_highlight)
        pvalue_highlight = (numpy.log10(pvalue_highlight) * -1.0)
        plot_list = []
        script = None
        div = None
        if plots and not self.class_keys:
            self.class_keys = self.load_lipid_classes()
        y_range = []
        for (group1, group2) in plots:
            ratio_name = self.calc_ratio(group1, group2)
            data = {}
            for key, row in self.rows.items():
                subclass_key = row['class']
                class_name = self.class_keys[subclass_key]['class']
                if class_name not in data:
                    data[class_name] = {
                            'lipid': [],
                            'log2': [],
                            'p': []
                    }
                # TODO: fix the lipid being used as labels
                data[class_name]['lipid'].append(key)
                data[class_name]['log2'].append(row['log_ratio[' + ratio_name +
                    ']'])
                data[class_name]['p'].append(row['log_p_value[' + ratio_name +
                ']'])
                y_range.append(row['log_p_value[' + ratio_name + ']'])
            p = figure(title = ratio_name, x_axis_label = 'log2(ratio)', y_axis_label = '-log10(p value)', width = 1000, height = 800, toolbar_location = "above")
            hover = HoverTool(tooltips=[
                ('name', "@lipid"),
                ('ratio', "@log2"),
                ('p value', "@p")
            ])
            p.add_tools(hover)
            # add highlight to interesting regions
            box_lt = BoxAnnotation(plot=p, bottom=pvalue_highlight, right=(ratio_highlight * -1), fill_alpha=0.1, fill_color='red')
            box_rt = BoxAnnotation(plot=p, bottom=pvalue_highlight, left=ratio_highlight, fill_alpha=0.1, fill_color='green')
            p.renderers.extend([box_lt, box_rt])

            palette_key = 0
            legend_items = []
            for class_name, source in data.items():
                class_points = p.circle('log2', 'p', size=10, color=d3['Category20'][self.MAX_CLASSES][palette_key], alpha=0.5, source = source)
                legend_items.append((class_name, [class_points]))
                palette_key += 1
            legend = Legend(
                    items = legend_items,
                    click_policy = 'hide',
                    location = (0, -30)
            )
            p.add_layout(legend, 'right')
            p.output_backend = 'svg'
            vol_file = 'volcano_' + ratio_name
            vol_file_svg = vol_file + '.svg'
            vol_file_png = vol_file + '.png'
            vol_path_svg = self.root_path + vol_file_svg
            vol_path_png = self.root_path + vol_file_png
            export_svgs(p, filename=vol_path_svg)
            export_png(p, filename=vol_path_png)
            # save paths to zip
            self.paths_to_zip[vol_file_svg] = vol_path_svg
            self.paths_to_zip[vol_file_png] = vol_path_png
            plot_list.append([p])
        if plot_list:
           script, div = components(gridplot(plot_list))
        return script, div
