import csv
import os
import numpy
from scipy.stats import ttest_ind
import zipfile
import re
from flask import request
from flask import current_app as app
from collections import OrderedDict
import logging
from bkcharts import Bar
from bokeh.layouts import gridplot
from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.models import HoverTool
from bokeh.embed import components
from bokeh.sampledata.autompg import autompg as df

class LipidAnalysis:

    ROUND_TO = 2
    POST_NORMAL_ROUND = 8

    def __init__ (self, paths, debug = False):
        self.paths = paths
        # debug adds cols to results to show pre normalized values
        self.debug = debug
        self.area_start = 'Area['
        self.group_area_start = 'GroupArea['
        self.groups = {}
        # cols should be the same in all files
        # col names will be taken from first file
        self.rows = self.get_rows_from_files(self.paths)
        # filled in calc_class_stats
        self.class_stats = {}
        self.subclass_stats = {}
        self.class_dict = {}
        self.subclass_dict = {}

        # file paths, eventualy these may not be hardcoded
        self.root_path = app.config['UPLOAD_FOLDER'] + '/'
        lipid_class_file = 'lipidKey.csv'
        self.lipid_class_path = app.config['BASE_DIR'] + '/' + lipid_class_file
        self.lipid_results_file = 'lipid_analysis.csv'
        self.lipid_results_path = self.root_path + self.lipid_results_file
        self.subclass_file = 'subclass_stats.csv'
        self.subclass_path = self.root_path + self.subclass_file
        self.class_file = 'class_stats.csv'
        self.class_path = self.root_path + self.class_file
        zip_file = 'lipid_results.zip'
        self.zip_path = self.root_path + zip_file

    def get_rows_from_files(self, paths):
        rows = OrderedDict()
        cols = []
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
                            row_cols = ln.split('\t')
                        else: # data lines
                            ln = ln.replace('\n', '')
                            row = ln.split('\t')
                            # remove trailing newline
                            row[(len(row) - 1)] = row[(len(row) - 1)].strip('\n')
                            row_d = OrderedDict(zip(row_cols, row))
                            # cols from all files must be the same
                            row_d = self.limit_row_cols(cols, row_d)
                            # calc retention time: average of GroupTopPos
                            ret_time = round(numpy.mean(self.list_col_type(row_d, 'GroupTopPos')), self.ROUND_TO)
                            row_d['ret_time'] = ret_time # add to row
                            row_d.move_to_end('ret_time', last=False)
                            # unique name for row LipidIon + ret_time
                            name = row_d['LipidIon'] + '_' + str(ret_time)
                            row_d['name'] = name
                            row_d.move_to_end('name', last=False)
                            if name in rows: # rare case
                                # if lipid has same name then keep the one with
                                # greater area
                                avg_areas = self.list_col_type(row_d, self.area_start)
                                avg_prev_areas = self.list_col_type(rows[name], self.area_start)
                                if avg_areas > avg_prev_areas:
                                    rows[name] = row_d
                            else: # add the new row
                                rows[name] = row_d
                if not cols:
                    cols = row_cols
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
        # make sure results are sorted by key
        res = [x for y, x in sorted(self.rows.items(), key=lambda t: t[0].lower())]
        self.write_csv(self.lipid_results_path, self.get_cols(), res)
        # create a zip file for lipids and stats
        z = zipfile.ZipFile(self.zip_path, "w")
        z.write(self.lipid_results_path, self.lipid_results_file)
        if os.path.exists(self.class_path):
            z.write(self.class_path, self.class_file)
        if os.path.exists(self.subclass_path):
            z.write(self.subclass_path, self.subclass_file)
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
            self.rows = subtracted

    def calculate_avg_blank(self, blank_cols, row):
            blank_vals = []
            for col in blank_cols:
                blank_vals.append(float(row[col]))
            return round(numpy.mean(blank_vals), self.ROUND_TO)

    def remove_columns(self, remove_cols):
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
            grps = re.search('(.*[+,-]).*', name)
            lipid_charge = grps.group(1)
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
            if row['Rej.'] == '0':
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
        group_pq_max = max(self.list_col_type(row, 'GroupPQ'))
        if group_pq_max <= group_pq_fil:
            return False
        group_sn_max = max(self.list_col_type(row, 'GroupS/N'))
        if group_sn_max <= group_sn_fil:
            return False
        if group_area_fil > 0: # all values are pos ints, so skip if 0
            group_area_max = max(self.list_col_type(row, 'GroupArea'))
            if group_area_max <= group_area_fil:
                return False
        if group_height_fil > 0: # all values are pos ints, so skip if 0
            group_height_max = max(self.list_col_type(row, 'GroupHeight'))
            if group_height_max <= group_height_fil:
                return False
        return True

    def normalize(self, form_data):
        # most common is to not normalize
        # or one can use avg intensity calculated from data
        # or input manual values
        if self.rows:
            normal = self.rows
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
                                if self.debug: # put old values in rows to debug
                                    normal[name][col + 'old'] = row[col]
                                    normal[name][col + 'div'] = float(form_data[form_name])
                                normal[name][col] = round(row[col] / float(form_data[form_name]), self.POST_NORMAL_ROUND)
                # use calculated intensity
                elif form_data['normalize'] == 'intensity':
                    intensities = self.calc_intensities(area_cols)
                    for name, row in normal.items():
                        for col in area_cols:
                            sam = self.get_sample_from_col(col)
                            # TODO: error case if no intensity?
                            if intensities[sam] > 0:
                                if self.debug:
                                    normal[name][col + 'old'] = row[col]
                                    normal[name][col + 'div'] = intensities[sam]
                                # TODO: 8 dec place
                                normal[name][col] = round(float(row[col])/intensities[sam],
                                self.POST_NORMAL_ROUND)
                self.recalc_cols()
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

    def recalc_cols(self):
        self.groups = self.get_groups()
        for name, row in self.rows.items():
            stats = OrderedDict()
            # for each group recalc the avg and std from areas
            for group, nums in self.groups.items():
                if group not in stats:
                    stats[group] = []
                for num in nums:
                    num_col = self.area_start + group + '-' + num + ']'
                    stats[group].append(float(row[num_col]))
            for group, val_lst in stats.items():
                self.rows[name]['GroupAVG[' + group + ']'] = round(numpy.mean(val_lst), self.POST_NORMAL_ROUND)
                self.rows[name]['GroupRSD[' + group + ']'] = round(numpy.std(val_lst), self.POST_NORMAL_ROUND)

    def calc_class_stats(self, opt_class_stats):
        sub_success = True
        class_success = True
        # if class stats are requested
        if opt_class_stats:
            self.class_keys = self.load_lipid_classes()
            class_stats = {}
            subclass_stats = {}
            for name, row in self.rows.items():
                # take subclass key from row
                subclass_key = row['Class']
                if subclass_key in self.class_keys:
                    # get corresponding names from class_keys
                    subclass_name = self.class_keys[subclass_key]['subclass']
                    class_name = self.class_keys[subclass_key]['class']
                    # populate new subclasses and classes
                    if subclass_name not in subclass_stats:
                        subclass_stats[subclass_name] = {}
                        if class_name not in class_stats:
                            class_stats[class_name] = {}
                    # add row to group stats for the class and subclass that
                    # correspond to the row
                    subclass_stats[subclass_name] = self.group_stats(row,
                            subclass_stats[subclass_name])
                    class_stats[class_name] = self.group_stats(row,
                            class_stats[class_name])
            # write files
            subclass_cols = self.stats_cols('subclass')
            self.subclass_stats, self.subclass_dict = self.format_stats('subclass', subclass_stats)
            sub_success = self.write_csv(self.subclass_path, subclass_cols, self.subclass_dict.values())
            class_cols = self.stats_cols('class')
            self.class_stats, self.class_dict = self.format_stats('class', class_stats)
            class_success = self.write_csv(self.class_path, class_cols, self.class_dict.values())
        return sub_success, class_success

    def stats_cols(self, cat):
        cols = [cat]
        grps = self.get_groups()
        for g in grps:
            cols.append(g + ' cnt')
            cols.append(g + ' avg')
            cols.append(g + ' std')
        return cols

    def format_stats(self, cat, stats):
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
                row[group + ' std'] = std
                stats[name][group]['std'] = std
            rows[name] = row
        return stats, rows

    def group_stats(self, row, grp_info):
        prefix = 'GroupArea'
        if not self.groups:
            self.groups = self.get_groups()
        for key in self.groups.keys():
            if key not in grp_info:
                # keep cnt and a list of group areas for group
                grp_info[key] = {'cnt': 0, 'grp_areas': []}
            areas = self.list_col_type(row, self.area_start + key)
            if max(areas) > 0.0:
                grp_info[key]['cnt'] += 1
            grp_info[key]['grp_areas'].append(numpy.mean(areas))
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
        data = {
                'lipid': [],
                'group': [],
                'nb': [],
                'avg': [],
                'sum': [],
                'std': []
        }
        for lipid, groups in self.class_stats.items():
            for group, stats in groups.items():
                data['lipid'].append(lipid)
                data['group'].append(group)
                data['nb'].append(stats['cnt'])
                data['avg'].append(stats['avg'])
                data['sum'].append(stats['sum'])
                data['std'].append(stats['std'])
        bar_cnt = Bar(data, title = 'nb of lipids', ylabel = 'nb of lipids', values = 'nb', label ='lipid', group = 'group')
        bar_avg = Bar(data, title = 'intensity', ylabel = 'sum of area per group', values = 'sum', label ='lipid', group = 'group')
        bars = gridplot(
                [bar_cnt, None],
                [bar_avg, None]
        )
        script, div = components(bars)
        return script, div

    def calc_ratio(self):
        for key, row in self.rows.items():
            ratio = float(row['GroupArea[s2]'])/float(row['GroupArea[s1]'])
            log_ratio = numpy.log2(ratio)
            self.rows[key]['ratio'] = ratio
            self.rows[key]['log_ratio'] = log_ratio
            s2 = self.list_col_type(row, self.area_start + 's2')
            s1 = self.list_col_type(row, self.area_start + 's1')
            t, p = ttest_ind(s2, s1)
            self.rows[key]['p_value'] = p

    def volcano_plot(self):
        self.calc_ratio()
        data = {
                'lipid': [],
                'log2': [],
                'p': []
        }
        for key, row in self.rows.items():
            data['lipid'].append(('Name', key))
            data['log2'].append(row['log_ratio'])
            data['p'].append(row['p_value'])
        source = ColumnDataSource(data = data)
        hover = HoverTool(tooltips=[
            ('name', "@lipid")
        ])
        p = figure(title='test', tools=[hover], x_axis_label = 'log2(ratio)', y_axis_label = 'p value')
        p.circle('log2', 'p', size=10, color="red", legend='Temp.', alpha=0.5,
                source = source)
        #hover = p.select_one(HoverTool)
        #hover.point_policy = "follow_mouse"
        script, div = components(p)
        return script, div
