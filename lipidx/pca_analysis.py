import numpy
import re
from collections import OrderedDict
from bokeh.plotting import figure
from bokeh.models import Legend
from bokeh.palettes import d3
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from bokeh.models import HoverTool


class PCAAnalysis:

    def __init__(self, path):
        self.path = path
        self.area_start = 'area['
        self.rows = self.get_rows_from_files()

    def get_rows_from_files(self):
        rows = OrderedDict()
        with open(self.path, 'r') as f:
            for i, ln in enumerate(f):
                if i == 0:
                    row_cols = []
                ln = ln.replace('\n', '')
                file_row = re.split(',', ln)
                if not row_cols:
                    row_cols = file_row
                else:  # data lines
                    row = file_row
                    row_d = OrderedDict(zip(row_cols, row))
                    name = row_d['name']
                    rows[name] = row_d
        return rows

    def pca(self):
        cols = self.get_sample_cols()
        data = self.get_sample_data(cols)
        sample_grps = self.get_sample_groups(cols)
        pca = PCA(n_components=2)
        # scale and center the data (though it seems fit_transform centers when
        # used alone)
        std_data = StandardScaler().fit_transform(data)
        data_pca = pca.fit_transform(std_data)
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
        source = {}
        for grp, data in sample_grps.items():
            indexes = data['index']
            source['pc_1'] = [x for i, x in enumerate(pca.components_[0]) if i in indexes]
            source['pc_2'] = [y for i, y in enumerate(pca.components_[1]) if i in indexes]
            source['name'] = data['sam']
            class_points = p.circle('pc_1', 'pc_2', size=10, color=d3['Category20'][20][palette_key], source = source)
            palette_key += 1
            legend_items.append((grp, [class_points]))
        hover = HoverTool(tooltips=[
            ('name', "@name")
        ])
        p.add_tools(hover)
        legend = Legend(
                items = legend_items,
                click_policy = 'hide',
                location = (0, -30)
        )
        p.add_layout(legend, 'right')
        return p

    def get_sample_cols(self):
        area_cols = self.get_cols(self.area_start)
        # if using template instead of lipid_analysis results then use all cols
        # but name
        if not area_cols:
            area_cols = list(next(iter(self.rows.values())).keys())[1:]
        return area_cols

    def get_sample_data(self, area_cols):
        sample_grps = {}
        data = []
        for row in self.rows.values():
            area_data = self.limit_row_cols(area_cols, row)
            area_data = [float(x) for x in area_data.values()]
            data.append(area_data)
            for col in area_cols:
                prefix = col.split('_')[0]
                if prefix not in sample_grps:
                    sample_grps[prefix] = []
                sample_grps[prefix].append(row[col])
        return data

    def get_sample_groups(self, area_cols):
        sample_grps = OrderedDict()
        for i, col in enumerate(area_cols):
            name = col.replace(self.area_start, '').replace(']', '')
            prefix = name.split('-')[0]
            if prefix not in sample_grps:
                sample_grps[prefix] = {'index': [], 'sam': []}
            sample_grps[prefix]['index'].append(i)
            sample_grps[prefix]['sam'].append(name)
        return sample_grps

    def limit_row_cols(self, cols, row):
        if cols:
            extra = set(row.keys()) - set(cols)
            for k in extra:
                del row[k]
        return row

    def get_cols(self, start=None):
        first = list(self.rows.keys())[0]
        keys = list(self.rows[first].keys())
        if start:
            keys = [i for i in keys if i.startswith(start)]
        return keys
