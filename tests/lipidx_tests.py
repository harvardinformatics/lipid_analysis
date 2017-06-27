import unittest
import os
import csv
import json
from lipidx.lipid_analysis import LipidAnalysis
from lipidx.forms import LipidAnalysisForm as form
from lipidx import app

class LididxTests(unittest.TestCase):

    def setUp(self):
        app.testing = True
        self.app = app.test_client()
        self.sample_data_dir = (os.path.abspath( os.path.dirname( __file__ ) ) +
        '/sample_data/')

    def tearDown(self):
        pass

    def test_lipid_analysis_init_cnt(self):
        res = False
        with app.app_context():
            la = LipidAnalysis([self.sample_data_dir + 'neg_short.txt', self.sample_data_dir + 'pos_short.txt'])
            # check that rows have expected count
            res = (len(la.rows) == 64)
        assert res

    def test_lipid_analysis_init_diff(self):
        res = False
        with app.app_context():
            la = LipidAnalysis([self.sample_data_dir + 'neg_short.txt', self.sample_data_dir + 'pos_short.txt'])
            expected = self.csv_to_row_dict(self.sample_data_dir +
            'init_expected.csv')
            res, msg = self.diff_dicts(expected, la.rows)
        self.assertTrue(res, msg)

    def test_group_ions_diff(self):
        res = False
        with app.app_context():
            la = self.get_inst_from_step('init_expected.csv')
            expected = self.csv_to_row_dict(self.sample_data_dir +
            'group_ions_expected.csv')
            la.group_ions(form.ION_GROUP_WITHIN_DEFAULT)
            res, msg = self.diff_dicts(expected, la.rows)
        self.assertTrue(res, msg)

    def test_filter_rows_default(self):
        res = False
        with app.app_context():
            la = self.get_inst_from_step('group_ions_expected.csv')
            expected = self.csv_to_row_dict(self.sample_data_dir +
            'filter_expected.csv')
            la.filter_rows(form.RET_TIME_DEFAULT,
                    form.GROUP_PQ_DEFAULT,
                    form.GROUP_SN_DEFAULT,
                    form.GROUP_AREA_DEFAULT,
                    form.GROUP_HEIGHT_DEFAULT
            )
            res, msg = self.diff_dicts(expected, la.rows)
        self.assertTrue(res, msg)

    def test_subtract_blank(self):
        res = False
        with app.app_context():
            la = self.get_inst_from_step('filter_expected.csv')
            expected = self.csv_to_row_dict(self.sample_data_dir +
            'subtract_expected.csv')
            la.subtract_blank('c', form.MULT_FACTOR_DEFAULT)
            res, msg = self.diff_dicts(expected, la.rows)
        self.assertTrue(res, msg)

    def test_remove_columns(self):
        res = False
        with app.app_context():
            la = self.get_inst_from_step('subtract_expected.csv')
            expected = self.csv_to_row_dict(self.sample_data_dir +
            'remove_expected.csv')
            la.remove_columns(', '.join(form.COLS_TO_REMOVE))
            test_cols = la.rows[next(iter(la.rows))].keys()
            expected_cols = expected[next(iter(expected))].keys()
            same, msg = self.diff_keys(expected_cols, test_cols)
        self.assertTrue(same, msg)

    def test_normalize(self):
        res = False
        with app.app_context():
            la = self.get_inst_from_step('remove_expected.csv')
            expected = self.csv_to_row_dict(self.sample_data_dir +
            'normalize_expected.csv')
            la.normalize({'normalize': 'intensity'})
            res, msg = self.diff_dicts(expected, la.rows)
        self.assertTrue(res, msg)

    def test_class_stats(self):
        res = False
        with app.app_context():
            # no normalization is the default, use remove_expected to init rows
            la = self.get_inst_from_step('remove_expected.csv')
            expected = self.csv_to_row_dict(self.sample_data_dir +
            'class_stats_expected.csv', 'class')
            la.calc_class_stats({'class_stats':'y'})
            res, msg = self.diff_dicts(expected, la.class_dict)
        self.assertTrue(res, msg)

    def test_subclass_stats(self):
        res = False
        with app.app_context():
            # no normalization is the default, use remove_expected to init rows
            la = self.get_inst_from_step('remove_expected.csv')
            expected = self.csv_to_row_dict(self.sample_data_dir +
            'subclass_stats_expected.csv', 'subclass')
            la.calc_class_stats({'class_stats':'y'})
            res, msg = self.diff_dicts(expected, la.subclass_dict)
        self.assertTrue(res, msg)

    # helper function get csv file into same format as rows for diff
    def csv_to_row_dict(self, csv_file, key_on = 'name'):
        row_dict = {}
        with open(csv_file, mode='r') as f:
            rows = csv.DictReader(f)
            row_dict = {r[key_on]: r for r in rows}
        return row_dict

    def diff_keys(self, exp_keys, test_keys):
        exp_keys = set(exp_keys)
        test_keys = set(test_keys)
        not_expected = test_keys - exp_keys
        missing = exp_keys - test_keys
        same = (missing == set({}) and not_expected == set({}))
        msg = ''
        if not same:
            msg = ('not expected: ' + json.dumps(list(not_expected)) + '\n missing: '
                    + json.dumps(list(missing)))
        return same, msg

    def diff_dicts(self, exp, test):
        same, msg = self.diff_keys(exp.keys(), test.keys())
        if same:
            same, msg = self.diff_rows(exp, test)
        return same, msg

    def diff_rows(self, exp, test):
        exp_keys = set(exp.keys())
        errors = {}
        for key, row in exp.items():
            for col, val in row.items():
                if col:
                    if val and col not in test[key]:
                        if col not in errors:
                            errors[col] = {'missing':[], 'diff':[]}
                        errors[col]['missing'].append(key)
                    elif col in test[key]:
                        val_str = str(val)
                        test_str = str(test[key][col])
                        if len(val_str) > 1 and len(test_str) > 1:
                            val_str = val_str[0:-2]
                            test_str = test_str[0:-2]
                        if test_str != val_str:
                            if col not in errors:
                                errors[col] = {'missing':[], 'diff':[]}
                            errors[col]['diff'].append(key + ' ' + str(val) + ' vs ' +
                                    str(test[key][col]))

        same = (errors == {})
        msg = ''
        if errors:
            msg = 'row errors: ' + json.dumps(errors)
        return same, msg

    def get_inst_from_step(self, prev_step_file):
        # pass empty path and then fill rows with prev step results
        la = LipidAnalysis([])
        rows = self.csv_to_row_dict(self.sample_data_dir + prev_step_file)
        la.rows = rows
        return la

if __name__ == '__main__':
    unittest.main()
