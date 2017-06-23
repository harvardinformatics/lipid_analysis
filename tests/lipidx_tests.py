import unittest
import os
import csv
import json
from lipidx.lipid_analysis import LipidAnalysis
from lipidx import constants as const
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
            'init_results.csv')
            res, msg = self.diff_dicts(expected, la.rows)
        self.assertTrue(res, msg)

    def test_filter_rows_default(self):
        res = False
        basedir = os.path.abspath( os.path.dirname( __file__ ) ) + '/'
        with app.app_context():
            la = LipidAnalysis([basedir + 'pos_short.txt', basedir + 'neg_short.txt'])
            la.filter_rows(const.RET_TIME_DEFAULT,
                    const.GROUP_PQ_DEFAULT,
                    const.GROUP_SN_DEFAULT,
                    const.GROUP_AREA_DEFAULT,
                    const.GROUP_HEIGHT_DEFAULT
            )
            # confirm filtered row count is correct
            res = (len(la.rows) == 24)
        assert res

    def test_subtract_blank(self):
        res = False
        basedir = os.path.abspath( os.path.dirname( __file__ ) ) + '/'
        with app.app_context():
            la = LipidAnalysis([basedir + 'pos_short.txt', basedir + 'neg_short.txt'])
            la.filter_rows(const.RET_TIME_DEFAULT,
                    const.GROUP_PQ_DEFAULT,
                    const.GROUP_SN_DEFAULT,
                    const.GROUP_AREA_DEFAULT,
                    const.GROUP_HEIGHT_DEFAULT
            )
            la.subtract_blank('c', const.MULT_FACTOR_DEFAULT)
            print(len(la.rows))
            # confirm filtered row count is correct
            res = (len(la.rows) == 24)
        assert res

    # helper function get csv file into same format as rows for diff
    def csv_to_row_dict(self, csv_file):
        row_dict = {}
        with open(csv_file, mode='r') as f:
            rows = csv.DictReader(f)
            row_dict = {r['name']: r for r in rows}
        return row_dict

    def diff_dicts(self, exp, test):
        exp_keys = set(exp.keys())
        test_keys = set(test.keys())
        not_expected = test_keys - exp_keys
        missing = exp_keys - test_keys
        same = (missing == set({}) and not_expected == set({}))
        msg = ''
        if same:
            same, msg = self.diff_rows(exp, test)
        else:
            msg = ('not expected: ' + json.dumps(list(not_expected)) + '\n missing: '
                    + json.dumps(list(missing)))
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
                    elif col in test[key] and str(test[key][col]) != str(val):
                        if col not in errors:
                            errors[col] = {'missing':[], 'diff':[]}
                        errors[col]['diff'].append(key)

        same = (errors == {})
        msg = ''
        if errors:
            msg = 'row errors: ' + json.dumps(errors)
        return same, msg

if __name__ == '__main__':
    unittest.main()
