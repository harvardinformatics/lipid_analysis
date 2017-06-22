import unittest
import os
from lipidx.lipid_analysis import LipidAnalysis
from lipidx import constants as const
from lipidx import app

class LididTests(unittest.TestCase):

    def setUp(self):
        app.testing = True
        self.app = app.test_client()

    def tearDown(self):
        pass

    def test_lipid_analysis_init(self):
        res = False
        basedir = os.path.abspath( os.path.dirname( __file__ ) ) + '/'
        with app.app_context():
            la = LipidAnalysis([basedir + 'pos_short.txt', basedir + 'neg_short.txt'])
            # check that rows have expected count
            res = (len(la.rows) == 64)
        assert res

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

if __name__ == '__main__':
    unittest.main()
