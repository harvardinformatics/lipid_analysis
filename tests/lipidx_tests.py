import unittest
import os
from lipidx.lipid_analysis import LipidAnalysis
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
            res = (len(la.rows) == 74)
        assert res

if __name__ == '__main__':
    unittest.main()
