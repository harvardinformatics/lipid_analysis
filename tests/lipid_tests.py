import lipid_analysis.lipid_app
import unittest

class LididTests(unittest.TestCase):

    def setUp(self):
        lipid_app.app.testing = True
        self.app = lipid_app.app.test_client()

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
