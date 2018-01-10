import os
import logging

BASE_DIR = os.path.abspath( os.path.dirname( __file__ ) ) + '/lipidx/'
ADMINS = frozenset(os.environ.get('LIPIDX_ADMIN_EMAILS').split(','))
WTF_CSRF_ENABLED = True
SECRET_KEY = os.environ.get('LIPIDX_KEY','n`}UyR_r+9w2]%xZ~H?FRz^')
loglevel = os.environ.get('LIPIDX_LOGLEVEL','ERROR')
if loglevel == 'DEBUG':
    DEBUG = True
WTF_CSRF_SECRET_KEY = os.environ.get('WTF_CSRF_KEY',"aEsu'a}-j\>rJ4'8MFz{<yn")
UPLOAD_FOLDER = os.path.join(BASE_DIR, 'files/')
FILE_FOLDER = os.path.join(BASE_DIR, 'files/')
ALLOWED_EXTENSIONS = set(['txt', 'pdf', 'png', 'jpeg', 'gif', 'doc', 'xls', 'csv'])
MAX_CONTENT_LENGTH = 16 * 1024 * 1024

