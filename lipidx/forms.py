from flask_wtf import FlaskForm
from wtforms import (IntegerField, FloatField, StringField,
    TextAreaField, RadioField, BooleanField)
from wtforms.validators import Optional
from flask_wtf.file import FileField, FileRequired


class ElseOptional(Optional):
    def __init__(self, attr, val, *args, **kwargs):
        self.attr = attr
        self.val = val
        super(ElseOptional, self).__init__(*args, **kwargs)

    def __call__(self, form, field):
        if getattr(form, self.attr).data != self.val:
            super(ElseOptional, self).__call__(form, field)

class LipidAnalysisForm(FlaskForm):
    ION_GROUP_WITHIN_DEFAULT = 0.9
    RET_TIME_DEFAULT = 3
    GROUP_PQ_DEFAULT = 0.8
    GROUP_SN_DEFAULT = 100
    GROUP_AREA_DEFAULT = 0
    GROUP_HEIGHT_DEFAULT = 0
    MULT_FACTOR_DEFAULT = 3

    COLS_TO_REMOVE = ['ARatio', 'HRatio', 'ADiff', 'HDiff', 'GroupHeight', 'HeightRSD',
        'Height', 'NormArea', 'NormHeight', 'Hwhm(L)', 'Hwhm(R)', 'AreaScore', 'DataId', 'Scan',
        'It.', 'z', 'Delta(Da)', 'mScore', 'Occupy']

    file_msg = 'Must submit a file to process'
    file1 = FileField('File 1', [FileRequired()])
    file2 = FileField('File 2')
    group_ions_within = FloatField('Group ions with same lipid charge within ret time', default = ION_GROUP_WITHIN_DEFAULT)
    retention_time_filter = IntegerField('Retention Time', default =
            RET_TIME_DEFAULT)
    group_pq_filter = FloatField('GroupPQ', default = GROUP_PQ_DEFAULT)
    group_sn_filter = IntegerField('GroupS/N', default = GROUP_SN_DEFAULT)
    group_area_filter = IntegerField('Group Area', default = GROUP_AREA_DEFAULT)
    group_height_filter = IntegerField('Group Height', default =
            GROUP_HEIGHT_DEFAULT)
    blank = StringField('Name of the blank (exp. c, s1, s2)')
    mult_factor = IntegerField('Blank Multiplication Factor', default =
            MULT_FACTOR_DEFAULT)
    remove_cols = TextAreaField('Columns to remove (comma seperated)', default =
            ', '.join(COLS_TO_REMOVE))
    normalize = RadioField('Normalization', choices = [('none', 'None'),
    ('intensity', 'Sum of area intensities'), ('values', 'Enter values')], default = 'none')
    normal_c = StringField('c')
    normal_s1 = StringField('s1')
    normal_s2 = StringField('s2')
    normal_s3 = StringField('s3')
    normal_s4 = StringField('s4')
    normal_s5 = StringField('s5')
    normal_s6 = StringField('s6')
    normal_s7 = StringField('s7')
    normal_s8 = StringField('s8')
    normal_s9 = StringField('s9')
    normal_s10 = StringField('s10')
    class_stats = BooleanField('Class stats', default = True)
    group1 = StringField('group1')
    group2 = StringField('group2')
    group3 = StringField('group3')
    group4 = StringField('group4')
    group5 = StringField('group5')
    group6 = StringField('group6')

class VolcanoForm(FlaskForm):
    file_msg = 'Must submit a file to process'
    file1 = FileField('File 1', [FileRequired()])
    group1 = StringField('group1')
    group2 = StringField('group2')
    group3 = StringField('group3')
    group4 = StringField('group4')
    group5 = StringField('group5')
    group6 = StringField('group6')
