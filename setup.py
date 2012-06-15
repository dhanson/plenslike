#!/usr/bin/env python

import glob

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('plenslike',parent_package,top_path)
    config.add_extension('_plenslike', ['plenslike/plenslike.c'], undef_macros=['NDEBUG'])
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup

    setup(name='plenslike',
          packages=['plenslike'],
          package_data={'plenslike': glob.glob("data/*")},
          configuration=configuration)
