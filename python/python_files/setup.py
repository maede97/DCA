from setuptools import setup,find_packages

import sys
if sys.version_info < (3,0):
  sys.exit('Sorry, Python < 3.0 is not supported')

setup(
  name        = 'PythonDCA',
  version     = '0.1',
  packages    = find_packages(),
)
