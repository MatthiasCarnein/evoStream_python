from distutils.core import setup, Extension

setup(name='evoStreamPkg', version='1.0',  \
      ext_modules=[Extension('evoStream', ['evoStreamWrapper.cpp', 'evoStream.cpp', 'MC.cpp'])])
