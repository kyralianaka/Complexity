# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 11:41:16 2019

"""
from distutils.core import setup,Command
from distutils.extension import Extension
import setuptools,platform,os

if platform.system() == 'Darwin':
    from distutils import sysconfig
    vars = sysconfig.get_config_vars()
    vars['LDSHARED'] == vars['LDSHARED'].replace('-bundle','-dynamicLib')
    os.environ["CC"] = "clang"
    os.environ["CCX"] = "clang"

#Changed because /lzc not cross-platform compatible
liblzc = Extension("lzc",sources= ["lzc/lzc.c"])

setup(name = 'complexity',
      version = '0.2.0',
      description = 'python implementation of sequence complexity measures',
      author = "Kyra Kadhim, Kevin Brown",
      author_email = "kadhimk@oregonstate.edu",
      packages = ['complexity','complexity.lzc','complexity.ncd'],
      package_dir = {'complexity': ''},
      ext_package = 'complexity',
      ext_modules = [liblzc],
      license = 'BSD-3',
      classifiers=[
            'License :: OSI Approved :: BSD-3 License',
            'Intended Audience :: Developers',
            'Intended Audience :: Scientists',
            'Programming Language :: Python',
            'Topic :: Dynamical Systems',
            'Topic :: Statistics'
      ],
     )
