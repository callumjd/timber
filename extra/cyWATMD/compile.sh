#!/bin/bash

cython example.pyx
python setup.py build_ext --inplace
