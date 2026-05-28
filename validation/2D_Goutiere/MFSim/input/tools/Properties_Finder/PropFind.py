#!/usr/bin/env python
# -*- coding:utf-8 -*-

#'''
#--------------------------------------
#Fluid Properties Finder (PropFind)
#
#Uses the cantera library to get the 
#physical properties of a fluid
#
#Python module Dependencies:
#`commentjson`
#`cantera`
#
#--------------------------------------
#'''

# Import System libraries
import cantera as cantera

# Import software libraries
import lib.SearchCantera as search


search.start()
search.input()
search.set_cantera()
search.get_cantera()
search.end()

