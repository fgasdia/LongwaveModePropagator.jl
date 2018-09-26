# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 23:14:18 2018

@author: forrest
"""

import pyproj

rx = (40.014984, -105.270546)
tx = (44.646394, -67.281069)

SPHERE_G = pyproj.Geod('+a=6366200 +b=6366200')

SPHERE_G.inv(tx[1], tx[0], rx[1], rx[0])

