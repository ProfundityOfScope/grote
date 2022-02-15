#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 13:27:46 2020

@author: bruzewskis
"""
from PyOPT import _Container

class Scan():
    
    def __init__(self, target, resource):
        
        self.name = target #target.name
        self.target = target
        self.resource = resource
        
    def __str__(self):
        return self.name
    
    def pprint(self, _prefix=''):
        
        return _prefix + self.name + '\n'
    
class Loop(_Container._Container):
    
    def __init__(self, name, loopIters=1):
        super().__init__(name, Scan)
        
        self.loopIters = loopIters
        
        
        