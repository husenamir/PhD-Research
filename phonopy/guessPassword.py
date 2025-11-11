# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 14:59:50 2023

@author: ahusen
"""

import datetime
import genetic

# **Helper functions for genetic engine**
def test_Hello_World():
  target = "Hello World!"
  guess_password(target)
  
def guess_password(target):
  geneset = " abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!."
  startTime = datetime. datetime. now()
  
  def fnGetFitness(genes):
      return get_fitness(genes, target)
  
  def fnDisplay(genes):  
      display(genes, target, startTime)
  optimalFitness = len(target)
  genetic. get_best(fnGetFitness, len(target), optimalFitness, geneset, fnDisplay)
  
# **Display**
def display(genes, target, startTime):
  timeDiff = datetime. datetime. now() - startTime
  fitness = get_fitness(genes, target)
  print("{}\t{}\t{}". format(genes, fitness, timeDiff))
  
# *Fitness Function**
def get_fitness(genes, target):
  return sum(1 for expected, actual in zip(target, genes) if expected == actual)

if __name__ == '__main__' :
  test_Hello_World()