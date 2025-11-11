# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 11:47:59 2023

@author: ahusen
"""

import random as rd
import numpy as np
import datetime as dt
geneSet = " abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!."
target = "Hello World!"

# ** Parents selection ** 
def generate_parent(length):
  genes = []
  while len(genes) < length:
      sampleSize = min(length - len(genes), len(geneSet))
      genes.extend(rd.sample(geneSet, sampleSize))
  return ''.join(genes)
#list.extend() to appends multiple items to a list
#string.join() uses the given string as a separator

# ** Fitness function **

def get_fitness(guess):
  return sum(1 for expected, actual in zip(target, guess) if expected == actual)
#zip() is a built-in function that makes it possible to iterate over two lists

# ** Mutation **

def mutate(parent):
  index = rd.randrange(0, len(parent)) #randrange(start, stop, int_step)
  childGenes = list(parent)
  newGene, alternate = rd.sample(geneSet, 2)
  childGenes[index] = alternate if newGene == childGenes[index] else newGene
  return ''.join(childGenes)

# ** Display **
def display(guess):
  timeDiff = dt.datetime.now() - startTime
  fitness = get_fitness(guess)
  print("{}\t{}\t{}". format(guess, fitness, timeDiff))

# ** Maain **
startTime = dt.datetime.now()
rd.seed(1)
bestParent = generate_parent(len(target))
bestFitness = get_fitness(bestParent)

# **Final Part**

while True:
  child = mutate(bestParent)
  childFitness = get_fitness(child)
  if bestFitness >= childFitness:
      continue
  display(child)
  if childFitness >= len(bestParent):
      break
  bestFitness = childFitness

  bestParent = child

