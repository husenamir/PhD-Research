# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 14:29:13 2023

@author: ahusen
"""

# **Creating Engine**


import random as rd

# ** Generating parents**
def _generate_parent(length, geneSet):
  genes = []
  while len(genes) < length:
      sampleSize = min(length - len(genes), len(geneSet))
      genes.extend(rd.sample(geneSet, sampleSize))
  return ''.join(genes)


# ** Mutation**
def _mutate(parent, geneSet):
  index = rd.randrange(0, len(parent)) #randrange(start, stop, int_step)
  childGenes = list(parent)
  newGene, alternate = rd.sample(geneSet, 2)
  childGenes[index] = alternate if newGene == childGenes[index] else newGene
  return ''.join(childGenes)

"""
Main Loop:
Parameters are:
• the function it calls to request the fitness for a guess,
• the number of genes to use when creating a new gene sequence,
• the optimal fitness value,
• the set of genes to use for creating and mutating gene sequences, and
• the function it should call to display, or report, each improvement found.
"""
def get_best(get_fitness, targetLen, optimalFitness, geneSet, display):
  rd. seed()
  bestParent = _generate_parent(targetLen, geneSet)
  bestFitness = get_fitness(bestParent)
  display(bestParent)
  if bestFitness >= optimalFitness:
      return bestParent
  while True:
      child = _mutate(bestParent, geneSet)
      childFitness = get_fitness(child)
      if bestFitness >= childFitness:
          continue
      display(child)
      if childFitness >= optimalFitness:
          return child
      bestFitness = childFitness
      bestParent = child