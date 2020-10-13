#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 15:17:10 2020

@author: viola-j
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import nbinom

# Main class for a cluster, essentially a tree with pointers.
# Only twist is that each patient has a "center" [x, y] and generation
# which is probably a mistake
    
class infected:
    def __init__(self):
        self.parent = None
        self.children = []
        self.center = np.array([0, 0])
        self.generation = 0
    def copy(self, maxinfected = None):
        newinfected = infected()
        newinfected.center = 1*self.center
        newinfected.generation = self.generation
        for n, child in enumerate(self.children):
            if maxinfected is None or n < maxinfected:
                newchild = child.copy(maxinfected)
                newchild.parent = newinfected
                newinfected.children.append(newchild)
        return newinfected
    def countinfected(self):
        if len(self.children) == 0:
            return 1
        else:
            nchildren = 0
            for child in self.children:
                nchildren += child.countinfected()
            return nchildren + 1 # +1 for self!

def law_nbinom(n = .1, p = 0.038461538461538464):
    return lambda : nbinom.rvs(n, p)

def newctrs(ctr, Npts, r, avoid = None):
    if avoid is None:
        N = Npts
        th0 = 0
        xs = [ctr[0]+r*np.cos(2*np.pi*j/N+th0) for j in range(N)]
        ys = [ctr[1]+r*np.sin(2*np.pi*j/N+th0) for j in range(N)]
    else:
        N = Npts+1
        th0 = np.imag(np.log(avoid[0] + 1j*avoid[1]))
        xs = [ctr[0]+r*np.cos(2*np.pi*j/N+th0) for j in range(1,N)]
        ys = [ctr[1]+r*np.sin(2*np.pi*j/N+th0) for j in range(1,N)]
    return [np.array([xs[j], ys[j]]) for j in range(len(xs))]

def infectnew(patient, r, law = None):
    if law is None:
        law = law_nbinom()
    nchildren = law()
    if patient.parent is None:
        avoid = None
    else:
        avoid = patient.parent.center - patient.center
    newcs = newctrs(patient.center, nchildren, r, avoid)
    for j, ctr in enumerate(newcs):
        # print(j, ctr) # old diagnostic
        kid = infected()
        # print(kid.children) # old diagnostic
        kid.parent = patient
        patient.children.append(kid)
        kid.center = ctr
        kid.generation = patient.generation + 1
    return True

# Make a cluster, default the negative binomial
def makecluster(patient0, maxgen, law = None, rfn = None):
    if law is None:
        law = law_nbinom()
    if rfn is None:
        rfn = lambda g : 2**(-g) # default length is 1.5**{-gen}
    gen = patient0.generation
    infectious = [patient0]
    while gen < maxgen and len(infectious) > 0:
        nextinfectious = []
        for patient in infectious:
            infectnew(patient, rfn(gen), law = law)
            for newpatient in patient.children:
                nextinfectious.append(newpatient)
        gen += 1
        infectious = nextinfectious
    return True

def cluster_getpatients(patient0):
    allpatients = [patient0]
    for patient in allpatients:
        for child in patient.children:
            allpatients.append(child)
    return allpatients

def cluster_scatter(patients, ax, sizefn = None, colorfn = None, zorder = 1, maxgen = None):
    if maxgen is None:
        maxgen = max([patient.generation for patient in patients])
    if colorfn is None:
        if maxgen == 0:
            colorfn = lambda g: plt.cm.cool(0.)
        else:
            colorfn = lambda g: plt.cm.cool(g/maxgen)
    if sizefn is None:
        sizefn = lambda g: 200*1.5**(-g) # default size is 200...
    # populate scatter data
    xs, ys, sizes, colors = [], [], [], []
    for patient in patients:
        xs.append(patient.center[0])
        ys.append(patient.center[1])
        sizes.append(sizefn(patient.generation))
        colors.append(colorfn(patient.generation))
    ax.scatter(xs, ys, sizes, colors, zorder = zorder)

def cluster_edges(patients, ax, color = "black", zorder = 0):
    for patient in patients:
        for child in patient.children:
            ax.plot([patient.center[0], child.center[0]], [patient.center[1], child.center[1]], color = color, zorder = zorder)

def generate_tree(maxgen, minpatients = 0, law = None, rfn = None):
    numpatients = -1
    while numpatients < minpatients:
        patient0 = infected()
        makecluster(patient0, maxgen, law, rfn)
        patients = cluster_getpatients(patient0)
        numpatients = len(patients)
    return patient0

def plot_truncated_tree(patient0, maxgen, sizefn = None, colorfn = None, file = None):
    patients = cluster_getpatients(patient0)
    fig = plt.figure()
    ax = fig.gca()
    cluster_edges(patients, ax)
    cluster_scatter(patients, ax, sizefn, colorfn, maxgen = maxgen)
    plt.gca().axis('equal')
    if file is None:
        file = 'myplot.svg'
    plt.savefig(file)

def summarize_tree(patient0):
    patients = cluster_getpatients(patient0)
    numpatients = len(patients)
    maxgenattained = max([patient.generation for patient in patients])
    return [numpatients, maxgenattained]
