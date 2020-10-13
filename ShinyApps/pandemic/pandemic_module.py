#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 15:17:10 2020

@author: viola-j
"""

import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import nbinom
from scipy.stats import poisson

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

# This is just for memory: since E(X) = k(1/p - 1)...
def nbinom_nR0_to_p(n, R0):
    return 1/(1+R0/n)

def law_nbinom(n = .1, p = 0.038461538461538464):
    return lambda : nbinom.rvs(n, p)

# This is the most frequently called law.
def law_nbinom_nR0(n = .1, R0 = 2.5):
    p = nbinom_nR0_to_p(n, R0)
    return lambda : nbinom.rvs(n, p)

def law_poisson(lam = 2.5):
    return lambda : poisson.rvs(lam)

# For people less experienced with Python, the "probability law" for a constant
def law_constant(n = 3):
    return lambda : n

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

# Not completed, I'm not sure that there's not information missing.
def cluster_saveparentage(patient0):
    patients = cluster_getpatients(patient0)
    whoparents = np.zeros((len(patients),), dtype = int)
    whoparents[0] = -1
    for npat, patient in enumerate(patients):
        for kpat in range(npat):
            if patient.parent is patients[kpat]:
                whoparents[npat] = kpat
    return whoparents

# Based on assumption that, as in saveparentage, whoparent is an increasing
# numpy array of numbers of parents
def cluster_reconstruct(whoparents, rfn = None):
    if rfn is None:
        rfn = lambda g : 2**(-g)
    patients = [infected() for j in range(len(whoparents))]
    for npat, patient in enumerate(patients):
        nchildren = np.sum([whoparents == npat])
        if patient.parent is None:
            avoid = None
        else:
            avoid = patient.parent.center - patient.center
        r = rfn(patient.generation)
        newcs = newctrs(patient.center, nchildren, r, avoid)
        whichchild = 0
        for kpat in range(npat+1, len(whoparents)):
            if whoparents[kpat] == npat:
                patients[kpat].parent = patient
                patient.children.append(patients[kpat])
                patients[kpat].generation = patient.generation + 1
                patients[kpat].center = newcs[whichchild]
                whichchild += 1
        if not whichchild == nchildren:
            print("Shit, I expected number of children to be {}, on patient number {}.".format(nchildren, npat))
    return patients

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

def cluster_markends(patients, ax, endgen, sizefn = None, zorder = 1.1):
    if sizefn is None:
        sizefn = lambda g: 200*1.5**(-g) # default size is 200...
    # populate scatter data
    xs, ys, sizes = [], [], []
    for patient in patients:
        if len(patient.children) == 0 and patient.generation < endgen:
            xs.append(patient.center[0])
            ys.append(patient.center[1])
            sizes.append(sizefn(patient.generation)/3)
    ax.scatter(xs, ys, sizes, color = "white", zorder = zorder)

def cluster_edges(patients, ax, color = "black", zorder = 0):
    for patient in patients:
        for child in patient.children:
            ax.plot([patient.center[0], child.center[0]], [patient.center[1], child.center[1]], color = color, zorder = zorder)

def main_cluster_draw(maxgen, ax = None, minpatients = 0, law = None, rfn = None, sizefn = None, colorfn = None):
    if ax is None:
        fig = plt.figure()
        ax = fig.gca()
    numpatients = -1
    while numpatients < minpatients:
        patient0 = infected()
        makecluster(patient0, maxgen, law, rfn)
        patients = cluster_getpatients(patient0)
        numpatients = len(patients)
    cluster_scatter(patients, ax, sizefn, colorfn)
    # cluster_markends(patients, ax, maxgen, sizefn)
    cluster_edges(patients, ax)
    plt.gca().axis('equal')
    plt.savefig('~/myplot.svg')
    maxgenattained = max([patient.generation for patient in patients])
    report = "drew {} patients over {} generations ({} possible)".format(len(patients), maxgenattained, maxgen)
    # print(report)
    return patient0

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

def test_nbinom(n = 1, p = 1):
    mean, var, skew, kurt = nbinom.stats(n, p, moments='mvsk')
    print("This negative binomial distribution has mean {:.5f} and variance {:.5f}.".format(mean, var))
    print("Here are a bunch of random draws!")
    return nbinom.rvs(n, p, size=100)

"""
#### Old versions ####

def nbinom_infect(patient, n, p, r):
    nchildren = nbinom.rvs(n, p)
    # print(nchildren) # old diagnostic
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

# Make a cluster with negative binomial
def makecluster(patient0, maxgen, n, p, law = None, rfn = None):
    if law = None:
    if rfn is None:
        rfn = lambda g : 2**(-g) # default length is 1.5**{-gen}
    gen = patient0.generation
    infectious = [patient0]
    while gen < maxgen and len(infectious) > 0:
        nextinfectious = []
        for patient in infectious:
            nbinom_infect(patient, n, p, rfn(gen))
            for newpatient in patient.children:
                nextinfectious.append(newpatient)
        gen += 1
        infectious = nextinfectious
    return True
"""
