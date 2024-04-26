#!/usr/bin/env python
import sympy
import numpy
import math
from sympy import nsimplify

""" 
Analytic solution to the n-site linear chain steady flow problem as described in 
   M Breuer, K Rosso, and J Blumberger, PNAS, 111(2), 2013 (10.1073/pnas.1316156111)

Usage: 
1. Go to the section "configuration" on top of this file. 
   a. Define the number of sites with the variable `hops` by giving it as many single-character sites
   as your system has.
   b. Define the input and output rate.
   c. Give the rates in the order of the sites. The first tuple element is the forward direction, 
   the second one the backward direction.
3. Run the script with
     $ python hopping.py
     
Note: Make sure that the python modules sympy and numpy are installed.

"""
kforward, kbackward, ps, replacements = dict(), dict(), dict(), dict()
# prepare symbols with their latex representation
kin, kout = sympy.symbols(r'k_\text{in} k_\text{out}')
j = sympy.symbols('J')
    
# helper functions
def build_jji(start, stop):
        """ Creates the base equation for a hop from site `start` to `stop`. """
        pj, pi = ps[stop], ps[start]
        kji = kforward[stop+start]
        kij = kbackward[start+stop]
        return sympy.Eq(j, kji*pi*(1-pj)-kij*pj*(1-pi))

def pplusN(start, stop, pstartN, pstartD, p0): # numerator of pplus
        # from p_start to p_stop, return p_stop
	kji = kforward[stop+start]
	kij = kbackward[start+stop]
	upper = kji*pstartN - kin*(1-p0)*pstartD
        #print "upper", upper
	return upper

def pplusD(start, stop, pstartN, pstartD, p0): # denominator of pplus
        # from p_start to p_stop, return p_stop
	kji = kforward[stop+start]
	kij = kbackward[start+stop]
	lower = kij*pstartD + pstartN*(kji - kij)
        # print "pstartD", pstartD
        # print "pstartN", pstartN
        # print "lower", lower
	return lower

# Function to solve maximun flux J_max with given sets of heme-heme ET rates.
def solve_flux(kfor, kback):
    # kji = [k21, k32, k43]; kij = [k12, k23, k34]

    # configuration
    #8 heme sites A, B, C, D, E, F, G, H.
    hopsconfig = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    hops = hopsconfig[:len(kfor)+1]
    #print hops, len(hops)
    #input and output rates are in excess; for maximum flux calculation.
    input_rate = 1e12
    output_rate = 1e12
    #forward and backward rates ((k_BA, k_AB), (k_CB, k_BC), (k_DC, k_CD))...
    #k21, k32, k43, k54, k65, k76, k87 = kfor[0], kfor[1], kfor[2], kfor[3], kfor[4], kfor[5], kfor[6]
    #k12, k23, k34, k45, k56, k67, k78 = kback[0], kback[1], kback[2], kback[3], kback[4], kback[5], kback[6]
    #rates = ((k21, k12), (k32, k23), (k43, k34), (k54, k45), (k65, k56), (k76, k67), (k87, k78)) 
    rates = tuple(zip(kfor, kback))
    #print rates, len(rates)
    
    assert(len(hops) == len(rates) + 1)
    
    #kforward, kbackward, ps, replacements = dict(), dict(), dict(), dict()
    replacements[kin] = input_rate
    replacements[kout] = output_rate
    for start, stop, rate in zip(hops[:-1], hops[1:], rates):
            kf = sympy.symbols(r'k_\text{%s%s}' % (stop, start))
            kforward['%s%s' % (stop, start)] = kf
            kb = sympy.symbols(r'k_\text{%s%s}' % (start, stop))
            kbackward['%s%s' % (start, stop)] = kb
            replacements[kf] = rate[0]
            replacements[kb] = rate[1]
    for hop in hops:
            ps[hop] = sympy.symbols(r'P_\text{%s}' % hop)
    
    # first and last half-complete flux relations
    population_expressionsN = [] # for numerator
    population_expressionsD = [] # for denominator
    population_expressionsN.append(1.0 - j/kin)
    population_expressionsD.append(1.0)
    #initial = sympy.Eq(j, kin*(1-ps[hops[0]]))
    #final = sympy.Eq(j, kout*ps[hops[-1]])
    
    # analytic part of repeated substitutions
    #substitutes = sympy.solve(initial, ps[hops[0]])
    #population_expressions.append(substitutes[0])
    #assert(len(substitutes) == 1)
    counter = 0
    for start, stop in zip(hops[:-1], hops[1:]):
            counter += 1
            population_expressionsN.append(nsimplify(pplusN(start, stop, population_expressionsN[-1], population_expressionsD[-1], population_expressionsN[0])))# pi, numerator
            population_expressionsD.append(nsimplify(pplusD(start, stop, population_expressionsN[-2], population_expressionsD[-1], population_expressionsN[0])))# pi, denominator
 	    #print nsimplify(population_expressionsN[-1])/nsimplify(population_expressionsD[-1])
            #base = build_jji(start, stop)
            #t = base.subs(ps[start], substitutes[0])
            #substitutes = sympy.solve(t, ps[stop])
	    #print substitutes[0]
            #population_expressions.append(substitutes[0])
            #assert(len(substitutes) == 1)
    #t = final.subs(ps[hops[-1]], population_expressions[-1])
    t = sympy.Eq(j, kout*population_expressionsN[-1]/population_expressionsD[-1])
    #jsolve = sympy.solve(t.subs(replacements), j)
    #print "jsolve", jsolve
    #print "jsolve[0]", jsolve[0]
    
    assert(type(t) == sympy.Eq)
    assert(type(t.rhs) == sympy.Mul)
    assert(t.lhs == j)
    assert(type(t.rhs.args[1]) == sympy.Pow)
    assert(t.rhs.args[0] == kout)
    
    #build characteristic polynomial
    #print "Polynomial equation of J built!"
    #print "t", t
    jpol = sympy.poly(j*population_expressionsD[-1]-kout*population_expressionsN[-1], j)
    #print jpol
    base_expression = jpol.args[0].subs(replacements)
    coefficients = sympy.poly(base_expression).coeffs()
    roots = numpy.roots(coefficients)
    
    #calculate populations on remaining solutions
    population_expressions = []
    for upper, lower in zip(population_expressionsN, population_expressionsD):
        population_expressions.append(nsimplify(upper/lower))
    
    #roots = roots[numpy.where(roots > 0)]
    roots = roots[numpy.where(abs(numpy.imag(roots)) <1.0)]
    roots=numpy.real(roots)
    #print "roots are:", roots
    plausible_js = []
    solutions = []
    for root in roots:
        plausible = True
        for population_expression, pvar in zip(population_expressions, ps):
            base_expression = population_expression.subs(replacements).subs(j, root)
            if float(base_expression) < -1.0e-3 or float(base_expression) > 1.00100:
 	        #print "root is %4.2e, probability not physical:" %root, float(base_expression)
                plausible = False
        if plausible:
            plausible_js.append(root)
            #print >>fo, math.log10(k_out_val), math.log10(root)
            for population_expression, pvar in zip(population_expressions, ps):
                base_expression = population_expression.subs(replacements).subs(j, root)
                solutions.append(float(base_expression))
            solutions.append(root)
    assert(len(plausible_js) == 1)
    return solutions

#a = solve_flux([13391333, 356671335, 14905524], [681249, 75912061, 20311347]) 
#print a
