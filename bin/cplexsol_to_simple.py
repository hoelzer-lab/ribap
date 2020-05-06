#!/usr/bin/env python3
# -*- coding: utf-8
import sys
import argparse

__author__ = 'pedro'
import xml.etree.ElementTree as et
import re

#x_A569t_B768h
EDGE_PATTERN = re.compile('^x_A([0-9]+)t_B([0-9]+)([ht])')

def readVariables(filename, gap=0, single=False):

    sol_list = list()
    objective = list()
    #print >> sys.stderr, "Opening ..."
    tree = et.parse(filename)
    root = tree.getroot()
    #print >> sys.stderr, "Parsing ..."
    best_obj = None
    solutions = root.iterfind('./CPLEXSolution') if root.tag == 'CPLEXSolutions' else [root]
    for sol in solutions:
        obj = float(sol[0].attrib['objectiveValue'])
        if best_obj is None:
            best_obj = obj
        # gap filter
        if gap>0:
            if abs(best_obj-obj)>gap*best_obj: continue

        objective.append(obj)
        res = {}
        sol_list.append(res)
        for var in sol.iterfind('./variables/variable'):

            name = var.get('name')
            value = int(round(float(var.get('value'))))
            if abs(value) < 0.01: continue
            res[name] = value
        if single:
            break
    return (sol_list, objective)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Outputs the O.F. value and non-zero vars from a ILP solution files.')
    parser.add_argument('-gap', default=0, type=float, help="Relative gap for the solution to be included, in relation to the best solution.")
    parser.add_argument('--single', help="Just one solution is output. Overrides gap")
    parser.add_argument("file", help="Gene similarity file(s)", nargs='+')
    param = parser.parse_args()
    for filename in param.file:
        #print >> sys.stderr, filename
        # if os.path.isfile("%s.simple.1" % filename): continue
        sol_list, objective = readVariables(filename, gap=param.gap, single=param.single)
       # print >> sys.stderr, "%d solutions. Writing ..." % len(sol_list)

        idx = 1
        for res,obj in zip(sol_list, objective):
            f = open("%s.simple.%d" % (filename,idx), 'w')
            idx += 1
            f.write("%f\n" % obj)
            for g1 in sorted(res):
                f.write("%s\t%d\n" % (g1, res[g1]))
            f.close()


