#!/usr/bin/env python
# coding: utf-8
#
# Bernardo M. Rocha, Fev, 2014
#

import sys

if __name__ == "__main__":

    f = open(sys.argv[1])
    l = f.readline()
    l = l.split()
    l = l[1:]
    n = len(l)
    for i in range(n):
        print('values[%d] = %e;' % (i, float(l[i])) )

