#!/usr/bin/env python

# bin/Binsanity
def get_cov_data(h):
    """Makes a dictionary of coverage values for affinity propagation"""
    all_cov_data = {}
    for line in open(str(h), "r"):
        line = line.rstrip()
        cov_data = line.split()
        all_cov_data[cov_data[0]] = cov_data
    return all_cov_data

