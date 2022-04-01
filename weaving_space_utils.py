#!/usr/bin/env python
# coding: utf-8

# Some strand label parser stuff
# This is much nicer in python than in R

import re

def parse_labels(ids):
  output = ids.split("|")
  if len(output) == 2:
    output.append("-")
  return(output)
    

def parse_strand_label(s):
  clean_s = re.sub("[(]+", "(", re.sub("[)]+", ")", s))
  is_combo = False
  output = []
  current = ""
  for c in clean_s:
    if is_combo:
      if c == ")":
        output.append(current)
        current = ""
        is_combo = False
      else:
        current = current + c
    else:
      if c == "(":
        is_combo = True
      else:
        output.append(c)
  return output

def get_strand_ids(strands_spec):
  return [parse_strand_label(x) for x in parse_labels(strands_spec)]

