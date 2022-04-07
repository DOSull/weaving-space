#!/usr/bin/env python
# coding: utf-8

# Some strand label parser stuff
# This is much nicer in python than in R

import re


def _parse_strand_label(s:str) -> list[str]:
    """Breaks a strand label specifiction in to a list of labels

    Args:
        s (str): see get_strand_ids() for details
    Returns:
        list[str]: list of strand labels.
    """    
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


def get_strand_ids(strands_spec: str) -> tuple[list[str]]:
    """Strands label string to split into strand labels.

    Args:
        strands_spec (str): string format "a|bc|(de)f" | separates strands in
            each direction and () designates combining labels into a single strand that will be sliced lengthwise. Example output:

                "a|bc|(de)f" -> (["a"], ["b", "c"], ["de", "f"])
            
            Superflous parentheses are removed, but no other error-checks are
            applied.

    Returns:
        tuple[str]: tuple of lists of labels for each set of strands.
    """    
    strand_ids = [_parse_strand_label(s) for s in strands_spec.split("|")]
    strand_ids = (strand_ids
                  if len(strand_ids) == 3
                  else strand_ids + [[""]])
    return tuple(strand_ids)

