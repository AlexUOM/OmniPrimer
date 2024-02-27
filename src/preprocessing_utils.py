"""Collection of preprocessing functions used to determine
   the flanking sequences around each exon"""

import itertools


def left_flanking(exon, intron, pos):
    """
    This function is used to process and prepare a 300 nt flanking sequence
    upstream the exon. The flanking sequence consists of any 300 nt found upstream
    the exon, regardless if these are intronic or exonic sequences.

    The default lenght is 300 nt, which is the reccommended size. This can be modified
    if necessary.
    """

    flanking_seq = [intron[pos]]
    while len("".join(flanking_seq)) <= 300:
        try:
            flanking_seq.append(exon[pos - 1])
            flanking_seq.append(intron[pos - 1])
            pos -= 1
        except:
            break
    flanking_seq = "".join(reversed(flanking_seq))
    flanking_seq = flanking_seq[-300:]
    return flanking_seq


def right_flanking(exon, intron, pos):
    """
    This function is used to process and prepare a 300 nt flanking sequence
    downstream the exon. The flanking sequence consists of any 300 nt found downstream
    the exon, regardless if these are intronic or exonic sequences.

    The default lenght is 300 nt, which is the reccommended size. This can be modified
    if necessary.
    """
    flanking_seq = [intron[pos + 1]]
    while len("".join(flanking_seq)) <= 300:
        try:
            flanking_seq.append(exon[pos + 1])
            flanking_seq.append(intron[pos + 2])
            pos += 1
        except:
            break
    flanking_seq = "".join(flanking_seq)
    flanking_seq = flanking_seq[:300]
    return flanking_seq


def letter_permutations():
    """
    This function is used to generate all possible permutations
    of alphabet letters.
    """

    for size in itertools.count(1):
        for s in itertools.product("abcdefghijklmnopqrstuvwxyz", repeat=size):
            yield "".join(s)
