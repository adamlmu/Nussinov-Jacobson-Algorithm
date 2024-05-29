import numpy as np
from draw_rna import ipynb_draw

def couple(pair):
    """
    Return the energy of the pair if RNA nucleotides are Watson-Crick or wobble base pairs
    """
    pairs = {("A", "U"): 2, ("U", "A"): 2,
             ("G", "C"): 3, ("C", "G"): 3,
             ("G", "U"): 1, ("U", "G"): 1}
    return pairs.get(pair, 0)

def fill(nm, rna):
    """
    Fill the matrix as per the Nussinov algorithm
    """
    minimal_loop_length = 3

    for k in range(1, len(rna)):
        for i in range(len(rna) - k):
            j = i + k

            if j - i >= minimal_loop_length:
                down = nm[i + 1][j] # 1st rule
                left = nm[i][j - 1] # 2nd rule
                diag = nm[i + 1][j - 1] + couple((rna[i], rna[j])) # 3rd rule

                rc = max([nm[i][t] + nm[t + 1][j] for t in range(i, j)]) # 4th rule

                nm[i][j] = max(down, left, diag, rc) # max of all
            
            else:
                nm[i][j] = 0

    return nm

def traceback(nm, rna, fold, i, j):
    """
    Traceback through complete Nussinov matrix to find optimal RNA secondary structure solution through max base-pairs
    """
    if i < j:
        if nm[i][j] == nm[i + 1][j]: # 1st rule
            traceback(nm, rna, fold, i + 1, j)
        elif nm[i][j] == nm[i][j - 1]: # 2nd rule
            traceback(nm, rna, fold, i, j - 1)
        elif nm[i][j] == nm[i + 1][j - 1] + couple((rna[i], rna[j])): # 3rd rule
            fold.append((i, j))
            traceback(nm, rna, fold, i + 1, j - 1)
        else:
            for k in range(i + 1, j):
                if nm[i][j] == nm[i, k] + nm[k + 1][j]: # 4th rule
                    traceback(nm, rna, fold, i, k)
                    traceback(nm, rna, fold, k + 1, j)
                    break

    return fold

def dot_write(rna, fold):
    dot = ["." for i in range(len(rna))]
    for s in fold:
        dot[min(s)] = "("
        dot[max(s)] = ")"
    return "".join(dot)

def init_matrix(rna):
    M = len(rna)

    # init matrix
    nm = np.zeros((M, M))

    return nm

if __name__ == "__main__":
    rna = "GGUCGGCAGUACCAAGUCGCGAAAGCGAUGGCCUUGCAAAGGGUAUGGUAAUAAGCUGCCCAC" # our RNA sequence
    nm = init_matrix(rna)

    nm = fill(nm, rna)

    fold = traceback(nm, rna, [], 0, len(rna) - 1)
    dot_bracket = dot_write(rna, fold)

    print("RNA Sequence:", rna)
    print("Dot-Bracket Notation:", dot_bracket)

    ipynb_draw.draw_struct(rna, dot_bracket)
