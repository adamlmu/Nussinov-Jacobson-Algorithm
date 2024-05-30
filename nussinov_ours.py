import numpy as np
import sys
from draw_rna import ipynb_draw

def couple(pair):
    pairs = {("G", "C"): 3, ("C", "G"): 3, ("A", "U"): 2, ("U", "A"): 2, ("G", "U"): 1, ("U", "G"): 1}
    return pairs.get(pair, 0)

def solv_matrix(seq):
    n = len(seq)
    m = np.zeros((n,n))
    for p in range(1, n):
        for i in range(n - p):
            j = i + p
            if(j - i >= 3):  # min length of a loop
                pair = (seq[i], seq[j])
                diagonal = m[i+1][j-1] + couple(pair)
                down = m[i+1][j]
                left = m[i][j-1]
                bifurcation = max([m[i][k] + m[k + 1][j] for k in range(i, j)])
                m[i][j] = max(diagonal, down, left, bifurcation)
            else:
                m[i][j] = 0
    return m

def traceback(m, seq, fold, i, j):
    pair = (seq[i], seq[j])
    if i < j:
        if m[i][j] == m[i][j-1]:
            traceback(m, seq, fold, i, j-1) 
        elif m[i][j] == m[i+1][j]:
            traceback(m, seq, fold, i+1, j)
        elif m[i][j] == m[i+1][j-1] + couple(pair):
            fold.append((i, j))
            traceback(m, seq, fold, i+1, j-1)
        else:
            for k in range(i+1, j):
                if m[i][j] == m[i][k] + m[k+1][j]:
                    traceback(m, seq, fold, i, k)
                    traceback(m, seq, fold, k+1, j)
                    break
    return fold

def dot_bracket_str(n, fold):
    dot = ["." for i in range(n)]
    for i, j in fold:
        dot[i] = "("
        dot[j] = ")"
    return "".join(dot)

if __name__ == "__main__":
    seq = "GGUCGGCAGUACCAAGUCGCGAAAGCGAUGGCCUUGCAAAGGGUAUGGUAAUAAGCUGCCCAC"
    n = len(seq)
    m = solv_matrix(seq)
    fold = traceback(m, seq, [], 0, n-1)
    dot_bracket = dot_bracket_str(n, fold)
    print("RNA seqence:\n", seq)
    print("RNA fold in dot-bracket representation:\n", dot_bracket)
    ipynb_draw.draw_struct(seq, dot_bracket)