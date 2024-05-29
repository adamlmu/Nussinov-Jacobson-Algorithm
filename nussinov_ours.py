import numpy as np

def solv_matrix(m, seq):
    min_loop = 3
    n = len(seq)
    pairs = {("A", "U"): 2, ("U", "A"): 2,
             ("G", "C"): 3, ("C", "G"): 3,
             ("G", "U"): 1, ("U", "G"): 1}
    
    for p in range(1, n):
            # in pass p, (i,j) moves diagonally down and to the right
            # p elements above the main diagonal.
            for i in range(n - p):
                # Fill M[i, j] with max score for interval [i,j] incl.
                j = i + p
                
                if(j-i >= min_loop):
                     pair = (seq(i),seq(j))
                     diagonal = m[i+1][j-1] + pairs.get(pair)
                     down = m[i+1][j]
                     left = m[i][j-1]
                     bifurcation = max([m[i][k] + m[k + 1][j] for k in range(i, j)])

                     m[i][j] = max(diagonal,down,left,bifurcation)

                else:
                     m[i][j] = 0
    return m
