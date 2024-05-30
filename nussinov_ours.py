
import numpy as np
from draw_rna import ipynb_draw

def energy_score(pair):
    pairs = {   ("G", "C"): 3, ("C", "G"): 3,
                ("A", "U"): 2, ("U", "A"): 2,
                ("G", "U"): 1, ("U", "G"): 1}
    return pairs.get(pair, 0)

def fill_matrix(seq):
    n = len(seq)
    m = np.zeros((n, n))
    for p in range(1, n):
        for i in range(n - p):
            j = i + p
            m[i][j] = calculate_ij_score(seq, m, i, j)
    return m

def calculate_ij_score(seq, m, i, j):
    if j - i >= 3:
        options = [m[i+1][j-1] + energy_score((seq[i], seq[j])),
                   m[i+1][j],
                   m[i][j-1]]
        for k in range(i+1, j):
            options.append(m[i][k] + m[k+1][j])

        return max(options)
    
    else:
        return 0

def traceback(m, seq, fold, i, j):
    pair = (seq[i], seq[j])

    if i < j:
        if m[i][j] == m[i][j-1]:
            traceback(m, seq, fold, i, j-1) 
        elif m[i][j] == m[i+1][j]:
            traceback(m, seq, fold, i+1, j)
        elif m[i][j] == m[i+1][j-1] + energy_score(pair):
            fold.append((i, j))
            traceback(m, seq, fold, i+1, j-1)
        else:
            for k in range(i+1, j):
                if m[i][j] == m[i][k] + m[k+1][j]:
                    traceback(m, seq, fold, i, k)
                    traceback(m, seq, fold, k+1, j)
                    break
    
    return fold

def generate_dot_bracket(n, fold):
    dot_bracket = ["." for i in range(n)]
    for start, end in fold:
        dot_bracket[start] = "("
        dot_bracket[end] = ")"

    return "".join(dot_bracket)

if __name__ == "__main__":
    seq = "GGUCGGCAGUACCAAGUCGCGAAAGCGAUGGCCUUGCAAAGGGUAUGGUAAUAAGCUGCCCAC"
    n = len(seq)

    m = fill_matrix(seq)
    fold = traceback(m, seq, [], 0, n - 1)
    dot_bracket_str = generate_dot_bracket(n, fold)

    print("RNA sequence:\n", seq)
    print("RNA fold in dot-bracket representation:\n", dot_bracket_str)
    ipynb_draw.draw_struct(seq, dot_bracket_str)
