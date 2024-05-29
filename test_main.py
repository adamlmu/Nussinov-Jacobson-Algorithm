from nussinov import *
import json

def parentheses_cleaning(dot_parentheses):
    dot_parentheses = dot_parentheses.replace('[', '(')
    dot_parentheses = dot_parentheses.replace(']', ')')
    return dot_parentheses

def count_parentheses(dot_par):
    dot_par = parentheses_cleaning(dot_par)
    count = 0
    for i in dot_par:
        if i == '(':
            count += 1
    return count

if __name__ == '__main__':
    
    # define constant variables
    GAMMA = set([('A', 'U'), ('U', 'A'), ('C', 'G'), ('G', 'C'), ('G', 'U'), ('U', 'G')])
    ell = 1

    # Read in test data
    with open("tests.json", "r") as test_file:
        data = json.load(test_file)

    # Comparing score
    for database in data:
        for seq in data[database]:
            sequence = data[database][seq]['sequence']
            dot_parentheses = data[database][seq]['dot-parentheses']

            ans, table = nussinov_table(sequence, GAMMA, ell)
            score = count_parentheses(dot_parentheses)

            print(seq)
            print(sequence)
            print(dot_parentheses)
            print("Our score:", ans)
            print("True score:", score)
            print()

# RCSB Protein Data Bank
# RNase P Database
# Rfam Database

