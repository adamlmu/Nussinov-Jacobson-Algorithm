import os, datetime, argparse
from nussinov import *
import pdb

def visualize_reconstruction(seq, table, record, ell=0, output_name='out.png'):
    '''
        Creating backtrace from the given DP table and reconstruction

        Reference for plotting:
            https://matplotlib.org/3.1.1/gallery/images_contours_and_fields/image_annotated_heatmap.html
    '''
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np

    pairs = set(record)
    left_p = set([r[0] for r in record])
    right_p = set([r[1] for r in record])
    n = len(seq)
    visited = np.zeros((n, n))
    arrows = [] # quadruples of (x, y, dx, dy)
    v_mark = 1

    stack = [(0,n-1)]
    visited[0, n-1] = v_mark
    while len(stack) > 0:
        v_mark += 1
        i,j = stack.pop()
        #  print(i, j) # color this
        if i + ell >= j:
            continue
        elif (i, j) in pairs:
            # arrow to (i+1, j-1)
            arrows.append((i, j, 1, -1))
            visited[i+1, j-1] = v_mark
            stack.append((i+1, j-1))
        elif i in left_p and j in right_p:
            for k in range(i+1+ell, j):
                if (i, k) in pairs:
                    # draw arrow to (i, k)
                    arrows.append((i, j, 0, k-j))
                    visited[i, k] = v_mark
                    stack.append((i, k))
                    # arrow to (k+1, j)
                    arrows.append((i, j, k+1-i, 0))
                    visited[k+1, j] = v_mark
                    stack.append((k+1, j))
                    break
        elif i in left_p:
            # arrow to (i, j-1)
            arrows.append((i, j, 0, -1))
            visited[i, j-1] = v_mark
            stack.append((i, j-1))
        else:
            # arrow to (i+1, j)
            arrows.append((i, j, 1, 0))
            visited[i+1, j] = v_mark
            stack.append((i+1, j))

    rep = to_dot_bracket(seq, record)
    ax = plt.gca()
    # revert the color for visited entries
    visited[visited > 0] = visited.max() + 1 - visited[visited>0]
    im = ax.imshow(visited, cmap='YlGn')
    threshold = im.norm(im.get_array().max())/2.
    table = np.array(table)
    ax.set_xticks(np.arange(n))
    ax.set_xticklabels([rep[i] + '\n' + seq[i] + str(i+1) for i in range(n)])
    ax.set_yticks(np.arange(n))
    ax.set_yticklabels([seq[i] + str(i+1) for i in range(n)])
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
    plt.setp(ax.get_xticklabels(), rotation=0, ha="right",
             rotation_mode="anchor")
    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)
    ax.set_xticks(np.arange(n+1)-.5, minor=True)
    ax.set_yticks(np.arange(n+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3, zorder=0)
    ax.tick_params(which="minor", bottom=False, left=False)

    for i in range(n):
        for j in range(n):
            text = ax.text(j, i, table[i, j],
                           ha='center', va='center',
                           color='w' if im.norm(visited[i, j])>threshold else 'black')

    for (x, y, dx, dy) in arrows:
        if dy > 0:
            dy -= 0.5
        elif dy < 0:
            dy += 0.5
        if dx > 0:
            dx -= 0.5
        elif dx < 0:
            dx += 0.5
        ax.arrow(y, x, dy, dx, width=0.01, head_width=0.2, alpha=0.5,
                 color='r', zorder=2)

    plt.savefig(output_name)
    plt.close()

if __name__ == '__main__':
    # read input
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', type=str,
                        required=True,
                        help='Input string'
                       )
    args = parser.parse_args()
    string = args.input
    from vis_config import PARAMS as param
    GAMMA = set([])
    for (a, b) in param['pairing']:
        GAMMA.add((a, b))
        GAMMA.add((b, a))
    print('Input: {}'.format(string))
    ell = param['hairpin']
    # setup output
    rule_str = '_'.join([a + '-' + b for (a,b) in param['pairing']]) + '_hairpin_{}'.format(ell)
    if param['output_dir'][-1] != '/':
        param['output_dir'] += '/'
    out_dir = param['output_dir'] + '{}_{}/'.format(string, rule_str)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    ans, table = nussinov_table(string, GAMMA, ell)
    ## enumerate solutions
    all_solns = reconstruction_all(string, table, GAMMA, ell=ell)
    print('Found {} solutions with score {}'.format(len(all_solns), ans))
    for i, record in enumerate(all_solns):
        if param['max_num_soln'] is not None and i >= param['max_num_soln']:
            break
        output = to_dot_bracket(string, record)
        output_name = os.path.join(out_dir, 'traceback_{}.png'.format(output))
        visualize_reconstruction(string, table, record, ell=ell, output_name=output_name)
        print('Solution {}: {}'.format(i+1, output))
