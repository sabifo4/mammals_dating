filename = '4705sp.tree'

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade


def remove_comments(clade: Clade):
    for c in clade.clades:
        remove_comments(c)
    delattr(clade, 'comment')


# for the posterior mean i.e. the existing branch lengths
all5000sp = Phylo.read(filename, format='nexus')
remove_comments(all5000sp.clade)
Phylo.write(all5000sp, '4705sp_mean.nwk', format='newick')

import re
comment_match = re.compile(r'{([.\dEe\-]+),\s?([\d.Ee\-]+)}')


def update_branch_length_from_comment(clade: Clade, low):
    comment = getattr(clade, 'comment')
    print(comment)
    c025, c975 = comment_match.findall(comment)[0]
    c025, c975 = float(c025), float(c975)
    if low:
        desired_node_age = c025
    else:
        desired_node_age = c975
    for child in clade.clades:
        if child.is_terminal():
            # update the height based on comment
            child.branch_length = desired_node_age
            print(child.name, 'is now', child.branch_length)
        else:
            update_branch_length_from_comment(child, low)
            distance_to_terminal_through_child = all5000sp.distance(child, child.get_terminals()[0])
            new_branch_length = desired_node_age - distance_to_terminal_through_child
            print(distance_to_terminal_through_child, desired_node_age, '->', child.branch_length, new_branch_length)
            # all descendants should have been done already, now only need to change this clade's branch length
            child.branch_length = new_branch_length


# for the 2.5% ci
all5000sp = Phylo.read(filename, format='nexus')
update_branch_length_from_comment(all5000sp.clade, True)
Phylo.write(all5000sp, '4705sp_ci025.nexus', format='nexus')
#---------------------------------------------------------------#
# SAC-211110 | Read file and, in case "\[" or "\]"
# are found, then remove them.
# E.g., [\[&95%={0.0694354, 0.280382},from="Marsupialia-time"\]]
# This seems to happen when using the WLS
with open('4705sp_ci025.nexus', 'r') as file :
  inp_file = file.read()

inp_file = inp_file.replace('\[', '')
inp_file = inp_file.replace('\]', '')

with open('4705sp_ci025.nexus', 'w') as file:
  file.write(inp_file)
#---------------------------------------------------------------#
remove_comments(all5000sp.clade)
Phylo.write(all5000sp, '4705sp_ci025.nwk', format='newick')

# for the 97.5% ci
all5000sp = Phylo.read(filename, format='nexus')
update_branch_length_from_comment(all5000sp.clade, False)
Phylo.write(all5000sp, '4705sp_ci975.nexus', format='nexus')
#---------------------------------------------------------------#
# SAC-211110 | Read file and, in case "\[" or "\]"
# are found, then remove them.
# E.g., [\[&95%={0.0694354, 0.280382},from="Marsupialia-time"\]]
# This seems to happen when using the WLS
with open('4705sp_ci975.nexus', 'r') as file :
  inp_file = file.read()

inp_file = inp_file.replace('\[', '')
inp_file = inp_file.replace('\]', '')

with open('4705sp_ci975.nexus', 'w') as file:
  file.write(inp_file)
#---------------------------------------------------------------#
remove_comments(all5000sp.clade)
Phylo.write(all5000sp, '4705sp_ci975.nwk', format='newick')
