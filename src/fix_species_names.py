#!/usr/bin/env python
"""
1. names.txt
2. RAxML_bestTree.BS_ML_GTRCAT
3. alignment.phylip
4. lineage.txt
"""

import glob, os, shutil, sys
import ete3

def get_names(filepath, quiet=True):
    old_names = { }
    # Read in all the names
    with open(filepath, 'r') as f:
        for line in f.readlines():
            names = line.strip()
            if ' ' not in names:
                old_names[names] = names
            else:
                old_names[names.split(' ')[0]] = names
    if not quiet:
        print('%d names found' % len(old_names))
    return old_names


def fix_names(species, filepath):
    old_names = get_names(filepath, quiet=False)
    changes = 0

    # Loop through each rename and handle as necessary
    for from_name, to_name in species.items():
        if from_name in old_names.keys():
            print('Rename %s to %s' % (from_name, to_name))
            if to_name in old_names.keys():
                print('Error: %s already exists!' % to_name)
                sys.exit(1)
            else:
                if ' ' not in old_names[from_name]:
                    old_names[to_name] = to_name
                    del old_names[from_name]
                    changes += 1
                else:
                    print('Error: compound taxa rename %s' % from_name)
                    sys.exit(1)
        else:
            print('Error: %s not found!' % from_name)

    print('%d names after %d replacements' % (len(old_names), changes))
    with open(os.path.basename(filepath), 'w') as f:
        for key in sorted(old_names):
            f.write('%s\n' % old_names[key])


def fix_alignment(species, filepath):
    """ Writes out as we go """
    # SAC-200330> Needed to add a new var for the path to output file to which lines with fixed species
    #             name are written because I cannot have two files being read and written simultaneously
    outname_raw = os.path.basename(filepath).split(".")
    outname = outname_raw[0] + "_out.phylip"
    outfile = os.path.basename(outname)	
    #outfile = os.path.basename(filepath) # SAC-200330> Commented. See explanation L.57-58 
    changes = 0
    with open(outfile, 'w') as f:
        for line in open(filepath, 'r'):
            parts = line.split(' ')		
            if parts[0] in species:
                #print('Renaming %s' % parts[0]) # SAC-200330> Commented. I want to see what it renames it to
                print('Renaming %s to %s' % (parts[0], species[parts[0]])) # SAC-200330> Update for L.68
                changes += 1
                f.write('%s %s' % (species[parts[0]], ''.join(parts[1:])))
            else:
                f.write('%s' % line)
    print('%d changes in alignment' % changes)


def fix_tree(species, filepath):
    tree = ete3.Tree(filepath)
    leaves = tree.get_leaves()
    print('%d leaves found in tree (includes outgroup)' % len(leaves))
    done = set()
    for leaf in leaves:
        if leaf.name in species:
            if leaf.name in done:
                print('Error! Already matched %s' % leaf.name)
                sys.exit(1)
            done.add(leaf.name)
            leaf.name = species[leaf.name]
    if len(species.keys()) != len(done):
        print('Error! Count mismatch. Only found %s' % done)
        sys.exit(1)
    print('%d taxa names changed' % len(done))
    tree.write(outfile=os.path.basename(filepath), format=5)


def fix_lineage(species, filepath):
    lineage = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            parts = line.split('|')
            lineage[parts[1]] = parts[0]

    print('%d lineage entries' % len(lineage))

    changes = 0
    for from_name, to_name in species.items():
        if from_name in lineage:
            if to_name in lineage:
                print('Error: %s already exists' % to_name)
                sys.exit(1)
            else:
                changes += 1
                lineage[to_name] = lineage[from_name]
                del lineage[from_name]
        else:
            print('Error: %s not found in lineage' % from_name)
            sys.exit(1)

    print ('%d taxa renamed' % changes)

    backup_files = glob.glob('lineage.txt.*.bak')
    if len(backup_files):
        next_backup = max([int(x.split('.')[2]) for x in backup_files]) + 1
    else:
        next_backup = 1

    shutil.copyfile(filepath, '%s.%d.bak' % (os.path.basename(filepath), next_backup))

    current_names = get_names('names.txt')
    expanded_names = {}
    for key, value in current_names.items():
        if ' ' in value:
            for v in value.split(' '):
                expanded_names[v] = value
    current_names.update(expanded_names)

    new_lineage_count = 0
    with open(os.path.basename(filepath), 'w') as f:
        for name in sorted(current_names.keys()):
            if name in lineage:
                f.write('%s|%s\n' % (lineage[name], name))
                new_lineage_count += 1
            else:
                print('Warning: %s not found' % name)
    print('Wrote %d lineage lines' % new_lineage_count)


def print_heading(word):
    print('\n%s\n%s\n' % (word, '-' * len(word)))


def run(rename_filename, tree_filepath, alignment_filepath, names_filepath, lineage_filepath):
    # load species renames
    with open(rename_filename, 'r') as f:
        lines = (l.split(',') for l in f.readlines())
        species = { old_name : new_name.strip() for old_name, new_name in lines }
        print("%d species to rename" % len(species))

    print_heading('names.txt')
    fix_names(species, names_filepath)

    print_heading('alignment')
    fix_alignment(species, alignment_filepath)

    print_heading('tree')
    fix_tree(species, tree_filepath)

    print_heading('lineage')
    fix_lineage(species, lineage_filepath)


if __name__ == '__main__':
    run(*sys.argv[1:])

