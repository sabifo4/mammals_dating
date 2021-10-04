import ete3
from Bio import Phylo

BASE_TREE = 'FigTree.tree'  # NOTE: you need to name the backbone tree!

# the base 72sp reference tree
# t72sp = Phylo.read(BASE_TREE, 'nexus')
# taxa72sp = set([x.name for x in t72sp.get_terminals()])

# outgroup taxa
outgroup = {'ornithorhynchus_anatinus', 'zaglossus_bruijni', 'tachyglossus_aculeatus'}

# replacement rules:
# - each subtree can have one or more rules specifying which clades to copy
# - 'source' specifies taxa whose common ancestor will be copied from subtree
# - 'target' specifies taxa whose common ancestor will be replace in 72sp tree with 'source'
# - 'placeholder' if True (default is False) will treat 'source' as a placeholder tip to insert as sibling of 'target'
on_laurasiatheria_therest = {
    'Laurasiatheria_chiro_subt1': [
        {
            'source': ['rhinopoma_hardwickii', 'myonycteris_angolensis'],
            'target': ['pteropus_vampyrus']
        }
    ],
    'Laurasiatheria_chiro_subt2': [
        {
            'source': ['nycteris_javanica', 'laephotis_wintoni'],
            'target': ['myotis_lucifugus']
        }
    ],
    'Laurasiatheria_cetartiodactyla': [
        {
            'source': ['vicugna_pacos', 'sus_scrofa', 'tursiops_truncatus', 'bos_taurus', 'ovis_aries', 'capra_hircus'],
            'target': ['vicugna_pacos', 'capra_hircus']
        }
    ],
}

direct_on_72sp = {
    'Afrotheria': [
        {
            'source': ['procavia_capensis', 'loxodonta_africana', 'echinops_telfairi'],
            'target': ['procavia_capensis', 'loxodonta_africana', 'echinops_telfairi']
        }
    ],
    'Euarchonta': [
        {
            'source': ['otolemur_garnettii', 'homo_sapiens', 'cynocephalus_volans'],  # cynocephalus_volans
            'target': ['otolemur_garnettii', 'homo_sapiens']
        },
        {
            'source': ['tupaia_belangeri', 'ptilocercus_lowii'],
            'target': ['tupaia_belangeri']
        },
    ],
    'Lagomorpha': [
        {
            'source': ['oryctolagus_cuniculus', 'ochotona_princeps'],
            'target': ['oryctolagus_cuniculus', 'ochotona_princeps']
        }
    ],
    'Marsupialia': [
        {
            'source': ['sarcophilus_harrisii', 'macropus_eugenii', 'monodelphis_domestica', 'lestoros_inca'],
            'target': ['monodelphis_domestica', 'notamacropus_eugenii', 'sarcophilus_harrisii'],
        },
        {
            'source': ['tachyglossus_aculeatus', 'zaglossus_bruijni', 'ornithorhynchus_anatinus'],
            'target': ['ornithorhynchus_anatinus']
        }
    ],
    'Xenarthra': [
        {
            'source': ['choloepus_hoffmanni', 'dasypus_novemcinctus'],
            'target': ['choloepus_hoffmanni', 'dasypus_novemcinctus']
        }
    ],
    'Rodentia_squirrel': [
        {
            'source': ['eliomys_quercinus', 'urocitellus_parryii_ablusus'],
            'target': ['ictidomys_tridecemlineatus']
        }
    ],
    'Rodentia_ctenohystrica_3': [
        {
            'source': ['laonastes_aenigmamus', 'phyllomys_aff_dasythrix'],
            'target': ['fukomys_damarensis', 'octodon_degus']
        }
    ],
    'Rodentia_subtree1': [
        {
            'source': ['pedetes_capensis_surdaster', 'peromyscus_boylii_rowleyi'],
            'target': ['rattus_norvegicus', 'dipodomys_ordii']
        }
    ],
    'Rodentia_subtree2': [
        {
            'source': ['rheomys_raptor', 'praomys_petteri'],
            'target': ['rattus_norvegicus', 'mus_musculus']
        }
    ],
    # laurasiatheria subtrees should have been grafted together
    'All_Laurasiatheria': [
        {
            'source': ['aonyx_cinerea', 'solenodon_paradoxus'],
            'target': ['sorex_araneus', 'capra_hircus']
        }
    ]
}


# # add annotation indicating source of 95% interval (and, therefore, node age)
# for n in t72sp.find_elements():
#     if hasattr(n, 'comment') and n.comment is not None:
#         n.comment = n.comment[:-1] + ',from="72sp"]'

# uncomment if you want separate trees for each step - useful for debugging
# t72sp = Phylo.read(BASE_TREE, 'nexus')
# for n in t72sp.find_elements():
#     if hasattr(n, 'comment') and n.comment is not None:
#         n.comment = n.comment[:-1] + ',from="72sp"]'


def do_one_rule_for_one_tree(backbone_tree: ete3.Tree, subtree: ete3.Tree, graft_rule):
    if graft_rule.get('placeholder', False):
        from Bio.Phylo.BaseTree import Clade
        ca_backbone = backbone_tree.common_ancestor(graft_rule['target'])
        parent = backbone_tree.get_path(ca_backbone)[-2]
        sibling_index = 0 if parent.clades[0] == ca_backbone else 1
        ca_subtree = Clade(branch_length=1, name=graft_rule['source'][0])
        ca_subtree.comment = '[&from="placeholder"]'
        new_branch_length = parent.clades[sibling_index].branch_length / 2
        parent.clades[sibling_index].branch_length = parent.clades[sibling_index].branch_length + new_branch_length
        parent.clades[sibling_index].branch_length = new_branch_length
        new_clade = Clade(branch_length=new_branch_length, clades=[parent.clades[sibling_index], ca_subtree])
        new_clade.comment = '[&from="placeholder"]'
        parent.clades[sibling_index] = new_clade
        distance_backbone = backbone_tree.distance(new_clade, graft_rule['target'][0])
        ca_subtree.branch_length = distance_backbone
    else:
        ca_subtree = subtree.common_ancestor(graft_rule['source'])
        ca_backbone = backbone_tree.common_ancestor(graft_rule['target'])
        # check all taxa in 72sp are present in subtree
        taxa_ca_backbone = set([x.name for x in ca_backbone.get_terminals()])
        taxa_ca_subtree = set([x.name for x in ca_subtree.get_terminals()])
        common_taxa = taxa_ca_backbone.intersection(taxa_ca_subtree)
        # if not all the taxa from the 72sp tree are present in the common taxa
        if common_taxa != taxa_ca_backbone:
            print(f'Backbone tree has {taxa_ca_backbone.difference(common_taxa)} but at not in subtree')
        # get the distance from tip to common ancestor in 72sp and subtree
        distance_backbone = backbone_tree.distance(ca_backbone, next(iter(common_taxa)))
        distance_subtree = subtree.distance(ca_subtree, next(iter(common_taxa)))

        # new stem length for 72sp tree
        new_stem_length = ca_backbone.branch_length + (distance_backbone - distance_subtree)

        print(f'MRCA age: {distance_backbone} -> {distance_subtree}')
        print(f'Adjusted branch length: {ca_backbone.branch_length} -> {new_stem_length}')
        print(f'Replacing {len(ca_backbone.get_terminals())} taxa with {len(ca_subtree.get_terminals())}')

        ca_backbone.branch_length = new_stem_length
        ca_backbone.comment = ca_subtree.comment
        ca_backbone.clades = ca_subtree.clades
        ca_backbone.name = None

    return backbone_tree

    # prefix names
    # for c in ca_72sp.get_terminals():
    #     c.name = f'{index}_{rule_index}_{c.name}'

    # check distances
    # for rule in rules_so_far:
    #     ca_72sp = t72sp.common_ancestor(rule)
    #     print('Node of MRCA', rule, ':', t72sp.distance(ca_72sp, rule[0]))
    #
    # print(f'Write {index}_{rule_index}')
    # Phylo.write(t72sp, f'5000sp_{index}_{rule_index}.tree', 'nexus')


def do_all_rules_on_tree(backbone_tree_file, prune_graft_rules, gather=False):
    if gather:
        filename = backbone_tree_file
        backbone_tree_file = Phylo.read(f'{backbone_tree_file}.tree', 'nexus')
        for n in backbone_tree_file.find_elements():
            if hasattr(n, 'comment') and n.comment is not None and "from=" not in n.comment:
                n.comment = n.comment[:-1] + f',from="{filename}"]'

    for tree_file, graft_rules in prune_graft_rules.items():
        backbone_tree_file = do_all_rules_for_one_tree(backbone_tree_file, tree_file, graft_rules, gather)
        print()

    return backbone_tree_file


def do_all_rules_for_one_tree(backbone_tree_file, tree_file, graft_rules, gather):
    subtree = Phylo.read(f'{tree_file}.tree', 'nexus')
    for n in subtree.find_elements():
        if hasattr(n, 'comment') and n.comment is not None and "from=" not in n.comment:
            n.comment = n.comment[:-1] + f',from="{tree_file}"]'

    # fix tilda names in subtree
    for terminal in subtree.get_terminals():
        if terminal.name.endswith('~'):
            terminal.name = terminal.name[:-1]

    total_subtree_taxa = len(subtree.get_terminals()) - len(outgroup)
    print(f'{"-" * 10} Subtree {tree_file}, {total_subtree_taxa} taxa {"-" * 10}')

    if isinstance(backbone_tree_file, str):
        backbone_tree_file = Phylo.read(f'{backbone_tree_file}.tree', 'nexus')

    rule_index = 0
    for graft_rule in graft_rules:
        if isinstance(backbone_tree_file, str):
            backbone_tree_file = Phylo.read(f'{backbone_tree_file}.tree', 'nexus')
        rule_index += 1
        print(f'- Applying rule {rule_index}')
        grafted_tree = do_one_rule_for_one_tree(backbone_tree_file, subtree, graft_rule)
        out_tree = f'Working_{tree_file}_{rule_index}.tree'
        Phylo.write(grafted_tree, out_tree, 'nexus')
        print()

    return backbone_tree_file


all_laurasiatheria = do_all_rules_on_tree('Laurasiatheria_therest', on_laurasiatheria_therest, gather=True)
Phylo.write(all_laurasiatheria, 'All_Laurasiatheria.tree', 'nexus')

all_mammals = do_all_rules_on_tree('Figtree', direct_on_72sp, gather=True)
Phylo.write(all_mammals, '4705sp.tree', 'nexus')

# taxa_in_1 = {t.name for t in Phylo.read('Rodentia_squirrel.tree', 'nexus').get_terminals()}
# taxa_in_2 = {t.name for t in Phylo.read('Rodentia_ctenohystrica_3.tree', 'nexus').get_terminals()}
# taxa_in_3 = {t.name for t in Phylo.read('Rodentia_subtree1.tree', 'nexus').get_terminals()}
# taxa_in_4 = {t.name for t in Phylo.read('Rodentia_subtree2.tree', 'nexus').get_terminals()}

# index = 0
# for tree_file, graft_rules in prune_graft_rules.items():
#     index += 1
#     subtree = Phylo.read(tree_file, 'nexus')
#     for n in subtree.find_elements():
#         if hasattr(n, 'comment') and n.comment is not None:
#             n.comment = n.comment[:-1] + f',from="subtree"]'
#
#     print(f'{index}. {tree_file}')
#
#     # fix tilda names in subtree
#     for t in subtree.get_terminals():
#         if t.name.endswith('~'):
#             t.name = t.name[:-1]
#
#     total_subtree_taxa = len(subtree.get_terminals()) - len(outgroup)
#
#     grafted_taxa_count = 0
#     rule_index = 0
#     rules_so_far = list()
#     for graft_rule in graft_rules:
#         do_graft()
#
#     print(f'{grafted_taxa_count}/{total_subtree_taxa} grafted')
#     print('-' * 80)

# afrotheria - monophyletic replace
# lagomorpha - monophyletic replace
# chiroptera - monophyletic replace
# marsupialia - monophyletic replace
# xenarthra - monophyletic replace
# euarchonta - separate (i) scandentia (ii) primates+dermoptera
# squirrel - single taxon replaced
# L.cetartiodactyla - monophyletic replace
# L.therest - separate (i) Lipotyphla (ii) Perissodactyla replace
# carnivora - (what to do about pholidota)
# outgroup


