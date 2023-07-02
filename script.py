import dendropy

# Given a tree determine taxa list [.]
# Prune tree to only contain said taxa [.]
# Run AU tests [ ]

# input_paths = [["constrained_data.fas", "gene_tree.txt", "concat_tree.txt", "constrained_tree.txt"],["constrained_data.fas", "gene_tree_2.txt", "concat_tree.txt", "constrained_tree.txt"],["constrained_data.fas", "gene_tree.txt", "gene_tree_2.txt", "concat_tree.txt", "constrained_tree.txt"]]
input_paths = [["AA_dataset.fas", "AA_gene_tree.txt.contree", "AA_concat.txt.contree", "AA_constrained.txt.contree"]]
test_names = ["AidanTrees"]

def concat_trees(trees):
    out = ""
    for tree in trees.values():
        out += tree.as_string(schema='newick') + "\n"
    return out


# Read Newick format tree file
def read_newick(filepath):
    tree = dendropy.Tree.get(path=filepath, schema='newick', preserve_underscores=True)
    print(f"Loaded tree: {filepath}")
    # print(tree.as_string(schema='newick'))
    # print(tree.as_ascii_plot())
    return tree

def find_common_taxa(taxa_sets):
    # Get the set of items present in the first list
    unique_items = set(taxa_sets[0])

    # Iterate over the remaining lists
    for taxa_set in taxa_sets[1:]:
        # Update the set of unique items by taking the intersection with each list
        unique_items.intersection_update(taxa_set)

    # Convert the set back to a list
    unique_items_list = list(unique_items)

    return unique_items_list

def get_labels_from_tree(tree):
    return[str(item.label) for item in tree.taxon_namespace]

# Input tree files
# treefiles = ["gene_tree.txt", "gene_tree_2.txt", "concat_tree.txt", "constrained_tree.txt"]
def find_pruned_trees(treefiles):
    tree_labels = {}
    trees = {}
    pruned_trees = {}

    alignment_file = treefiles[0]
    treefiles.pop(0)

    for filepath in treefiles:
        tree = read_newick(filepath)
        tree_labels[filepath] = get_labels_from_tree(tree)
        trees[filepath] = tree

    # print(tree_labels.values())
    common_taxa = find_common_taxa(list(tree_labels.values()))
    print("Found taxa in common:")
    # print(common_taxa)

    for filepath in treefiles:
        
        # find taxa to remove
        taxa_to_be_pruned = [item for item in tree_labels[filepath] if item not in common_taxa]

        # Make pruned version
        pruned_trees[filepath] = dendropy.Tree(trees[filepath])

        # Remove taxa from copy
        pruned_trees[filepath].prune_taxa_with_labels(taxa_to_be_pruned)

        # Save copy to file
        pruned_trees[filepath].write(path=filepath + "_pruned.tre", schema='newick')
        print("Pruned " + filepath)
        # print(pruned_trees[filepath].as_ascii_plot())


    print(common_taxa)

    print("Generating alignment file")
    pruned_alignment = ""
    with open(alignment_file, "r") as f:  
        for line in f:
            if line.startswith(">"):
                label = line.lstrip(">")
                label = label.rstrip("\n")
                if label in common_taxa:
                    pruned_alignment += ">" + label
                else:
                    label = ""
            
            elif label:
                pruned_alignment += "\n" + line

    with open(alignment_file + "_pruned.fas", "w") as f:
        f.write(pruned_alignment)

    print("Generated alignment file")

    return pruned_trees


# START SCRIPT


for i in range(len(input_paths)):
    with open(test_names[i] + "_pruned.treels", 'w') as f:

        print(">>> Finding Pruned Trees for " + str(", ".join(input_paths[i])))
        pruned_trees = find_pruned_trees(input_paths[i])
        f.write(concat_trees(pruned_trees))

print("Script complete")








# # Load gene tree
# gene_tree = read_newick("gene_tree.txt")
# gene_tree_taxa = get_taxa_labels(gene_tree)
# print("Gene Tree")
# print(gene_tree.as_ascii_plot())
# print("Gene tree labels")
# print(gene_tree.as_string(schema='newick'))
# print(gene_tree_taxa)

# # Load concat tree
# concat_tree = read_newick("concat_tree.txt")
# concat_tree_taxa = get_taxa_labels(concat_tree)
# print("Concat tree labels")
# print(concat_tree.as_string(schema='newick'))
# print(concat_tree_taxa)

# # Create a concat copy 'Pruned'
# print("Concat Tree (Prune copy)")
# pruned_tree = dendropy.Tree(concat_tree)
# pruned_taxa = pruned_tree.taxon_namespace
# print(pruned_tree.as_ascii_plot())

# # Determine labels to remove
# # gene_tree_taxa_labels = [item.label for item in gene_tree_taxa]
# # concat_tree_taxa_labels = [item.label for item in concat_tree_taxa]

# # taxa_to_be_pruned = [item for item in pruned_taxa if item.label not in gene_tree_taxa_labels]
# # print("Labels to remove...")
# # prune_tree_taxa_labels = [str(item.label) for item in taxa_to_be_pruned]
# # print(prune_tree_taxa_labels)

# pruned_tree.prune_taxa(taxa_to_be_pruned) # Prune tree

# # Print pruned tree to screen
# print("Pruned Tree:")
# print(get_taxa_labels(pruned_tree))
# print(pruned_tree.as_string(schema='newick'))
# print(pruned_tree.as_ascii_plot())

# # Write to file
# pruned_tree.write(path='pruned_tree.nwk', schema='newick')

# print("Script Complete")