import dendropy

#################
### FUNCTIONS ###
#################

def concat_trees(trees):
    """
    Collect all the pruned trees together into a single object
    This is then used to write the .treels file

    > trees : List of trees in a newick format
    """
    out = ""
    for tree in trees.values():
        out += tree.as_string(schema='newick') + "\n"
    return out


def read_newick(filepath):
    """
    Puts the Newick tree read in from a .contree file into a computer friendly tree object

    > filepath : the name of the raw unpruned file containing a newick format tree
    """
    tree = dendropy.Tree.get(path=filepath, schema='newick', preserve_underscores=True)
    print(f">>> Loaded tree: {filepath}")
    # print("printing tree.as_string from read_newick file:",tree.as_string(schema='newick'))
    # print(tree.as_ascii_plot())
    return tree


def find_common_taxa(taxa_sets):
    """
    Take the taxa_set and identify the common taxa across all sets.
    common taxa then written to a list called unique_items_list
    List is used to prune trees to common taxa

    > taxa_sets : list of tree labels extracted from the newick tree
    """
    # An object containing sets of all the taxa found in each of the input newicks
    unique_items = set(taxa_sets[0]) 

    # Iterate over the remaining lists
    for taxa_set in taxa_sets[1:]:
        # .intersection_update updates unique_items with common taxa from all lists
        unique_items.intersection_update(taxa_set) 

    # Convert the set back to a list
    unique_items_list = list(unique_items)
    print(">>> Number of common taxa:",len(unique_items_list))
    return unique_items_list


def get_labels_from_tree(tree):
    """
    Function which returns the labels of the newick tree as an object

    > tree : A newick format tree 
    """
    return[str(item.label) for item in tree.taxon_namespace]


def find_pruned_trees(treefiles):
    """
    Prune the trees inputted, and save them to new files.
    Also prune the aligned dataset to a new file

    > treefiles : 
    """
    tree_labels = {}
    trees = {}
    pruned_trees = {}

    #alignment_file = treefiles[0] # Always put the aligned and trimmed data first
    for f in treefiles:
        if f.endswith(".fas"):
            alignment_file = f
            treefiles.remove(f)

    #treefiles.pop(0) # remove it from the treefiles list

    for filepath in treefiles:
        tree = read_newick(filepath)
        tree_labels[filepath] = get_labels_from_tree(tree)
        trees[filepath] = tree

    # print(tree_labels.values())
    common_taxa = find_common_taxa(list(tree_labels.values()))

    for filepath in treefiles:
        
        # find taxa to remove based on the retrieved common taxa
        taxa_to_be_pruned = [item for item in tree_labels[filepath] if item not in common_taxa]

        # Make pruned version
        pruned_trees[filepath] = dendropy.Tree(trees[filepath])

        # Remove taxa from copy
        pruned_trees[filepath].prune_taxa_with_labels(taxa_to_be_pruned)

        # Save copy to file
        pruned_trees[filepath].write(path= "0PRUNED_" + filepath + ".tre", schema='newick')
        print(">>> Pruned " + filepath)

    # Prune the alignment
    print(">>> Generating alignment file")
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

    print(">>> Pruned concat dataset alignment file")

    return pruned_trees


def check_inputs(file_list):
    """
    Checks that the inputted files are valid

    > file_list : 1 .fas file and all .contree/.tree files being pruned
    """
    count = 0
    fas_files = []
    for f in file_list:
        for i in f:
            print(i)
            if i.endswith(".fas"):
                count += 1
                fas_files.append(i)
    if count > 1:
        raise ValueError(f"We can only have one '.fas' file, but you have submitted {count}: {fas_files}")

##############
### SCRIPT ###
##############

# List all the files being looked at
input_paths = [["gene_data.fas", "z_coi_gene_tree.txt.contree", "concat_tree.txt.contree", "constrained_tree_2.txt.contree"]]
check_inputs(input_paths)

# This is the name of the .treels file which will contain all the pruned trees (extensions will be added automatically)
test_names = ["0_all_trees"]

for i in range(len(input_paths)):
    with open(test_names[i] + "_pruned.treels", 'w') as f:

        print(">>> Finding Pruned Trees for: " + str(", ".join(input_paths[i])), "\n")
        pruned_trees = find_pruned_trees(input_paths[i])
        f.write(concat_trees(pruned_trees))

print(">>>> Script complete <<<")
