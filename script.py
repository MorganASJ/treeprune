import dendropy
import os
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
    print("Number of common taxa:",len(unique_items_list), "\n")
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
        pruned_trees[filepath].write(path= "0PRUNED_" + filepath + ".pruned.tre", schema='newick')
        print("Pruned tree: " + filepath)

    # Prune the alignment
    print("\nGenerating alignment file")
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

    with open("0PRUNED_" + alignment_file + ".pruned.fas", "w") as f:
        f.write(pruned_alignment)

    print("Pruned alignment:", alignment_file)

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
            if i.endswith(".fas"):
                count += 1
                fas_files.append(i)
    if count > 1:
        raise ValueError(f"We can only have one '.fas' file, but you have submitted {count}: {fas_files}")
    else:
        print("Files have been accepted. Script will now continue. \n")
    

def get_files(tree_extension):
    """
    Retrieves all files which can be used in analysis from the directory
    this script was run from

    > tree_extension : string containing the file type for the newick trees (".contree" for consensus", .treefile" for ML )
    
    """
    inputs = []
    input_path = []
    list_of_files = os.listdir()

    for f in list_of_files:
        if f.endswith(".fas") or f.endswith(tree_extension): # may have to change the contree ending
            inputs.append(f)
    inputs = sorted(inputs) # sort alphatbeticall to always have the same order
    
    inputs = [item for item in inputs if "PRUNED" not in item]
            
    input_path.append(inputs)
    print("The following files have been found:\n", ">>>", input_path)
    return input_path


def write_log(name_of_log, tree_files, alignment_files):

    with open(name_of_log, 'w') as file:
        file.write("Order of treefiles pruned, order for AU test:")
        for f in tree_files:
            file.write(f)
            file.write("\n")
        
        file.write("Alginment pruned:")
        file.write(alignment_files[0])
        file.write("\n")
    print("Written to: ", name_of_log)



##############
### SCRIPT ###
##############

# List all the files being looked at
# input_paths = [["gene_data.fas", "z_coi_gene_tree.txt.contree", "2_concat_tree.txt.contree", "1_constrained_tree_2.txt.contree"]]
# check_inputs(input_paths)
# print(input_paths)

# Automatically locate the files in path
input_paths = get_files(".contree")
print("input paths: -------- ", input_paths)
alignment_paths = []
tree_paths = []
for f in input_paths:
        for i in f:
            if i.endswith(".fas"):
                alignment_paths.append(i)
            if i.endswith(".contree"): #CHANGE TO SUIT
                tree_paths.append(i)

check_inputs(input_paths)


# This is the name of the .treels file which will contain all the pruned trees (extensions will be added automatically)
test_names = ["all_trees"]

for i in range(len(input_paths)):
    with open(test_names[i] + ".pruned.treels", 'w') as f:
        pruned_trees = find_pruned_trees(input_paths[i])
        f.write(concat_trees(pruned_trees))


write_log("LOGNAME.log",tree_paths, alignment_paths)

print("\n>>>> Script complete <<<\n")
