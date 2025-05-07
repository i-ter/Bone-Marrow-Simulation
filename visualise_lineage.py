import networkx as nx
import matplotlib.pyplot as plt
import re

LINEAGE_TREE = {
       'HSC': ['MPP1'],
       'MPP1': ['MPP2'],
       'MPP2': ['MPP3'],
       'MPP3': ['MPP4'],
       'MPP4': ['MPP5'],
       'MPP5': ['CMP', 'CLP'],
       'CMP': ['MEP', 'GMP'],
       'GMP': ['Myeloblast1'],
       'MEP': ['Megakaryocyte', 'Erythroblast1'],
       'Megakaryocyte': ['Platelet'],
       'Erythroblast1': ['Erythroblast2'],
       'Erythroblast2': ['Erythroblast3'],
       'Erythroblast3': ['Erythroblast4'],
       'Erythroblast4': ['Erythroblast5'],
       'Erythroblast5': ['Erythroblast6'],
       'Erythroblast6': ['Erythroblast7'],
       'Erythroblast7': ['Erythroblast8'],
       'Erythroblast8': ['RBC'],
       'Myeloblast1': ['Myeloblast2'],
       'Myeloblast2': ['Myeloblast3'],
       'Myeloblast3': ['Myeloblast4'],
       'Myeloblast4': ['Myeloblast5'],
       'Myeloblast5': ['Myeloblast6'],
       'Myeloblast6': ['Myeloblast7'],
       'Myeloblast7': ['Myeloid'],
       'CLP': ['Lymphocyte1'],
       'Lymphocyte1': ['Lymphocyte2'],
       'Lymphocyte2': ['Lymphocyte3'],
       'Lymphocyte3': ['Lymphocyte4'],
       'Lymphocyte4': ['Lymphocyte5'],
       'Lymphocyte5': ['Lymphocyte6'],
       'Lymphocyte6': ['Lymphocyte7'],
       'Lymphocyte7': ['Bcell'],
   }

def visualize_lineage_terminal(tree):
    """
    Visualize the lineage tree in the terminal using ASCII characters.
    """
    def print_tree(node, prefix='', is_last=True):
        # Print the current node
        connector = '└── ' if is_last else '├── '
        print(f"{prefix}{connector}{node}")
        
        # Prepare prefix for children
        new_prefix = prefix + ('    ' if is_last else '│   ')
        
        # Get children of the current node
        children = tree.get(node, [])
        for i, child in enumerate(children):
            is_last_child = i == len(children) - 1
            print_tree(child, new_prefix, is_last_child)

    # Start with the root node (HSC)
    if 'HSC' in tree:
        print("--------------------------------")
        print("Cell Lineage Tree")
        print_tree('HSC')
        print("--------------------------------")
    else:
        print("Root node 'HSC' not found in tree.")

if __name__ == "__main__":

    lineage_tree = LINEAGE_TREE
    visualize_lineage_terminal(lineage_tree)
