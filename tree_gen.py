from dendropy import Tree
from dendropy.simulate import treesim

def extended_birth_death_tree(
        *args,
        extend_relative: float = 0,
        extend_absolute: float = 0,
        **kwargs
) -> Tree:
    """Generate a birth-death tree and extend it by a specified amount.

    This function accepts arbitrary positional and keyword arguments that will
    be passed to treesim.birth_death_tree to generate the initial birth-death
    tree.

    Additionally, this function accepts arguments to specify the amount of time
    by which to extend the initial birth-death tree. This amount may be
    specified in absolute time units or in terms of a proportion of the initial
    tree's leaf distance. If both are specified, the amount of time by which to
    extend the tree is computed by adding the absolute amount to the amount
    computed using the relative proportion.

    Parameters:
        extend_relative (float): The proportion of the tree length to add.
        extend_absolute (float): An absolute amount to add to leaf edges.

    Returns:
        A birth-death tree generated and extended with the specified parameters.
    """
    
    tree = treesim.birth_death_tree(*args, **kwargs)
    leaves = list(tree.leaf_node_iter())
    extend_amount = extend_absolute + \
        extend_relative*leaves[0].distance_from_root()
    for leaf in leaves:
        leaf.edge.length += extend_amount
    return tree
