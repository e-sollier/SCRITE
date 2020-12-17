# Functions used for initializing and optimizing cell lineage trees.
# The tree consists of the possible mutation sites. The cells are attached to the tree.
# A cell attached to the root has no mutation in any of the mutation sites.
# A cell attached to another part of the tree has the mutation it is attached to and all
# mutations of its ancestors.

# Trees are represented as parent vectors: parVec[i] is the index of the parent of the node indexed i.

import numpy as np
import random


def prüfer_to_parent(code):
    """
    Converts Prüfer code to parent vector
    Args:
        code    - prüfer code (list)
        codelen - length of Prüfer code (int)
        
    Returns:
        par_vec - parent vector (list)
    """
    codelen = len(code)
    root = codelen + 2        # same as node count
    par_vec = [0] * (codelen + 1)
    baum = []

    for s in range(codelen):
        comb = code + baum
        for c in range(codelen + 1):
            if not c in comb:
                baum.append(c)
                par_vec[c] = code.pop(0)
                break

    # the last two remaining nodes treated seperately
    last = []
    for l in range(root):
        if not l in baum:
            last.append(l)

    par_vec[last[0]] = last[1]
   
    return par_vec

                
def getRandParentVec(n,seed=None):
    """
    Creates a random parent vector -> This is used to start the tree optimization with a random tree
    Arg:
        n - length of parent vector (int)
        
    Returns:
        Parent vector (list)
    """
    if seed is not None:
        random.seed(seed)
        
    randCode = []
    
    for i in range(n-1):                         # length of Prüfer code
        randCode.append(random.randint(0, n))               # random Prüfer code with n+1 nodes
        
    return prüfer_to_parent(randCode)


def parentVector2ancMatrix(parVec):
    """
    Determines ancestor matrix from parent vector.
    Args:
        parVec    - parent vector (list)
        
    Returns:
        ancMatrix - ancestor matrix (numpy array)
    """
    n = len(parVec)
    ancMatrix = np.zeros((n,n))
    for j in range(n):
        ancMatrix[j][j] = 1     # mutation counted as its own ancestor
     
    for i in range(n):
        a = i
        while a < n:
            if parVec[a] < n:
                ancMatrix[parVec[a]][i] = 1
            a = parVec[a]
        
    return ancMatrix

def parent2children(parVec,reorder = True):
    """
    Converts a parent, where parVec[i] is the parent of i, into a list of children lists, where children[i] is the list of children of i.
    If reorder is True, the list of children is ordered depending on the size of their subtrees.
    """
    nb_mut = len(parVec)

    #Find the children of each node
    children = [[] for v in range(nb_mut+1)]
    for q in range(nb_mut):
        children[parVec[q]].append(q)
    
    #Reorder the children
    if reorder:
        def subtree_size_rec(node,children,subtree_sizes):
            if children[node]==[]:
                subtree_sizes[node] = 0
                return 1
            else:
                for child in children[node]:
                    subtree_sizes[node] += subtree_size_rec(child,children,subtree_sizes)
                return subtree_sizes[node]+1
        subtree_sizes = [0] * (nb_mut+1)
        _ = subtree_size_rec(nb_mut,children,subtree_sizes)
        for q in range(nb_mut+1):
            children[q] = sorted(children[q], key = lambda x :subtree_sizes[x],reverse=False)
    return children

def DFS(parVec):
    """
    Performs a Depth First Traversal of the tree, and returns the (pre)order of the nodes defined by this traversal.
    """
    def DFS_rec(children,root):
        """Recursive DFS of the subtree whose root is 'root' """
        l = [root]
        if children[root]==[]:
            return l
        else:
            for child in children[root]:
                l = l + DFS_rec(children,child)
            return l

    children = parent2children(parVec)
    nb_mut = len(parVec)
    return DFS_rec(children,nb_mut)
    

def proposeNewTree(moveProbsParams, ancMatrix, currTreeParentVec):
    """
    Propose a new mutation tree.
    Args:
        moveProbsParams   - determines the weights of the three move types (prune&re-attach, swap node labels, swap subtrees) (list)
        ancMatrix         - ancestor matrix of current parent vector (numpy array)
        currTreeParentVec - parent vector of current tree (list)
        
    Returns:
        propTreeParentVec - parent vector of proposal tree (list)
    """
    num_mut = len(currTreeParentVec)
    moveType = random.choices([0,1,2], weights = (moveProbsParams[0], moveProbsParams[1], moveProbsParams[2]), k = 1)[0]
    
    if (moveType == 2):  # swap two subtrees in different lineages
        swapNodes = np.random.choice(num_mut, 2, replace=False)
        if (ancMatrix[swapNodes[1]][swapNodes[0]] == 0) and (ancMatrix[swapNodes[0]][swapNodes[1]] == 0):
            propTreeParentVec =  currTreeParentVec[:]
            propTreeParentVec[swapNodes[1]] =  currTreeParentVec[swapNodes[0]]
            propTreeParentVec[swapNodes[0]] =  currTreeParentVec[swapNodes[1]]  
        else:
            moveType = 0
            
    if (moveType == 0):     # prune and re-attach 
        nodeToMove = random.randrange(num_mut)   # pick a node to move with its subtree
        possibleParents = list(np.where(ancMatrix[nodeToMove,:]==0)[0]) + [num_mut]# possible attachment points (cannot attach a node to one of its descendants), root (num_mut) is also possible parent
        possibleParents.remove(currTreeParentVec[nodeToMove])
        if possibleParents==[]:
            moveType=1
        else:           
            newParent = random.choice(possibleParents)   # randomly pick a new parent among available nodes                            
            propTreeParentVec =  currTreeParentVec[:]
            propTreeParentVec[nodeToMove] = newParent
        
        
    if (moveType == 1):         # swap two nodes
        switchNodes = np.random.choice(num_mut, 2, replace=False)
        propTreeParentVec =  currTreeParentVec[:]
        for j in range(num_mut):
            if((currTreeParentVec[j] == switchNodes[0]) and (j != switchNodes[1])):   # change the parent of the children
                propTreeParentVec[j] = switchNodes[1]
                
            if((currTreeParentVec[j] == switchNodes[1]) and (j != switchNodes[0])):
                propTreeParentVec[j] = switchNodes[0]
                        
        propTreeParentVec[switchNodes[1]] = currTreeParentVec[switchNodes[0]]   # switch the nodes
        propTreeParentVec[switchNodes[0]] = currTreeParentVec[switchNodes[1]]
        
        if(propTreeParentVec[switchNodes[1]] == switchNodes[1]):     # if one is the parent of the other
            propTreeParentVec[switchNodes[1]] = switchNodes[0]
        
        if(propTreeParentVec[switchNodes[0]] == switchNodes[0]):
            propTreeParentVec[switchNodes[0]] = switchNodes[1]

    return propTreeParentVec
