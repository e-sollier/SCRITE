# Functions used to create the desired output files

import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
import seaborn as sns
from SCRITE.scores import calculate_llr_mut, log_likelihood_attachments
from SCRITE.trees import parentVector2ancMatrix, DFS, parent2children



def get_likelihood_attachments(parVec,params,ref,alt,indices_LOH):
    """
    Given a mutation tree, compute the likelihoods of the attachments of the cells to the mutations.
    Args:
        parVec: mutation tree given as a list of parents
        params                  - overdispersion_wt, overdispersion_mut, dropout, sequencing_error_rate (dict)
        ref                     - reference read counts (list)
        alt                     - alternative read counts (list)
    """
    llr_mut = calculate_llr_mut(params, ref,alt,indices_LOH=indices_LOH)

    # Compute the attachment probabilities by doing a depth first traversal of the mutation tree
    scores = log_likelihood_attachments(parVec,llr_mut)
    scores = np.exp(scores - np.max(scores,axis=0)) #offset the probabilities for numerical stability. This offset will cancel out.
    return scores


def get_optimal_attachments(parVec, params,ref,alt,indices_LOH):
    """
    Determines the optimal attachment points of the cells to the mutation tree
    Args:
        parVec                  - parent vector (list)
        params                  - overdispersion_wt, overdispersion_mut, dropout, sequencing_error_rate (dict)
        ref                     - reference read counts (list)
        alt                     - alternative read counts (list)
        
        
    Returns:
        attachmentPoints        - optimal attachment points (list)
    """    
    likelihood_attachments = get_likelihood_attachments(parVec,params,ref,alt,indices_LOH)
    optimal_attachments = np.argmax(likelihood_attachments,axis=0)

    return optimal_attachments
  
    

def graphviz(parVec,params,gene_names,ref,alt,indices_LOH,includeCells=False):
    """
    Create a graphviz file of the mutation tree
    Args:
        parVec                  - parent vector of tree (list)
        params                  - overdispersion_wt, overdispersion_mut, dropout, sequencing_error_rate (dict)
        gene_names              - The names of the mutation sites (list)
        ref                     - Counts for the ref allele
        alt                     - Counts for the alt allele
        indices_LOH             - Indices corresponding to LOH events
        includeCells            - if True, include the cells in the mutation tree
        
    Returns:
        gv                      - graphviz file (string)
    """    
    num_mut,num_cells = np.shape(alt)
    gene_names.append("Root")
    
    gv = "digraph G {\n"
    gv += "node [color=deeppink4, style=filled, fontcolor=white];\n"

    for i in range(num_mut):
        # the parantheses around the gene names help displaying them if they contain special characters
        gv += ("\"" + gene_names[parVec[i]] + "\"" + " -> "  + "\"" + gene_names[i]  + "\"" + ";\n" ) 

    if includeCells:
        gv += "node [color=lightgrey, style=filled, fontcolor=black];\n"
                                
        attachmentPoints = get_optimal_attachments(parVec, params, ref,alt,indices_LOH)
        
        for y, a in enumerate(attachmentPoints):
            gv += "\"" + gene_names[a] + "\"" + " -> s"  + str(y + 1) + ";\n"

    gv += "}\n"
    return gv



def getPosteriorMutationProbabilities(parVecSample, paramsSample,ref,alt,indices_LOH):
    """
    Args:
        parVecSample            - Sample of parent vectors (list)
        paramsSample            - Sample of params: overdispersion_wt, overdispersion_mut, dropout, prior_p_mutation (dict)
        ref                     - reference read counts (array)
        alt                     - alternative read counts (array)
        
    Returns:
        attachmentPoints        - optimal attachment points (list)
    """    
    num_mut,num_cells = np.shape(ref)
    mutation_probabilities = np.zeros((num_mut,num_cells))
    # For each tree in the sample
    for sample_index in range(len(parVecSample)):
        params = paramsSample[sample_index]
        parVec = parVecSample[sample_index]
        ancMatrix = parentVector2ancMatrix(parVec)

        likelihood_attachments = get_likelihood_attachments(parVec,params,ref,alt,indices_LOH)
        denominators = np.sum(likelihood_attachments,axis=0)
        prob_mutations = np.zeros((num_mut,num_cells))
        for j in range(num_mut):
            attachments_with_mutation = np.concatenate([ancMatrix[j,:],[0]]) # the descendants of the mutation have the mutation
            prob_mutations[j,:] = np.dot(attachments_with_mutation.T,likelihood_attachments) / denominators # sum the probabilities of the attachments where mutation j occurs
        prob_mutations = np.nan_to_num(prob_mutations,nan=0.5) # nan when denominator is 0, which occurs when all attachments have very low probability
        mutation_probabilities = mutation_probabilities + prob_mutations

    mutation_probabilities = mutation_probabilities / len(parVecSample)
    return mutation_probabilities

def reorder_mutations_cells(parVec,params,ref,alt,indices_LOH):
    """ 
    Reorder the mutations and the cells to reflect clustering induced by the mutation tree.
    The order of the mutations is defined by a Depth First traversal of the mutation tree,
    and the order of the cells is based on their best attachment in the mutation tree.
    """
    DFS_order = DFS(parVec)
    attachments = get_optimal_attachments(parVec,params,ref,alt,indices_LOH)
    order_cells = []
    for mutation in DFS_order:
        #Find all cells attached to this mutation
        cells_attached = list(np.where(attachments == mutation)[0])
        order_cells = order_cells + cells_attached
    return DFS_order[1:],order_cells


def create_dendrogram_from_tree(parVec):
    """Converts a tree (given as a parent vector) to a dendrogram. It is used to visualize the tree next to the heatmap."""

    def recursive_dendrogram(node,children_list,cluster_count,depth):
        children = children_list[node]
        clusters = []
        new_clusters = []
        size=1
        for child in children:
            cluster_id,cluster_size,clusters_subtree,cluster_count = recursive_dendrogram(child, children_list,cluster_count,depth+1)
            clusters = clusters + clusters_subtree
            new_clusters.append((cluster_id,cluster_size))
            size+=cluster_size
        new_clusters.append((node,1))
        if len(new_clusters)>1:
            for i in range(len(new_clusters)-1):
                clusters.append([new_clusters[i][0],new_clusters[i+1][0],-depth,new_clusters[i][1] + new_clusters[i+1][1]])
                cluster_count+=1
                new_clusters[i+1] =  (cluster_count, new_clusters[i][1] + new_clusters[i+1][1])
            clusters[-1][2] = 0.5 - depth
            return (cluster_count,size,clusters,cluster_count)
        else:
            return (node,size,clusters,cluster_count)

    n = len(parVec)
    children_list = parent2children(parVec)

    cluster_id,size,clusters,cluster_count = recursive_dendrogram(n,children_list,n-1,1)
    Z = np.array(clusters[:-1])
    Z[:,2] = Z[:,2] - np.min(Z[:,2])
    return Z

def get_cells_annotations(metadata,cell_names):
    if metadata is None:
        return None
    l = []
    for x in cell_names:
        if x in metadata.index:
            if metadata.loc[x,"neoplastic"]=="Neoplastic":
                l.append("red")
            else:
                l.append("blue")
        else:
            l.append("grey")
    cells_annotations = pd.DataFrame({"Neoplasticity":np.array(l)},index = cell_names)
    return cells_annotations



def plot_mutation_probabilities(mutation_probabilities,tree,params,df_ref,df_alt,cells_annotations=None,output_file=None):
    """
    Plot the mutation probabilities for each cell and each mutation, and reorder the cells and mutations based on the tree.
    Also plot a dendrogram next to the heatmap to show the tree.
    """
    dendrogram = create_dendrogram_from_tree(tree)
    bool_LOH = [name[:3]=="LOH" for name in df_ref.index]
    indices_LOH = np.where(bool_LOH)[0]
    order_mutations, order_cells = reorder_mutations_cells(tree,params,np.array(df_ref),np.array(df_alt),indices_LOH)

    cell_names = mutation_probabilities.columns
    cell_names = np.array(cell_names)[order_cells]
    mut_names = mutation_probabilities.index
    reordered_p_mut = mutation_probabilities.loc[:,cell_names]


    if cells_annotations is not None:
        cells_annotations = cells_annotations.loc[cell_names,:]

    fontsize = 3 if len(mut_names) >=120 else 6
    yticks = 1 if len(mut_names)<=240 else 0
    if cells_annotations is None:
        cg = sns.clustermap(reordered_p_mut,row_linkage = dendrogram, col_cluster=False, cmap='RdYlBu_r',vmin=0,vmax=1,cbar_kws= {"label":"Mutation Probabilities"},
        yticklabels=yticks,xticklabels=0, annot_kws={"size": 4})
    else:
        cg = sns.clustermap(reordered_p_mut,row_linkage = dendrogram, col_cluster=False, cmap='RdYlBu_r',col_colors=cells_annotations,
        cbar_kws= {"label":"Mutation Probabilities"},
        yticklabels=yticks,xticklabels=0, annot_kws={"size": 4},vmin=0,vmax=1)
        
    cg.ax_heatmap.invert_yaxis()
    cg.ax_row_dendrogram.invert_yaxis()
    cg.ax_heatmap.set_yticklabels(cg.ax_heatmap.get_ymajorticklabels(), fontsize = fontsize)
    cg.ax_heatmap.set_xlabel("Cells")
    cg.ax_heatmap.set_ylabel("Mutations")
    if output_file is None:
        plt.show()
    else:
        plt.savefig(output_file,dpi=200)