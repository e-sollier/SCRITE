# Functions to convert read counts to mutation probabilities and to calculate the paramter and tree log-scores

import math
import numpy as np
from scipy.special import loggamma
from SCRITE.trees import DFS, parentVector2ancMatrix, parent2children


# Convert read counts to mutation probabilities

def calculate_llr_mut(params,ref,alt,indices_LOH=[]):
    """
    Args:
        params                  - dictionary of parameters (overdispersion_wt, overdispersion_mut, dropout, sequencing_error_rate)
        alt                     - alternative read counts (list)
        ref                     - wildtype/reference read counts (list)
        
    Returns:
        llr_mut                 - log likelihood ratio of the data given that the cell has a mutation compared to not having the mutation
    """
    alpha_hom = params["overdispersion_hom"] * params["sequencing_error_rate"]
    beta_hom = params["overdispersion_hom"] * (1 - params["sequencing_error_rate"])

    alpha_het = params["overdispersion_het"] * (1/2 - 1/3 * params["sequencing_error_rate"]) 
    beta_het = params["overdispersion_het"] * (1 - (1/2 - 1/3 * params["sequencing_error_rate"]))
    

    gamma_const_hom = loggamma(params["overdispersion_hom"]) - loggamma(alpha_hom) - loggamma(beta_hom)
    gamma_const_het = loggamma(params["overdispersion_het"]) - loggamma(alpha_het) - loggamma(beta_het)

    S = alt
    C = alt+ref

    #Vectorize computations
    log_gamma_s_alpha_hom = loggamma(S + alpha_hom)
    log_gamma_c_s_beta_hom = loggamma(ref + beta_hom)
    log_gamma_c_overdisp_hom = loggamma(C + params["overdispersion_hom"])
    logp_hom = log_gamma_s_alpha_hom + log_gamma_c_s_beta_hom - log_gamma_c_overdisp_hom + gamma_const_hom

    log_gamma_s_beta_hom = loggamma(S + beta_hom)
    log_gamma_c_s_alpha_hom = loggamma(ref + alpha_hom)
    p_mut_homozygous = np.exp(log_gamma_s_beta_hom + log_gamma_c_s_alpha_hom + gamma_const_hom - log_gamma_c_overdisp_hom) # dropout of the ref allele

    log_gamma_s_alpha_het = loggamma(S+alpha_het)
    log_gamma_c_s_beta_het = loggamma(ref + beta_het)
    log_gamma_c_overdisp_het = loggamma(C + params["overdispersion_het"])
    logp_het_without_dropout = log_gamma_s_alpha_het + log_gamma_c_s_beta_het + gamma_const_het - log_gamma_c_overdisp_het
    logp_het_without_dropout[np.abs(logp_het_without_dropout)<10e-12]=0 # faster
    p_het_without_dropout = np.exp(logp_het_without_dropout)
    
    p_het = params["dropout"]/2 * np.exp(logp_hom) + params["dropout"]/2 * p_mut_homozygous + (1-params["dropout"]) * p_het_without_dropout

    #avoid log(0)
    nan_indices = (p_het<=10e-300)
    p_het[nan_indices] = np.nan
    logp_het = np.log(p_het)
    logp_het[nan_indices] = logp_het_without_dropout[nan_indices]

    llr = logp_het - logp_hom
    llr[indices_LOH,:]*=-1 # For alleles which are lost, we go from heterozygous to homozygous, so the log likelihood ratio is reversed.
    return llr


def log_pdf(a, b, x):
    """
    Logarithm of probability density function (pdf) of the beta distribution at point x
    Args:
        a - alpha (float)
        b - beta (float)
        x - range: (0,1) / probability density is determined at point x (float)
        
    Returns:
        Logarithm of probability density function (float)
    """
    return math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + (a - 1) * math.log(x) + (b - 1) * math.log(1 - x)


def log_likelihood_attachments(parVec,llr_mut):
    num_mut,num_cells = np.shape(llr_mut)
    
    DFS_order = DFS(parVec)                         
    log_likelihoods = np.zeros((num_mut+1,num_cells))

    for k in range(1, num_mut + 1):        
        node = DFS_order[k]
        # As we go down in the tree, a mutation replaces a reference, so the log likelihood ratio is added
        log_likelihoods[node,:] = log_likelihoods[parVec[node],:] + llr_mut[node,:]

    return log_likelihoods


def log_scoretree(llr_mut, parVec, marginalization):
    """
    Log score of the tree
    Args:
        llr_mut         - log likelihood ratio of the data given that the cell has a mutation compared to not having the mutation (matrix)
        parVec          - parent vector of tree (list)
        marginalization - if true the attachment points are marginalized (bool)
        
    Returns:
        log_score       - logarithmic tree score (float)
    """     
    log_likelihoods = log_likelihood_attachments(parVec,llr_mut)


    if marginalization == False:
        log_score = np.sum(np.max(log_likelihoods,axis=0))
        
    if marginalization == True:
        log_score = np.sum(np.max(log_likelihoods,axis=0)) + np.sum(np.log(np.mean(np.exp(log_likelihoods - np.max(log_likelihoods,axis = 0)),axis=0)))

    return log_score



def mutation_probabilities_tree(llr_mut,parVec):
    """Compute the mutation probabilities given a tree"""
    num_mut,num_cells = np.shape(llr_mut)
    ancMatrix = parentVector2ancMatrix(parVec)

    scores = log_likelihood_attachments(parVec,llr_mut)
    likelihood_attachments = np.exp(scores - np.max(scores,axis=0)) #offset the probabilities for numerical stability. This offset will cancel out.
    denominators = np.sum(likelihood_attachments,axis=0)
    mutation_probabilities = np.zeros((num_mut,num_cells))
    for j in range(num_mut):
        attachments_with_mutation = np.concatenate([ancMatrix[j,:],[0]]) # the descendants of the mutation have the mutation
        mutation_probabilities[j,:] = np.dot(attachments_with_mutation.T,likelihood_attachments) / denominators # sum the probabilities of the attachments where mutation j occurs
    mutation_probabilities= np.nan_to_num(mutation_probabilities,nan=0.5) # nan when denominator is 0, which occurs when all attachments have very low probability
    return mutation_probabilities


def log_scoreparams(params, prior_params):
    """
    Args:
        params               - overdispersion_wt, overdispersion_mut, dropout, prior_p_mutation (dict)
        prior_params         - priors for the parameters
        
    Returns:
        Logarithmic parameters score (float)
    """
    log_score = 0
    for param in params:
        log_score += log_pdf(prior_params[param]["alpha"],prior_params[param]["beta"],params[param]/prior_params[param]["max"])      
    
    return log_score
