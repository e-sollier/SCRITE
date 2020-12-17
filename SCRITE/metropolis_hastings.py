# Functions for running the parameter and tree optimization using a metropolis-hastings algorithm
import math
import random
from random import gauss
import numpy as np
import pandas as pd
from SCRITE.trees import getRandParentVec, proposeNewTree, parentVector2ancMatrix
from SCRITE.scores import calculate_llr_mut, log_scoretree, log_scoreparams, mutation_probabilities_tree



def acceptance(x_logs, x_new_logs):
    """
    If the new score is better, the move is accepted.
    If it is worse, the probability of acceptance depends on how much worse it is.
    Args:
        x_logs      - previous log score (float)
        x_new_logs  - new log score (float)
        
    Returns:
        Accepted or not (bool)
    """
    if x_new_logs > x_logs:
        return True
    
    else:
        accept = np.random.uniform(0,1)
        return (accept < math.e**(x_new_logs - x_logs))
      

def sample_multivariate_normal(params,prior_params):
    """
    Used to draw samples from a multivariate normal distribution
    Args:
        params    - previous parameters (dict)
        prior_params   - prior for the parameters
        
    Returns:
        newParams: proposed parameters
    """
    newParams = {}
    for param in params:
        newParams[param] = params[param] + gauss(0,0.001) * (prior_params[param]["max"] - prior_params[param]["min"])
    return newParams
  
    

def runMCMC(MCMC_params,params,prior_params,df_ref,df_alt,seed=0):
    """
    Runs the Metropolis-Hastings algorithm.
    It samples from the posterior paramter distributions / optimizes the parameters and mutation tree
    Args:
        MCMC params             - parameters for the MCMC: number of loops, number of reps etc...
        loops                   - number of loops within a MCMC (int)
        params                  - initial values for overdispersion_hom, overdispersion_het, dropout, prior_p_mutation (dict)
        prior_params            - min, max, alpha and beta of prior parameter distributions 
        df_ref                  - dataframe of reference read counts
        df_alt                  - dataframe of alternative read counts
        
    Returns:
        sample                  - all samples after burn-in of current log-score, current tree log-score, current params and curent parent vector (list)
        bestTree                - Tree which achieved the highest score
        bestParams              - Paramameters which achieved the highest score
        bestScore               - Highest score obtained among all samples
        df_mut_prob             - Posterior mutation probabilities
    """
    np.random.seed(seed)
    ref=np.array(df_ref)
    alt=np.array(df_alt)
    bool_LOH = [name[:3]=="LOH" for name in df_ref.index]
    indices_LOH = np.where(bool_LOH)[0] 
    
    burnIn = MCMC_params["loops"] * MCMC_params["burnInPhase"]
    eps = 0.00000000001
    bestTree= []
    sample = []
    bestScore = bestTreeLogScore = -10000000000
    bestParams = params
    mutation_probabilities_sample = []
    
    for r in range(MCMC_params["reps"]):       # starts over, but keeps sampling, bestScore, bestTreeLogScore

        currTreeParentVec = getRandParentVec(ref.shape[0])     # start MCMC with random tree
        currTreeAncMatrix =  parentVector2ancMatrix(currTreeParentVec)
        currParams = params
        llr_mut = calculate_llr_mut(params,ref,alt,indices_LOH)
        curr_llr_mut = llr_mut
        currTreeLogScore = log_scoretree(llr_mut, currTreeParentVec, MCMC_params["marginalization"])
        currParamsLogScore = log_scoreparams(currParams, prior_params)
        currScore = currTreeLogScore + currParamsLogScore
        
        if currScore > bestScore:
            bestScore = currScore
            bestTreeLogScore = currTreeLogScore
            bestParams = currParams

        moveAcceptedParams = 0
        totalMovesParams = 0
        moveAcceptedTrees = 0
        totalMovesTrees =  0
            
        for l in range(MCMC_params["loops"]):
            if(l % 10000 == 0):
                print("At mcmc repetition " , r + 1 , "/" , MCMC_params["reps"] , ", step " , l , " best tree score: " , bestTreeLogScore,"\n", sep = "")
                
            rand = np.random.uniform(0,1)
            if rand < MCMC_params["probUpdateParams"]:         # true if this move changes parameters, not the tree
                totalMovesParams += 1
                propParams = sample_multivariate_normal(currParams, prior_params)

                #if the proposed parameters are out of range (or close to), they are not considered
                out_of_range = False
                for param in propParams:
                    out_of_range = out_of_range or propParams[param] < prior_params[param]["min"] + 0.00001
                    out_of_range = out_of_range or propParams[param] > prior_params[param]["max"] - 0.00001
                if out_of_range:
                    continue
                propParamsLogScore = log_scoreparams(propParams, prior_params)
                llr_mut = calculate_llr_mut(propParams, ref,alt,indices_LOH)
                propTreeLogScore = log_scoretree(llr_mut, currTreeParentVec, MCMC_params["marginalization"])
                propScore = propTreeLogScore + propParamsLogScore

                if acceptance(currScore, propScore):  # the proposed move is accepted
                    moveAcceptedParams += 1
                    currTreeLogScore  = propTreeLogScore
                    currParams = propParams.copy()
                    currParamsLogScore = propParamsLogScore
                    currScore = propScore
                    curr_llr_mut = llr_mut
        
            else: # if the move changes the tree not the parameter
                totalMovesTrees += 1
                propTreeParentVec = proposeNewTree(MCMC_params["moveProbs"], currTreeAncMatrix, currTreeParentVec)
                propTreeLogScore = log_scoretree(curr_llr_mut, propTreeParentVec, MCMC_params["marginalization"])
                if acceptance(currTreeLogScore, propTreeLogScore):
                     # the proposed tree is accepted
                    moveAcceptedTrees += 1
                    currTreeAncMatrix = parentVector2ancMatrix(propTreeParentVec)
                    currTreeParentVec = propTreeParentVec                                
                    currTreeLogScore  = propTreeLogScore                     
                    currScore = currTreeLogScore + currParamsLogScore
        
                    
                
            if(l >= burnIn and l % MCMC_params["sampleStep"] == 0):
                sample.append([currScore, currTreeLogScore, currParams, currTreeParentVec])
                mutation_probabilities_sample.append(mutation_probabilities_tree(curr_llr_mut,currTreeParentVec))
                
            if(currScore > bestScore + eps):
                bestTree = currTreeParentVec        
                bestTreeLogScore = currTreeLogScore
                bestScore = currScore
                bestParams = currParams.copy()

    noStepsAfterBurnin = MCMC_params["reps"] * (MCMC_params["loops"] - burnIn)

    print( "best log score for tree: " , bestTreeLogScore)
    print( "total #steps after burn-in: ", noStepsAfterBurnin)
    print( "percentage of new Parameters accepted:", (moveAcceptedParams / (totalMovesParams+1)) * 100, "%")
    print( "percentage of Tree moves accepted:", (moveAcceptedTrees / totalMovesTrees) * 100, "%")
    
    if(MCMC_params["probUpdateParams"] != 0):
        print( "best value for overdispersion_hom: " , bestParams["overdispersion_hom"])
        print( "best value for overdispersion_het: " , bestParams["overdispersion_het"])
        print( "best value for  dropout: " , bestParams["dropout"])
        print( "best value for  sequencing_error_rate: " , bestParams["sequencing_error_rate"])
        print( "best log score for (Tree, Params): " , bestScore)
    df_mut_prob = pd.DataFrame(np.mean(mutation_probabilities_sample,axis=0),index = df_ref.index,columns = df_ref.columns)
    return sample, bestTree, bestParams, bestScore, df_mut_prob
