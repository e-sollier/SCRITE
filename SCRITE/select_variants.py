import numpy as np
import pandas as pd

def filter_min_distance(loci,df_ref,df_alt,min_dist=10000):
    """Only keep loci which are at least min_dist bases apart. This is used to make sure that there is at most one locus per gene,
    so that the allelic dropouts are independant.
    When one locus has to be removed, remove the one with the lowest coverage."""
    loci_kept = []
    last_locus_added = ("0",-10000)
    last_coverage=0
    for locus in loci:
        _,chromosome, pos = locus.split("_")
        coverage = np.sum((df_ref.loc[locus,:] + df_alt.loc[locus,:])>0)
        pos = int(pos)
        if chromosome != last_locus_added[0] or abs(pos-last_locus_added[1])>min_dist:
            last_locus_added = (chromosome,pos)
            loci_kept.append(locus)
            last_coverage=coverage
        elif last_coverage*1.7<coverage: # if this new position is close to the previous one but has a much higher coverage, use this one instead
            _ = loci_kept.pop()
            last_locus_added = (chromosome,pos)
            loci_kept.append(locus)
            last_coverage=coverage

    return loci_kept


def select_loci_LOH(df_ref,df_alt,SNPs=None,minPercentageCellsAllele=10):
    """
    Select loci for which one allele might be lost in some cells.
    Args:
        df_ref: dataframe of the counts for the ref allele (mutations in rows and cells in columns)
        df_alt: dataframe of the counts for the alt allele
        SNPs: set of positions corresponding to common SNPs (optional). If it is given, only these positions will be considered for LOH
        minPercentageCellsAllele: minimum percentage of the cells with at least one read at this position, for each allele
    Returns;
        df_ref_selected: dataframe of the counts for the ref allele, for the selected loci
        df_alt_selected: dataframe of the counts for the alt allele (that might be lost), for the selected loci
    """
    ref_selected = []
    alt_selected = []
    loci_selected = []
    nb_cells = df_ref.shape[1]

    for locus in df_ref.index:
        if SNPs is None or locus in SNPs:
            chr,_ = locus.split("_")
            ref = df_ref.loc[locus,:]
            alt = df_alt.loc[locus,:]
            cells_alt = np.sum(alt>0)
            cells_ref = np.sum(ref>0)
            if cells_ref / float(nb_cells) >= minPercentageCellsAllele/100 and cells_alt / float(nb_cells) >= minPercentageCellsAllele/100 and (not chr in ["chrM","chrX"]):
                #only consider loci for which many cells have reads
            
                cells_alt_noref = np.sum((alt>4) & (ref==0)) # number of cells where the alt allele is expressed but not the ref allele
                cells_ref_noalt = np.sum((ref>4) & (alt==0)) # number of cells where the ref allele is expressed but not the alt allele

                if cells_alt_noref >= 1.6 * cells_ref_noalt + 4 and cells_alt/nb_cells >= 2*minPercentageCellsAllele/100:
                    # the ref allele might be lost in some cells
                    ref_expressed = np.where(ref>0)[0]
                    ratio_alt_ref = np.mean(alt[ref_expressed] / ref[ref_expressed])
                    # make sure that the ref allele is not consistently expressed less than the alt allele
                    if ratio_alt_ref <= 2: 
                        ref_selected.append(alt)
                        alt_selected.append(ref) #the allele that is lost is considered as the alt allele
                        loci_selected.append("LOH_"+locus)
                elif cells_ref_noalt >= 1.6 * cells_alt_noref + 4 and cells_ref/nb_cells >= 2*minPercentageCellsAllele/100: 
                    # the alt allele might be lost in some cells
                    alt_expressed = np.where(alt>0)[0]
                    ratio_ref_alt = np.mean(ref[alt_expressed] / alt[alt_expressed])
                    # make sure that the alt allele is not consistently expressed less than the ref allele
                    if ratio_ref_alt <= 2:
                        ref_selected.append(ref)
                        alt_selected.append(alt) #the allele that is lost is considered as the alt allele
                        loci_selected.append("LOH_"+locus)

    df_ref_selected = pd.DataFrame(ref_selected, index = loci_selected, columns = df_ref.columns)
    df_alt_selected = pd.DataFrame(alt_selected, index = loci_selected, columns = df_alt.columns)

    # Filter loci which are too close from each other, so that allelic dropouts remain independant
    loci_selected_dist_filtered = filter_min_distance(loci_selected,df_ref_selected,df_alt_selected,min_dist=10000)
    df_ref_selected = df_ref_selected.loc[loci_selected_dist_filtered,:]
    df_alt_selected = df_alt_selected.loc[loci_selected_dist_filtered,:]
    if df_alt_selected.shape[0]>400:
        #if too many sites were selected, use stricter filtering criterions
        pos_selected = [x[4:] for x in loci_selected]
        return select_loci_LOH(df_ref.loc[pos_selected,:],df_alt.loc[pos_selected,:],SNPs=SNPs,minPercentageCellsAllele = minPercentageCellsAllele+2)
    return df_ref_selected,df_alt_selected

    
def select_loci_classified(df_ref,df_alt,neoplastic_cells,regular_cells,SNPs=None,minPercentageCellsExpressed=8):
    """
    Select loci where there might be a LOH, when we already know which cells are 
    neoplastic and which cells are healthy.
    Args:
        df_ref: dataframe of the counts for the ref allele (mutations in rows and cells in columns)
        df_alt: datafrale if the counts for the alt allele
        neoplastic_cells: cells identified as neoplastic, which might contain LOH events or somatic mutations
        regular_cells: cells identified as healthy (ie immune cells...) which should not contain LOH events nor somatic mutations
        SNPs: list of positions corresponding to common SNPs (optional). If it is given, only these positions will be considered for LOH, 
            and only the other positions will be considered for somatic mutations
        minPercentageCellsExpressed: minimum percentage of the cells with at least one read at this position
    Returns:
        df_ref_selected: dataframe of the counts for the ref allele, for the selected loci
        df_alt_selected: dataframe of the counts for the alt allele (that might be lost), for the selected loci
    """
    ref_selected = []
    alt_selected = []
    loci_selected = []

    nb_cells = df_ref.shape[1]
    nb_neoplastic = len(neoplastic_cells)
    nb_regular = len(regular_cells)

    for locus in df_ref.index:
        chr,_ = locus.split("_")
        ref = df_ref.loc[locus,:]
        alt = df_alt.loc[locus,:]
        cells_expressed = np.sum((alt+ref)>0)
        if cells_expressed / float(nb_cells) >= minPercentageCellsExpressed/100 and (not chr in ["chrM","chrX"]):
            #only consider loci for which many cells have reads
            ref_neoplastic = np.sum(df_ref.loc[locus,neoplastic_cells]>0)
            ref_regular = np.sum(df_ref.loc[locus,regular_cells]>0)
            alt_neoplastic = np.sum(df_alt.loc[locus,neoplastic_cells]>0)
            alt_regular = np.sum(df_alt.loc[locus,regular_cells]>0)

            if (SNPs is None or not locus in SNPs) and \
                ref_regular / nb_regular >= 0.7* minPercentageCellsExpressed / 100 and alt_regular < 0.03 * ref_regular and \
                alt_neoplastic>max(5,0.10 * ref_neoplastic):
                #Somatic mutation
                #Ref allele is expressed in many regular cells, alt allele is not expressed in regular cells
                #and alt allele is expressed at least 10% as much as the ref allele (to account for the fact that not all tumor cells may have the mutation)
                ref_selected.append(ref)
                alt_selected.append(alt)
                loci_selected.append("MUT_"+locus)

            elif SNPs is None or locus in SNPs:
                #LOH
                alt_noref_regular = np.sum((df_alt.loc[locus,regular_cells]>4) & (df_ref.loc[locus,regular_cells]==0))
                alt_noref_neoplastic = np.sum((df_alt.loc[locus,neoplastic_cells]>4) & (df_ref.loc[locus,neoplastic_cells]==0))
                ref_noalt_regular = np.sum((df_ref.loc[locus,regular_cells]>4) & (df_alt.loc[locus,regular_cells]==0))
                ref_noalt_neoplastic = np.sum((df_ref.loc[locus,neoplastic_cells]>4) & (df_alt.loc[locus,neoplastic_cells]==0))
                if ref_noalt_neoplastic > 1.5 * alt_noref_neoplastic and ref_noalt_regular <= 2.5 * alt_noref_regular:
                    # Alt allele lost in neoplastic cells
                    if alt_regular / (nb_regular+1) >= 0.5*minPercentageCellsExpressed /100 and ref_neoplastic / (nb_neoplastic+1) >= 0.5*minPercentageCellsExpressed /100:
                        alt_expressed = np.where(alt>0)[0]
                        ratio_ref_alt = np.mean(ref[alt_expressed] / alt[alt_expressed])
                        # make sure that the alt allele is not consistently expressed less than the ref allele
                        if ratio_ref_alt <= 2: 
                            ref_selected.append(ref)
                            alt_selected.append(alt)
                            loci_selected.append("LOH_"+locus)
                elif alt_noref_neoplastic > 1.5 * ref_noalt_neoplastic and alt_noref_regular <= 2.5 * ref_noalt_regular:
                    # Ref allele lost in neoplastic cells
                    if ref_regular / (nb_regular+1) >= 0.5*minPercentageCellsExpressed /100 and alt_neoplastic / (nb_neoplastic+1) >= 0.5*minPercentageCellsExpressed /100:
                        ref_expressed = np.where(ref>0)[0]
                        ratio_alt_ref = np.mean(alt[ref_expressed] / ref[ref_expressed])
                        # make sure that the ref allele is not consistently expressed less than the alt allele
                        if ratio_alt_ref <= 2: 
                            ref_selected.append(alt)
                            alt_selected.append(ref)
                            loci_selected.append("LOH_"+locus)

    df_ref_selected = pd.DataFrame(ref_selected, index = loci_selected, columns = df_ref.columns)
    df_alt_selected = pd.DataFrame(alt_selected, index = loci_selected, columns = df_alt.columns)
    # Filter loci which are too close from each other, so that allelic dropouts remain independant
    loci_selected_dist_filtered = filter_min_distance(loci_selected,df_ref_selected,df_alt_selected,min_dist=10000)
    df_ref_selected = df_ref_selected.loc[loci_selected_dist_filtered,:]
    df_alt_selected = df_alt_selected.loc[loci_selected_dist_filtered,:]

    if df_alt_selected.shape[0]>600:
        #if too many sites were selected, use stricter filtering criterions
        pos_selected = [x[4:] for x in loci_selected]
        return select_loci_classified(df_ref.loc[pos_selected,:],df_alt.loc[pos_selected,:],neoplastic_cells,regular_cells,SNPs,minPercentageCellsExpressed = minPercentageCellsExpressed+2)
    return df_ref_selected,df_alt_selected

def identify_neoplastic_regular(df_mut_prob):
    """
    Given a dataframe of mutation probabilities, identify some cells as healthy and some as neoplastic.
    The cells with more mutations are considered neoplastic.
    Not all cells are classified as neoplastic or regular: those for which it's unclear are not included
    Args:
        df_mut_prob: dataframe of posterior mutation probabilities
    Returns:
        neoplastic_cells: array of cell IDs for the cells identified as neoplastic
        regular_cells: array of cell IDs for the cells identified as regular
    """
    number_events_per_cell = np.sum(df_mut_prob>0.7,axis=0)
    sorted_nb_events = sorted(list(number_events_per_cell))
    sum_probs_per_cell = np.sum(df_mut_prob,axis=0)
    sorted_sum_probs = sorted(list(sum_probs_per_cell))

    width = len(sorted_nb_events) //5
    low_nb_events = np.mean(sorted_nb_events[:width])
    high_nb_events = np.mean(sorted_nb_events[-width:])
    regular_threshold_nb_events = 0.7 * low_nb_events + 0.3 * (low_nb_events+high_nb_events)/2
    neoplastic_threshold_nb_events = 0.1 * high_nb_events + 0.9*(low_nb_events+high_nb_events)/2

    low_sum_probs = np.mean(sorted_sum_probs[:width])
    high_sum_probs = np.mean(sorted_sum_probs[-width:])
    regular_threshold_sum_probs = 0.8 * low_sum_probs + 0.2 * (low_sum_probs+high_sum_probs)/2
    neoplastic_threshold_sum_probs = 0.1 * high_sum_probs + 0.9*(low_sum_probs+high_sum_probs)/2
    
    regular_cells = np.where((number_events_per_cell <= regular_threshold_nb_events) & (sum_probs_per_cell<=regular_threshold_sum_probs))[0]
    neoplastic_cells = np.where((number_events_per_cell >= neoplastic_threshold_nb_events) & (sum_probs_per_cell>=neoplastic_threshold_sum_probs))[0]
    regular_cells = np.array(df_mut_prob.columns[regular_cells])
    neoplastic_cells = np.array(df_mut_prob.columns[neoplastic_cells])

    return neoplastic_cells,regular_cells




def select_events_present(df_mut_prob,tree):
    """
    Select the events (mutations/LOH) which occured in some cells, and returned the tree which contains only those events.
    """
    def find_ancestor_in_list(parent,tree,list_selected):
        """
        """
        if parent in list_selected:
            return parent
        if parent==len(tree):
            return len(tree) #root
        else:
            return find_ancestor_in_list(tree[parent],tree,list_selected)

    selected_events = []
    old_indices = []
    for i,event in enumerate(df_mut_prob.index):
        event_exist_in_some_cells = (np.sum(df_mut_prob.loc[event,:]>0.7)>=4 ) or (np.sum(df_mut_prob.loc[event,:]>0.6)>15 )
        if event_exist_in_some_cells:
            selected_events.append(event)
            old_indices.append(i)

    new_tree = []
    for i in old_indices:
        parent_oldindex = find_ancestor_in_list(tree[i],tree,old_indices)
        for j in range(len(old_indices)):
            if parent_oldindex ==old_indices[j]:
                new_tree.append(j)
        if parent_oldindex==len(tree):
            new_tree.append(len(old_indices))

    return selected_events,new_tree