import pandas as pd
import numpy as np
import scipy.stats as st
from statsmodels.stats.weightstats import ztest


first_cell_type_expressions_path = input('Path to first cell expression values: ')
second_cell_type_expressions_path = input('Path to second cell expression values: ')
save_results_table = input('Path to output file: ')

first_table = pd.read_csv(first_cell_type_expressions_path, index_col = 0)
second_table = pd.read_csv(second_cell_type_expressions_path, index_col = 0)


def check_intervals_intersect(first_ci, second_ci):   
    are_intersect = False
    left_first, right_first = first_ci 
    left_second, right_second = second_ci
    
    if left_first <= right_second and left_second <= right_first:
        are_intersect = True
    else:
        are_intersect = False
    return are_intersect


def check_dge_with_ci(first_table, second_table):
    ci_test_results = []
    all_genes = list(first_table.columns[:-1])
    for gene in all_genes:
        first_ci = st.t.interval(alpha=0.95,
              df=len(first_table[gene]) - 1,
              loc=np.mean(first_table[gene]),
              scale=st.sem(first_table[gene]))
        second_ci = st.t.interval(alpha=0.95,
              df=len(second_table[gene]) - 1,
              loc=np.mean(second_table[gene]),
              scale=st.sem(second_table[gene]))
        intersect = check_intervals_intersect(first_ci, second_ci)
        ci_test_results.append(intersect)
    return ci_test_results


def check_dge_with_ztest(first_table, second_table):
    z_test_results = []
    all_genes = list(first_table.columns[:-1])
    for gene in all_genes:
        _, p_value = ztest(
            first_table[gene],
            second_table[gene]
            )
        if p_value < 0.05:
            z_result = True
        else:
            z_result = False
        z_test_results.append(z_result)
    return z_test_results

def check_pvalues(first_table, second_table):
    z_test_pval = []
    all_genes = list(first_table.columns[:-1])
    for gene in all_genes:
        stat, p_value = ztest(
            first_table[gene],
            second_table[gene]
            )
        
        z_test_pval.append(p_value)
    return z_test_pval


def mean_diff(first_table, second_table):
    mean_diff = []
    all_genes = list(first_table.columns[:-1])
    for gene in all_genes:
        mean_first = np.mean(first_table[gene])
        mean_second = np.mean(second_table[gene])
        mean_diff.append(mean_second - mean_first)
    return mean_diff



results = {
    'Genes': list(first_table.columns[:-1]),
    'ci-test_results': check_dge_with_ci(first_table, second_table),
    'z-test_results': check_dge_with_ztest(first_table, second_table),
    'z-test_p-value'": check_pvalues(first_table, second_table),
    'mean_diff': mean_diff(first_table, second_table)
}


results = pd.DataFrame(results)
results.head()

results.to_csv(save_results_table)
