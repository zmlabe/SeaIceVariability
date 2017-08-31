"""
Function calculates statistical significance
 
Notes
-----
    Author : Zachary Labe
    Date   : 13 August 2017
    
Usage
-----
    [1] calc_indttest(varx,vary)
"""
    
def calc_indttest(varx,vary):
    """
    Function calculates statistical difference for 2 related
    sample t-test

    Parameters
    ----------
    varx : 3d array
    vary : 3d array
    
    Returns
    -------
    stat = calculated t-statistic
    pvalue = two-tailed p-value

    Usage
    -----
    stat,pvalue = calc_ttest(varx,vary)
    """
    print('\n>>> Using calc_ttest function!')
    
    ### Import modules
    import numpy as np
    import scipy.stats as sts
    
    ### 2-independent sample t-test
    stat,pvalue = sts.ttest_rel(varx,vary,nan_policy='omit')
    
    ### Significant at 95% confidence level
    pvalue[np.where(pvalue >= 0.05)] = np.nan
    pvalue[np.where(pvalue < 0.05)] = 1.
    
    print('*Completed: Finished calc_ttest function!')
    return stat,pvalue
