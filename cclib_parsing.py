from energy_cal import *
def cclib_dic(dimer_dic,job):
    """
    Function takes a dictionary of the dimer and the job name to parse the job file with cclib,
    then creates a dictionary where the keys are the names of the dimers, and the values are cclib objects.

    dimer_dic: dictionary (not .json)
    job: str, name of the job (e.g. 'pop_opt_n','ver_a')

    ## improvement: add a progress bar.
    """
    ccdata_list = []

    for name in dimer_dic.keys():

        data = log_data(name,job)
        ccdata_list.append(data)

    cclib_dic = { k : v for k, v in zip(dimer_dic.keys(),ccdata_list) }

    return cclib_dic
