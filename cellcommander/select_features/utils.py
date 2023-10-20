import numpy as np
import logging

logger = logging.getLogger("cellcommander")

def is_counts(data):
    if np.all(data >= 0) and np.all(data.astype(int) == data):
        logger.info("The matrix contains count data.")
        return True
    else:
        logger.info("The matrix does not contain count data.")
        return False

def is_mostly_counts(data, percent=0.9):
    """Want to check if some percent of the data is counts

    Args:
        data (_type_): _description_
    """
    if np.all(data >= 0) and np.all(data.astype(int) == data):
        logger.info("The matrix contains all count data.")
        return True
    elif np.sum(data >= 0) / data.size >= percent and np.sum(data.astype(int) == data) / data.size >= percent:
        greater_than_0 = (np.sum(data >= 0) / data.size)*100
        int_equals = (np.sum(data.astype(int) == data) / data.size)*100
        logger.info(f"The matrix contains mostly count data. {greater_than_0}% of the data is greater than 0 and {int_equals}% of the data is equal to its integer value.")
        return True
    else:
        logger.info("The matrix does not contain mostly count data.")
        return False