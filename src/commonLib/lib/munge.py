import numpy as np
from datetime import datetime

def merge_indicator_rename(series, labels:list):
    """ Renames a pandas merge indicator column with labels """
    if len(labels) != 3:
        raise ValueError('Number of labels for renaming merge indicator must be 3')
    out = np.select(
        [series == 'left_only', series == 'right_only', series == 'both'],
        labels,
        None
    )
    return(out)

def get_timestamp():
    now = datetime.now()
    return(now.strftime("%Y%m%dT%H%M%S+%f"))
