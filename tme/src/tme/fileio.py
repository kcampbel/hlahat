import os
def check_paths_exist(paths):
    exists = dict()
    for ii in paths:
        if os.path.exists(ii):
            exists[ii] = True
        else:
            exists[ii] = False
    return exists

