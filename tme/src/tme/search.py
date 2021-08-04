import re
import os, fnmatch

def search_list(pattern, input):
    r = re.compile(pattern, flags = re.I)
    mat = filter(r.search, input)
    return(list(mat))

def search_dict(pattern, input, item = "keys"):
    newlist = list()
    r = re.compile(pattern)
    #r = re.compile(pattern, flags = re.I)
    if item == 'keys':
        for i in input.keys():
            newlist.append(i)
    if item == 'values':
        for i in input.values():
            newlist.append(i)
    mat = filter(r.search, newlist)
    return(list(mat))

def replace_element_dict(x, dictionary):
    """ Searches list for keys in dictionary and replaces matching elements with value
    """
    out = list(x)
    for pattern in dictionary.keys():
        sl = search_list(pattern, out)
        if sl:
            element = out.index(sl[0])
            out[element] = dictionary[pattern]
    return(out)

def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below
    supplied root directory.'''
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)
