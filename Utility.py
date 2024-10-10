import os


def is_nonzero_file(fpath):  
    # Check if a file exists and is not empty
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0






    