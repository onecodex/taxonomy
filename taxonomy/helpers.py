import hashlib
import gzip
import json
import logging
import os
import platform
import sys


logger = logging.getLogger('taxonomy')
logger.setLevel(logging.INFO)
ch = logging.StreamHandler(sys.stdout)
ch.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
logger.addHandler(ch)


def read_json(f):
    if isinstance(f, dict):
        return f
    try:
        if isinstance(f, (file, gzip.GzipFile)):
            return json.load(f)
    except NameError:  # "file" doesn't exist in python 3
        from io import IOBase
        if isinstance(f, (IOBase, gzip.GzipFile)):
            return json.load(f)
    if os.path.splitext(f)[1] == '.gz':
        input_json = json.load(gzip.open(f, mode='r'))
    else:
        input_json = json.load(open(f, mode='r'))
    return input_json


def md5_for_file(filename, block_size=2**20):
    md5 = hashlib.md5()
    with open(filename, mode='rb') as file_handle:
        while True:
            data = file_handle.read(block_size)
            if not data:
                break
            md5.update(data)
    return md5.hexdigest()


def creation_time(file_path):
    """
    Try to get the date that a file was created, falling back to when it was
    last modified if that isn't possible.
    See http://stackoverflow.com/a/39501288/1709587 for explanation.
    """
    if platform.system() == 'Windows':
        return os.path.getctime(file_path)
    else:
        stat = os.stat(file_path)
        try:
            return stat.st_birthtime
        except AttributeError:
            # We're probably on Linux. No easy way to get creation dates here,
            # so we'll settle for when its content was last modified.
            return stat.st_mtime
