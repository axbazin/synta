#!/usr/bin/env python3
#coding:utf-8

import logging
import gzip
from io import TextIOWrapper

def is_compressed(file_or_file_path):
    """
        Checks is a file, or file path given is compressed or not
    """
    file = file_or_file_path
    if type(file) == str:
        file = open(file, "rb")
    else:
        try:
            file = open(file.name, "rb")
        except:
            return False
    if file.read(2).startswith(b'\x1f\x8b'):
        return True
    file.close()
    return False


def read_compressed_or_not(file_or_file_path):
    """
        Reads a file, compressed or not.
        Copied from http: //www.github.com/ggautreau/PPanGGOLiN.git's utils.py.
    """
    file = file_or_file_path
    if isinstance(file, str):
        file = open(file, "rb")
    else:
        try:
            file = open(file.name, "rb")
        except AttributeError:
            return(file)
    if file.read(2).startswith(b'\x1f\x8b'):
        file.seek(0)
        logging.getLogger().info("Uncompressing the file: '" + file.name + "' ...")
        return(TextIOWrapper(gzip.open(filename=file, mode="r")))
    else:
        file.close()
        file = open(file.name, "r")
        return(file)


def write_compress_or_not(file_path, compress):
    """
        Returns a file-like object, compressed or not.
    """
    if compress:
        return gzip.open(file_path + ".gz", mode="wt")
    else:
        return open(file_path, "w")

if __name__ == "__main__":
    pass