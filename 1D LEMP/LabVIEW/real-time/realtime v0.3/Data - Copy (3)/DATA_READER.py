# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 16:37:21 2024

@author: Alexander
"""
import numpy as np
import os


# def file_next(base_name, extension, increment=1, start=1, flist=set()):
#     index = start
#     while True:
#         file_name = f"{base_name}{index}.{extension}"
#         if os.path.exists(file_name):
#             flist.add(file_name)
#         else:
#             return flist
#         index += increment

def list_files_in_directory(directory):
    files = []
    for filename in os.listdir(directory):
        if os.path.isfile(os.path.join(directory, filename)):
            files.append(filename)
    return files

current_directory = os.getcwd()
files = list_files_in_directory(current_directory)
print("Files in directory:")
for file in files:
    print(file)