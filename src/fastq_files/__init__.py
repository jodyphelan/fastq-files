"""Package to work with fastq files"""

import os
import subprocess as sp
import re
from collections import defaultdict
from typing import List
import hashlib
import logging

def get_md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

__version__ = "0.0.2"

class SingleSample:
    """
    Class to represent a sample with a single fastq file
    
    Parameters
    ----------
    prefix : str
        Sample ID
    r1 : str
        Path to the R1 file
    
    Attributes
    ----------
    prefix : str
        Sample ID
    r1 : str
        Path to the R1 file
    """
    def __init__(self,prefix: str,r1: List[str]) -> None:
        self.prefix = prefix
        self.r1 = r1
        self.md5r1 = []



    def __repr__(self) -> str:
        return f"Sample(prefix={self.prefix},r1={self.r1})"

    def combine_reads(self):
        """
        Function to combine R1 and R2 files into a single file
        """
        sp.call(f"ln -s {self.r1} {self.prefix}.fastq.gz",shell=True)

    def calculate_md5(self):
        """
        Function to calculate md5 checksums for R1 and R2 files
        """
        for r in self.r1:
            md5sum_file = f"{r}.md5"
            if not os.path.exists(md5sum_file):
                md5 = get_md5(r)
                with open(md5sum_file,"w") as f:
                    f.write(f"{md5}")
            else:
                with open(md5sum_file,"r") as f:
                    md5 = f.read().strip()
            self.md5r1.append(md5)

class PairedSample:
    """
    Class to represent a sample with R1 and R2 files
    
    Parameters
    ----------
    prefix : str
        Sample ID
    r1 : List[str]
        List of R1 files
    r2 : List[str]
        List of R2 files
    
    Attributes
    ----------
    prefix : str
        Sample ID
    r1 : List[str]
        List of R1 files
    r2 : List[str]
        List of R2 files
    multiple : bool
        Whether there are multiple R1 and R2 files
    """
    def __init__(self,prefix: str,r1: List[str],r2: List[str]) -> None:
        self.prefix = prefix
        self.r1 = sorted(r1)
        self.r2 = sorted(r2)
        self.md5r1 = []
        self.md5r2 = []
        
        if len(r1)!=len(r2):
            raise ValueError("Number of R1 and R2 files do not match for sample %s" % prefix)
        if len(r1)==1:
            self.multiple = False
        else:
            self.multiple = True

    def __repr__(self) -> str:
        return f"Sample(prefix={self.prefix},r1={self.r1},r2={self.r2})"

    def combine_reads(self):
        """
        Function to combine R1 and R2 files into a single file
        """
        if len(self.r1)>1:
            sp.call(f"cat {' '.join(self.r1)} > {self.prefix}_1.fastq.gz",shell=True)
            sp.call(f"cat {' '.join(self.r2)} > {self.prefix}_2.fastq.gz",shell=True)
        else:
            sp.call(f"ln -s {self.r1[0]} {self.prefix}_1.fastq.gz",shell=True)
            sp.call(f"ln -s {self.r2[0]} {self.prefix}_2.fastq.gz",shell=True)
    
    def calculate_md5(self):
        """
        Function to calculate md5 checksums for R1 and R2 files
        """
        for r in self.r1+self.r2:
            md5sum_file = f"{r}.md5"
            if not os.path.exists(md5sum_file):
                md5 = get_md5(r)
                with open(md5sum_file,"w") as f:
                    f.write(f"{md5}")
            else:
                with open(md5sum_file,"r") as f:
                    md5 = f.read().strip()
            if r in self.r1:
                self.md5r1.append(md5)
            if r in self.r2:
                self.md5r2.append(md5)



def get_single_samples(files: List[str], r1_suffix: str) -> List[SingleSample]:
    """
    Function to sort out single files from a list of files

    Parameters
    ----------
    files : List[str]
        List of files to sort out
    r1_suffix : str
        Suffix for R1 files

    Returns
    -------
    List[Sample]
        List of Sample objects
    """
    prefixes = defaultdict(list)

    print(files)

    for f in files:
        tmp1 = re.match("%s$" % r1_suffix,f)
        p = None
        if tmp1:
            p = tmp1.group(1).split("/")[-1]
            prefixes[p].append(f)

    runs = []
    for p,vals in prefixes.items():
        vals.sort()
        
        runs.append(
            SingleSample(p,vals)
        )

    return runs


def get_paired_samples(files: List[str], r1_suffix: str,r2_suffix: str) -> List[PairedSample]:
    """
    Function to sort out paired files from a list of files

    Parameters
    ----------
    files : List[str]
        List of files to sort out
    r1_suffix : str
        Suffix for R1 files
    r2_suffix : str
        Suffix for R2 files

    Returns
    -------
    List[Sample]
        List of Sample objects
    """
    prefixes = defaultdict(lambda:{"r1":[],"r2":[]})

    for f in files:
        tmp1 = re.match("%s$" % r1_suffix,f)
        tmp2 = re.match("%s$" % r2_suffix,f)
        p = None
        if tmp1:
            p = tmp1.group(1).split("/")[-1]
            prefixes[p]['r1'].append(f)
        elif tmp2:
            p = tmp2.group(1).split("/")[-1]
            prefixes[p]['r2'].append(f)

    runs = []
    for p,vals in prefixes.items():

        if len(vals['r1'])!=len(vals['r2']):
            raise ValueError(f"Number of R1 and R2 files for sample {p} do not match")
        vals['r1'].sort()
        vals['r2'].sort()
        runs.append(
            PairedSample(p,vals['r1'],vals['r2'])
        )
    return runs



def find_paired_fastq_samples(directories: List[str], r1_pattern: str, r2_pattern:str) -> List[PairedSample]:
    """
    Find fastq files in a directory and return a
    list of tuples with the sample name and the
    path to the fastq files from both pairs.

    Parameters
    ----------
    directories : List[str]
        List of directories to search for fastq files
    r1_pattern : str
        Regex pattern for R1 files
    r2_pattern : str
        Regex pattern for R2 files

    Returns
    -------
    List[Sample]
        List of Sample objects
    """
    files = []
    for d in directories:
        for a,b,c in os.walk(d):
            for f in c:
                files.append(f"{os.path.abspath(a)}/{f}")
    fastq_files = get_paired_samples(files,r1_pattern,r2_pattern)

    return fastq_files

def find_single_fastq_samples(directories: List[str], r1_pattern: str) -> List[SingleSample]:
    """
    Find fastq files in a directory and return a
    list of tuples with the sample name and the
    path to the fastq files from both pairs.

    Parameters
    ----------
    directories : List[str]
        List of directories to search for fastq files
    r1_pattern : str
        Regex pattern for R1 files

    Returns
    -------
    List[Sample]
        List of Sample objects
    """
    files = []
    for d in directories:
        print(d)
        for a,b,c in os.walk(d):
            print(a,b,c)
            for f in c:
                files.append(f"{os.path.abspath(a)}/{f}")
    fastq_files = get_single_samples(files,r1_pattern)

    return fastq_files



# def hamming(s1: str, s2: str) -> int:
#     return sum(c1 != c2 for c1, c2 in zip(s1, s2))

# def guess_paired_files(files: List[str]) -> List[Sample]:
#     pairs = set()
#     for fi in files:
#         for fj in files:
#             if fi != fj:
#                 if hamming(fi,fj)==1:
#                     pairs.add(tuple(sorted([fi,fj])))

#     for pair in pairs:
