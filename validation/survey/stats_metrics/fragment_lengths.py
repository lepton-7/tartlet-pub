# %%
import pysam
import pandas as pd

from glob import glob
from argparse import ArgumentParser
from tart.utils.mpi_context import BasicMPIContext

# %%

flagmask = 67

# %%

def main():
    # Determine MPI context
    mp_con = BasicMPIContext()
    comm = mp_con.comm
    rank = mp_con.rank

