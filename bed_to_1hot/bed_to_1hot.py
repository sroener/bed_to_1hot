"""Main module."""

import pandas as pd
import numpy as np
from pyfaidx import Fasta
import h5py

BEDCOLS = [
    "chrom",
    "chromStart",
    "chromEnd",
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "itemRGB",
    "blockCount",
    "blockSizes",
    "blockStarts",
]


def generate_colnames(df, labelnum=0):  # need to be adjusted for GC content
    """Generates column names for an input dataframe.

    Arguments:
        df [dataframe] -- [dataframe for which column names are generated]

    Returns:
        [list] -- [list of generated column names]
    """
    colnames = []
    for field in range(len(df.columns) - labelnum):
        colnames.append(BEDCOLS[field])
    for label in range(labelnum):
        colnames.append(f"label_{label+1}")
    return colnames


def read_bed_file(path, labelnum=0):
    """reads .bed file into pandas dataframe

    Args:
        path (str): path to bed file

    Returns:
        [pandas dataframe]: dataframe containing all fields contained in bed file
    """

    bed_df = pd.read_table(path, sep="\t", header=None)
    colnames = generate_colnames(bed_df, labelnum)
    bed_df.columns = colnames
    return bed_df


def one_hot_encoding(sequence):
    """One hot encoding of a DNA sequence.

    Args:
        sequence (str): [DNA sequence in IUPAC format]

    Returns:
        [numpy array]: [one hot encoded DNA sequence]
    """

    mydict = {
        "A": np.asarray([1, 0, 0, 0]),
        "a": np.asarray([1, 0, 0, 0]),
        "C": np.asarray([0, 1, 0, 0]),
        "c": np.asarray([0, 1, 0, 0]),
        "G": np.asarray([0, 0, 1, 0]),
        "g": np.asarray([0, 0, 1, 0]),
        "T": np.asarray([0, 0, 0, 1]),
        "t": np.asarray([0, 0, 0, 1]),
        "Y": np.asarray([0, 1, 0, 1]),
        "y": np.asarray([0, 1, 0, 1]),
        "R": np.asarray([1, 0, 1, 0]),
        "r": np.asarray([1, 0, 1, 0]),
        "S": np.asarray([0, 1, 1, 0]),
        "s": np.asarray([0, 1, 1, 0]),
        "W": np.asarray([1, 0, 0, 1]),
        "w": np.asarray([1, 0, 0, 1]),
        "K": np.asarray([0, 0, 1, 1]),
        "k": np.asarray([0, 0, 1, 1]),
        "M": np.asarray([1, 1, 0, 0]),
        "m": np.asarray([1, 1, 0, 0]),
        "B": np.asarray([0, 1, 1, 1]),
        "b": np.asarray([0, 1, 1, 1]),
        "D": np.asarray([1, 0, 1, 1]),
        "d": np.asarray([1, 0, 1, 1]),
        "H": np.asarray([1, 1, 0, 1]),
        "h": np.asarray([1, 1, 0, 1]),
        "V": np.asarray([1, 1, 1, 0]),
        "v": np.asarray([1, 1, 1, 0]),
        "N": np.asarray([0, 0, 0, 0]),
        "n": np.asarray([0, 0, 0, 0]),
        "-": np.asarray([0, 0, 0, 0]),
    }

    nuc_list = list()
    for nuc in list(sequence):
        nuc_list.append(mydict[nuc])
    result = np.stack(np.asarray(nuc_list, dtype="int8"))
    return result


def bed_encoding(bed_df, reference):
    """One Hot encoding of regions in a bed file.

    Args:
        bed_df (pandas Dataframe): [Dataframe of bed file]
        reference (fasta file): [Genome reference file in fasta format]

    Returns:
        [numpy array]: [One hot encoded DNA sequences from bed file]
    """

    fasta = Fasta(reference, as_raw=True)
    seq_list = list()
    for _, i in bed_df.iterrows():
        seq_list.append(one_hot_encoding(fasta[i[0]][i[1] : i[2]]))
    result = np.stack(seq_list)
    return result


def train_validate_test_split(bed_df, v_holdout, t_holdout, ref):
    """[summary]

    Args:
        bed_df ([type]): [Dataframe of bed file]
        v_holdout ([type]): [Holdout chromosomes for validation]
        t_holdout ([type]): [Holdout chromosomes for testing]
        ref ([type]): [Genome reference file in fasta format]

    Returns:
        [numpy arrays]: [train, validate, test data with data(x) and label(y) for each subset]
    """
    train_df = bed_df.loc[~bed_df["chrom"].isin(t_holdout + v_holdout)]
    validate_df = bed_df.loc[bed_df["chrom"].isin(v_holdout)]
    test_df = bed_df.loc[bed_df["chrom"].isin(t_holdout)]
    x_train = bed_encoding(train_df, ref)
    y_train = np.asarray(
        train_df.loc[:, train_df.columns.str.contains("label")], dtype="int8"
    )  # .flatten()
    x_val = bed_encoding(validate_df, ref)
    y_val = np.asarray(
        validate_df.loc[:, validate_df.columns.str.contains("label")], dtype="int8"
    )  # .flatten()
    x_test = bed_encoding(test_df, ref)
    y_test = np.asarray(
        test_df.loc[:, test_df.columns.str.contains("label")], dtype="int8"
    )  # .flatten()
    return x_train, y_train, x_val, y_val, x_test, y_test


def save_to_hd5(out_file, x_train, y_train, x_val, y_val, x_test, y_test):
    """Writes train,validate,test arrays to HD5 format"""
    data = h5py.File(out_file, "w")
    train_data = data.create_group("train_data")
    train_data.create_dataset("x_train", data=x_train)
    train_data.create_dataset("y_train", data=y_train)
    val_data = data.create_group("val_data")
    val_data.create_dataset("x_val", data=x_val)
    val_data.create_dataset("y_val", data=y_val)
    test_data = data.create_group("test_data")
    test_data.create_dataset("x_test", data=x_test)
    test_data.create_dataset("y_test", data=y_test)
    data.close()


def bed_to_1hot(input_file, output_file, reference, label_num, v_holdout, t_holdout):
    """Takes a bedfile, queries a fasta file and saves the one hot encoded sequences to hd5."""
    bed_df = read_bed_file(input_file, label_num)
    # print("generating data split")
    x_train, y_train, x_val, y_val, x_test, y_test = train_validate_test_split(
        bed_df=bed_df, ref=reference, v_holdout=v_holdout, t_holdout=t_holdout
    )
    save_to_hd5(
        out_file=output_file,
        x_train=x_train,
        y_train=y_train,
        x_val=x_val,
        y_val=y_val,
        x_test=x_test,
        y_test=y_test,
    )
