#!/usr/bin/env python

"""Tests for `bed_to_1hot` package."""

# pylint: disable=redefined-outer-name

import pytest
import numpy as np
import pandas as pd
import h5py

from bed_to_1hot import bed_to_1hot


@pytest.fixture(scope="module")
def dna_reference_location():
    """returns location of referencen fasta"""
    return "/home/user/git_private/data/reference/hs38.fa"


@pytest.fixture(scope="module")
def full_bed_location():
    """returns path to example full spec bed file"""
    return "tests/test_data/full_bed.bed"


@pytest.fixture(scope="module")
def full_spec_bed(full_bed_location):
    """returns dataframe with information from full spec bed file"""
    yield bed_to_1hot.read_bed_file(full_bed_location)


@pytest.fixture(scope="module")
def example_bed():
    """yields pandas dataframe with example bed content"""
    yield pd.read_csv("tests/test_data/encoding_test.bed", sep="\t", header=None)


@pytest.fixture(scope="module")
def example_bed_ref_array():
    """returns one hot encoded reference array for example bed"""
    return np.load("tests/test_data/example_bed_ref_array.npy")


@pytest.fixture(scope="module")
def dna_sequence():
    """returns DNA sequence"""
    return "ACGTYRSWKMBDHVN-"


@pytest.fixture(scope="module")
def dna_sequence_ref_array():
    """returns one hot encoded reference array for dna_sequence"""
    return np.load("tests/test_data/dna_sequence_ref_array.npy")


@pytest.fixture(scope="module")
def example_bed_l2():
    """yields pandas dataframe with example bed content"""
    yield bed_to_1hot.read_bed_file(
        "tests/test_data/encoding_test_2label.bed", labelnum=2
    )


@pytest.fixture(scope="module")
def example_bed_l2_h5():
    """yields HD5 instance with train, validation and test data for example_bed_l2"""
    yield h5py.File("tests/test_data/example_test_2label.h5", "r")


def test_csv_reader_header_fields(full_spec_bed):
    """
    Happy Path test to make sure the processed data contains the right header fields"""
    data = full_spec_bed
    print(data.head())
    header_fields = list(data.columns)
    assert header_fields == [
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


def test_bed_reader_data_contents(full_spec_bed):
    """
    Happy Path Test to examine that each row
    had the appropriate data type per field
    """

    data = full_spec_bed

    # Check row types
    for _, row in data.iterrows():
        assert isinstance(row["chrom"], str)
        assert isinstance(row["chromStart"], int)
        assert isinstance(row["chromEnd"], int)
        assert isinstance(row["name"], str)
        assert isinstance(row["score"], int)
        assert isinstance(row["strand"], str)
        assert isinstance(row["thickStart"], int)
        assert isinstance(row["itemRGB"], int)
        assert isinstance(row["blockCount"], int)
        assert isinstance(row["blockSizes"], str)
        assert isinstance(row["blockStarts"], str)

    # Basic data checks
    assert len(data) == 100  # We have collected 100 rows
    # first line
    assert data.iloc[0]["chrom"] == "chr21"
    assert data.iloc[0]["chromStart"] == 9928613
    assert data.iloc[0]["name"] == "uc002yip.1"
    # last line
    assert data.iloc[-1]["chrom"] == "chr21"
    assert data.iloc[-1]["chromStart"] == 26018661
    assert data.iloc[-1]["name"] == "uc002ylw.1"


def test_one_hot_encoding(dna_sequence, dna_sequence_ref_array):
    """tests proper one hot encoding of IUPAC DNA sequences"""

    upper = dna_sequence.upper()
    lower = dna_sequence.lower()

    ref_array = dna_sequence_ref_array

    np.testing.assert_array_equal(bed_to_1hot.one_hot_encoding(upper), ref_array)
    np.testing.assert_array_equal(bed_to_1hot.one_hot_encoding(lower), ref_array)


def test_bed_to_1hot(example_bed, dna_reference_location, example_bed_ref_array):
    """test DNA encoding from bed file"""

    bed = example_bed
    ref = dna_reference_location
    ref_array = example_bed_ref_array
    encoded_array = bed_to_1hot.bed_encoding(bed, ref)
    # test array dimensions
    assert len(encoded_array) == 10
    for seq in encoded_array:
        assert len(seq) == 10
    # test encoding
    np.testing.assert_array_equal(encoded_array, ref_array)


@pytest.mark.parametrize(
    "v_holdout,t_holdout",
    [
        ("chr1,chr5", "chr9"),
    ],
)
def test_train_validate_test_split(
    example_bed_l2, t_holdout, v_holdout, dna_reference_location, example_bed_l2_h5
):
    """test dataset split in training, validation and test set"""
    bed_df = example_bed_l2
    ref = dna_reference_location
    t_holdout = t_holdout.split(",")
    v_holdout = v_holdout.split(",")
    data = example_bed_l2_h5
    (
        x_train,
        y_train,
        x_val,
        y_val,
        x_test,
        y_test,
    ) = bed_to_1hot.train_validate_test_split(bed_df, v_holdout, t_holdout, ref)
    # test array dimensions
    assert x_train.shape == (6, 10, 4)
    assert y_train.shape == (6, 2)
    assert x_val.shape == (2, 10, 4)
    assert y_val.shape == (2, 2)
    assert x_test.shape == (2, 10, 4)
    assert y_test.shape == (2, 2)
    # test array content
    np.testing.assert_array_equal(np.array(data["train_data"]["x_train"]), x_train)
    np.testing.assert_array_equal(np.array(data["train_data"]["y_train"]), y_train)
    np.testing.assert_array_equal(np.array(data["val_data"]["x_val"]), x_val)
    np.testing.assert_array_equal(np.array(data["val_data"]["y_val"]), y_val)
    np.testing.assert_array_equal(np.array(data["test_data"]["x_test"]), x_test)
    np.testing.assert_array_equal(np.array(data["test_data"]["y_test"]), y_test)
