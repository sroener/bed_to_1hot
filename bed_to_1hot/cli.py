"""Console script for bed_to_1hot."""

import sys
import click

from bed_to_1hot.bed_to_1hot import bed_to_1hot


@click.command()
@click.option(
    "-i",
    "--bed-file",
    "bed_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input bed file with labeled regions",
)
@click.option(
    "-r",
    "--reference",
    "reference",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input genome reference in fasta format",
)
@click.option(
    "-o",
    "--output_file",
    "output_file",
    required=True,
    type=click.Path(writable=True),
    help="Path to output file.",
)
@click.option(
    "-n",
    "--label_num",
    type=click.INT,
    default=1,
    help="Number of separate label columns.",
)
@click.option(
    "--v_holdout",
    type=click.STRING,
    help="Holdout chromosomes for validation set. Expected format: chr1",
)
@click.option(
    "--t_holdout",
    type=click.STRING,
    help="Holdout chromosomes for test set. Expected format: chr1 or chr1,chr2",
)
def cli(bed_file, output_file, reference, label_num, v_holdout, t_holdout):
    """Console script for bed_to_1hot_hd5."""

    t_holdout = t_holdout.split(",")
    v_holdout = v_holdout.split(",")

    bed_to_1hot(
        input_file=bed_file,
        output_file=output_file,
        reference=reference,
        label_num=label_num,
        v_holdout=v_holdout,
        t_holdout=t_holdout
    )
    # print(type(bed_file))
    # print(bed_file, output_file, reference, label_num, v_holdout, t_holdout)


if __name__ == "__main__":
    sys.exit(cli())  # pragma: no cover
