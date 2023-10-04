import click

from glob import glob
from pathlib import Path
from subprocess import run, PIPE


@click.command()
@click.option(
    "-i",
    "--ref-dir",
    required=True,
    help="Directory with fasta files to build indexes. Fasta files must have the .fna extension.",
)
@click.option(
    "-p",
    "--threads",
    default=4,
    show_default=True,
    help="Number of threads to run hisat2-build with.",
)
def main(ref_dir, threads):
    arr = glob(f"{ref_dir}/*.fna")

    for fasta in arr:
        fasta_name = Path(fasta).name[:-4]
        index_dir = f"{fasta[:-4]}_index"

        Path(index_dir).mkdir(exist_ok=True, parents=True)

        run(
            [
                "hisat2-build",
                "-p",
                f"{threads}",
                f"{fasta}",
                f"{index_dir}/{fasta_name}_index",
            ],
            stdout=PIPE,
            stderr=PIPE,
        )
