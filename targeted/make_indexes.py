import click

from glob import glob
from pathlib import Path
from subprocess import run, PIPE
from tart.utils.helpers import print


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

    print(f"Found these reference files:")
    for i in arr:
        print(i)

    print("-------------------------------------------------")
    print("-------------------------------------------------\n")

    for fasta in arr:
        fasta_name = Path(fasta).name[:-4]
        index_dir = f"{fasta[:-4]}_index"

        Path(index_dir).mkdir(exist_ok=True, parents=True)

        print(f"Generating index for {fasta}")

        try:
            call = run(
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
            print("Success")

        except:
            print("Failed: ", call.stderr)
