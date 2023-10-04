import click
from targeted import (
    align_reads,
    filter_BAM_plots,
    make_BAMs,
    make_indexes,
    make_reference_sequences,
    parse_BAMs,
    switch_loc_in_ref,
)


@click.group()
def main():
    pass


main.add_command(align_reads.main)
main.add_command(filter_BAM_plots.main)
main.add_command(make_BAMs.main)
main.add_command(make_indexes.main)
main.add_command(make_reference_sequences.main)
main.add_command(parse_BAMs.main)
main.add_command(switch_loc_in_ref.main)
