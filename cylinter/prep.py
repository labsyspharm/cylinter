import sys
import argparse
import pathlib
from subprocess import call


def main(argv=sys.argv):

    parser = argparse.ArgumentParser(
        description='Prepare an input directory for CyLinter analysis.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '-t', action='store_true', help='optional flag for TMA data'
        )
    parser.add_argument(
        'source_dir', type=str,
        help='path to mcmicro output directory'
    )
    parser.add_argument(
        'dest_dir', type=path_resolved,
        help='path to CyLinter input directory'
    )
    args = parser.parse_args()

    call([f'sh {sys.prefix}/prep_subprocess.sh {args.t} {args.source_dir} {args.dest_dir} {sys.prefix}/config.yml'], shell=True)

    return 0


def path_resolved(path_str):
    """Return a resolved Path for a string."""
    path = pathlib.Path(path_str)
    path = path.resolve()
    return path
