import sys
import argparse
import pathlib
import logging
# import pandas as pd
from .config import Config
from . import pipeline, components

logger = logging.getLogger(__name__)


def main(argv=sys.argv):

    epilog = 'Pipeline modules:\n'
    epilog += '\n'.join(f"    {n}" for n in components.pipeline_module_names)
    parser = argparse.ArgumentParser(
        description='Perform CyLinter analysis on a data file.',
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        'config', type=path_resolved,
        help='Path to the configuration YAML file'
    )
    parser.add_argument(
        '--module', type=str,
        help='Pipeline module at which to begin processing (see below'
        ' for ordered list of modules)'
    )
    args = parser.parse_args(argv[1:])
    if not validate_paths(args):
        return 1
    if args.module and args.module not in components.pipeline_module_names:
        print(
            f"cylinter: error: argument --module: invalid choice '{args.module}'",
            file=sys.stderr
        )
        return 1

    logging.basicConfig(
        level=logging.INFO,
        format='%(levelname)s: %(message)s'
    )

    logger.info("Reading configuration file")
    config = Config.from_path(args.config)
    create_output_directory(config)

    logger.info("Executing pipeline")
    pipeline.run_pipeline(config, args.module)

    logger.info("Finished")

    return 0


def path_resolved(path_str):
    """Return a resolved Path for a string."""
    path = pathlib.Path(path_str)
    path = path.resolve()
    return path


def validate_paths(args):
    """Validate the Path entries in the argument list."""
    ok = True
    if not args.config.exists():
        print(
            f"Config path does not exist:\n     {args.config}\n",
            file=sys.stderr
        )
        ok = False
    return ok


def create_output_directory(config):
    """Create the output directory structure given the configuration object."""
    config.out_dir.mkdir(parents=True, exist_ok=True)
