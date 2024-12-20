#!/usr/bin/env python

# Disable false positive/opinionated lints
# pylint: disable = invalid-name, consider-using-f-string

"""
CRISPR Pooled Screen YAML Metadata Validator
--------------------------------------------
Takes the YAML output from one of the metadata parsers (such as the
Illumina Sample Sheet parser, with an appropriate schema) -- or an
otherwise suitable YAML input -- with the sample FASTQs to validate and
transform it into a form that the CRISPR Pooled Screen Nextflow pipeline
expects.
"""

# NOTE This is one big script, rather than being nicely modularised, to
# satisfy the constraints imposed by Nextflow. (With apologies to future
# maintainers!)

import argparse
import logging
import os
import re
import sys
import time
from collections import Counter
from collections.abc import Hashable, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import yaml
from schema import And, Optional, Or, Schema, SchemaError, Use  # type: ignore

## INPUT YAML SCHEMA ###################################################


def _MAYBE(obj: object = str) -> Or:
    # Optional values
    return Or(obj, None)


def _CAST(to: Callable[[str], Any] | None = None) -> And:
    # Cast a string input into something else
    return And(str, Use(to or (lambda x: x)))


def _VALUE_UNIT_SCHEMA(to: Callable[[str], Any] = str):
    # Run the value through a cast function
    return {
        "value": _CAST(to),
        "unit": _MAYBE(),
    }


def _maybe_numeric(value: str | None) -> int | float | str | None:
    # Attempt to cast a value into a numeric type, starting with
    # integers before moving to floating point, otherwise pass-through
    # the original value; unless it's empty or None, in which case,
    # return None
    if value is None or value == "":
        return value

    try:
        return int(value)

    except ValueError:
        try:
            return float(value)

        except ValueError:
            return value


_EXPERIMENT_SCHEMA = {
    "experiment_name": str,
    "experiment_date": str,
    "user": str,
    "point_of_contact": _MAYBE(),
    "biology_poc": _MAYBE(),
    "analysis_poc": _MAYBE(),
    "other_poc": _MAYBE(),
    "project_name": _MAYBE(),
    "project_code": _MAYBE(),
    "department": _MAYBE(),
    "guide_library": _MAYBE(),
    Optional(str): object,  # Everything else
}


_SAMPLE_SCHEMA = {
    "aliquot": str,
    And("condition", Use(lambda _: "conditions")): {
        str: _VALUE_UNIT_SCHEMA(_maybe_numeric),
    },
    "group": str,
    "replicate": _CAST(_maybe_numeric),
    "representation": _CAST(int),
    "sample": str,
    "reference": _MAYBE(
        {
            str: And(
                # We just want the value; the unit is discarded
                _VALUE_UNIT_SCHEMA(),
                Use(lambda vu: vu["value"] or None),
            )
        },
    ),
    Optional("is_ref"): _CAST(int),
    Optional(str): object,  # Everything else
}

SCHEMA = {
    Optional("experiment", default={}): _EXPERIMENT_SCHEMA,
    "samples": [_SAMPLE_SCHEMA],
}

## SAMPLE HANDLING #####################################################

# Pattern for '{<TAG>}' template tags
_TEMPLATE = re.compile(r"{.+?}")

# Pattern for Illumina FASTQ file names
_ILLUMINA_FASTQ = re.compile(
    r"""
        (?P<name>      # 'name' capture group
            [\w-]+     # At least one alphanumeric character, underscore or dash
        )

        _S             # Sequence marker

        (?P<sequence>  # 'sequence' capture group
            \d+        # At least one digit
        )

        (?:
            _L         # Lane marker

            (?P<lane>  # 'lane' capture group
                \d+    # At least one digit
            )+
        )?             # Optional

        _R?            # Read marker

        (?P<read>      # 'read' capture group
            \d+
        )

        (?:
            _\d+       # Optional suffix
        )?

        (?:
            \.fastq    # Filename extension
        )
    """,
    re.VERBOSE,
)


@dataclass
class SampleFastq:
    """
    Illumina FASTQ file, with its parts extracted
    """

    path: Path

    name: str
    read: str
    sequence: str
    lane: str | None

    def __init__(self, fastq: Path):
        self.path = fastq

        if (match := _ILLUMINA_FASTQ.search(fastq.name)) is None:
            logging.critical("Could not parse '%s' as an Illumina FASTQ", fastq)
            sys.exit(1)

        self.name = match["name"]
        self.read = match["read"]
        self.sequence = match["sequence"]
        self.lane = match["lane"]

    def __hash__(self) -> int:
        return hash(self.path)


class SampleError(Exception):
    """
    Raised when a sample cannot be validated
    """


def get_fastqs(fastq_dir: Path) -> Iterator[SampleFastq]:
    """
    Fetch all FASTQ files in the search path, parsed into their
    component parts (per the standard Illumina filename format)

    @param   fastq_dir  FASTQ search path
    @return  Iterator of parsed FASTQ files
    """
    yield from [
        SampleFastq(fastq.relative_to(fastq_dir)) for fastq in fastq_dir.glob("**/*.fastq*")
            if fastq.is_file()
    ]


def _serialise_condition(value: dict[str, Any]) -> str:
    match value:
        case {"value": value, "unit": None}:
            return str(value)

        case {"value": value, "unit": unit}:
            return f"{value}{unit}"

        case _:
            # This should never happen
            raise NotImplementedError


def build_sample(
    sample_metadata: dict, fastqs: dict[str, set[SampleFastq]]
) -> tuple[str, dict]:
    """
    Build the sample record from the raw input

    @param   sample_metadata  Validated sample input dictionary
    @param   fastqs           FASTQs, bucketed by name
    @return  Tuple of sample name and munged output dictionary
    """
    sample = sample_metadata.pop("sample")

    # group
    group = sample_metadata.pop("group")

    if not group:
        group = "_".join(
            f"{condition}_{_serialise_condition(value)}"
            for condition, value in sample_metadata["conditions"].items()
        )

    elif _TEMPLATE.search(group):
        try:
            # FIXME In general, the template values here are not
            # necessarily going to be strings, or serialise to strings
            # without spaces, which will thus trigger the no spaces
            # failure (see below). This is presumably left to the user
            # to provide valid input.
            group = group.format(**sample_metadata)

        except KeyError as e:
            logging.critical("Group for sample '%s' expects a %s value", sample, e)
            raise SampleError from e

    if any(char.isspace() for char in group):
        logging.critical("Group for sample '%s' contains spaces: '%s'", sample, group)
        raise SampleError

    sample_metadata["group"] = group

    # group_rep
    sample_metadata["group_rep"] = "{group}.rep_{replicate}".format(**sample_metadata)

    # reads
    try:
        for fastq in fastqs[sample]:
            reads = sample_metadata.setdefault(f"read{fastq.read}", [])
            reads.append(str(fastq.path))

    except KeyError as e:
        logging.critical("No FASTQ files found for sample '%s'", sample)
        raise SampleError from e

    return sample, sample_metadata


def _merkle_hash(value: Any) -> int:
    # Extend Python's native hash functionality to dictionaries, by
    # calculating their Merkle hash. That is, the hash of the tuple of
    # (key, value) pairs, sorted by key, which themselves have been
    # hashed (recursively).
    match value:
        case dict():
            return _merkle_hash(
                tuple(
                    _merkle_hash(
                        (
                            _merkle_hash(key),
                            _merkle_hash(value[key]),
                        )
                    )
                    for key in sorted(value)
                )
            )

        case Hashable():
            return hash(value)

        case _:
            raise TypeError(f"Cannot calculate the Merkle hash of '{value}'")


def check_groups_vs_conditions(munged: dict[str, dict]) -> bool:
    """
    Check each group has a single, unique set of conditions

    @param   munged  Munged sample dictionary
    @return  Whether the check passes
    """
    # Mapping from group name to conditions hashes
    group_conditions: dict[str, set[int]] = {}

    for sample in munged.values():
        conditions = group_conditions.setdefault(sample["group"], set())
        conditions.add(_merkle_hash(sample["conditions"]))

    unique_conditions = set()
    for conditions in group_conditions.values():
        if len(conditions) != 1:
            # Each group should have a single set of conditions
            return False

        unique_conditions |= conditions

    # Each group should have a unique condition; i.e., the number of
    # groups and the number of unique conditions must be equal
    return len(group_conditions) == len(unique_conditions)


## HELPERS #############################################################


class _UTCLoggingFormatter(logging.Formatter):
    converter = time.gmtime

    def __init__(self):
        super().__init__(
            fmt="{asctime}\t{levelname}\t{message}",
            datefmt="%Y-%m-%dT%H:%M:%SZ",
            style="{",
        )


def _configure_logging():
    logger = logging.getLogger()
    logger.setLevel("DEBUG" if "DEBUG" in os.environ else "INFO")

    handler = logging.StreamHandler()
    handler.setFormatter(_UTCLoggingFormatter())
    logger.addHandler(handler)


def _parse_arguments(*args: str) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__.strip(),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "fastq_dir",
        metavar="FASTQ_DIR",
        type=Path,
        help="directory containing respective FASTQ files",
    )

    parser.add_argument(
        "--input",
        type=argparse.FileType(mode="rt"),
        default="-",
        help="input metadata (YAML; defaults to stdin)",
    )

    parser.add_argument(
        "--output",
        type=argparse.FileType(mode="xt"),
        default="-",
        help="output file (YAML; defaults to stdout)",
    )

    return parser.parse_args(args)


## ENTRYPOINT ##########################################################


def main(*argv: str) -> int:
    """Entrypoint"""
    _configure_logging()
    args = _parse_arguments(*argv)

    if not args.fastq_dir.is_dir():
        logging.critical("FASTQ directory does not exist or is not a directory")
        sys.exit(1)

    logging.info("Looking up FASTQ files from %s", args.fastq_dir)
    fastqs: dict[str, set[SampleFastq]] = {}
    for fastq in get_fastqs(args.fastq_dir):
        # Bucket FASTQs by name
        fq = fastqs.setdefault(fastq.name, set())
        fq.add(fastq)

    logging.info("Reading input YAML from %s", args.input.name)
    raw = yaml.safe_load(args.input)

    # Validate the raw input and munge it into the expected shape
    try:
        validated = Schema(SCHEMA).validate(raw)
        munged = {
            "experiment": validated["experiment"],
            "samples": {},
        }

    except SchemaError as e:
        logging.debug("Schema validation failure: %s", e)
        logging.critical("Could not validate input for pipeline")
        sys.exit(1)

    all_ok = True

    samples: Counter[str] = Counter()
    for sample in validated["samples"]:
        try:
            sample_name, sample_munged = build_sample(sample, fastqs)

            samples.update([sample_name])
            munged["samples"][sample_name] = sample_munged

        except SampleError:
            all_ok = False

    if samples.total() > len(samples):
        logging.critical("Input contains duplicate sample names:")
        for sample_name, count in samples.items():
            if count > 1:
                logging.critical("* %s", sample_name)

        all_ok = False

    if not check_groups_vs_conditions(munged["samples"]):
        logging.critical("Number of groups does not match the number of conditions")
        all_ok = False

    if not all_ok:
        sys.exit(1)

    logging.info("Writing validated pipeline metadata YAML to %s", args.output.name)
    yaml.safe_dump(munged, args.output, explicit_start=True, sort_keys=False)
    args.output.flush()

    logging.info("Pipeline metadata validation complete")
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
