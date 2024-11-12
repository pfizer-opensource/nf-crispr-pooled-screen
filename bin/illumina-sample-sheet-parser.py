#!/usr/bin/env python

# Disable false positive/opinionated lints
# pylint: disable = invalid-name, too-few-public-methods, too-many-locals

"""
Illumina Sample Sheet Parser
----------------------------
Takes an Illumina Sample Sheet (CSV) and converts it to a standardised
YAML format, driven by a mapping schema file. The resulting YAML output
is then easier to handle in downstream processing.

---

Mapping Schema Format
---------------------
The mapping schema takes the form of a YAML file that maps output keys
to sections within the Illumina Sample Sheet. Sections in the input that
are not mapped to an output are discarded.

Each section has the following attributes:
  section               section name
  type                  type of CSV data represented by the section
  case-sensitive        whether section name and key/header matching
                        should be case-sensitive (defaults to false)
  delimiter             delimiter to delineate structure metadata (see
                        below; defaults to '|' character)
  required              list of required keys/headers
  optional              list of optional keys/headers

If required keys/headers are not present in the input, the process will
fail. If optional keys/headers are not present in the input, they will
be given null values in the output. Unspecified key/header mappings will
be passed into the output verbatim.

CSV data types:
  tabular               tabular CSV data, with a header as the first row
  key-value             two-column CSV data, with keys in the first
                        column and respective values in the second (all
                        other columns are ignored)

Section names, required and optional keys/headers accept both simple
strings and "names with aliases". The latter is represented in YAML as:

  name: canonical name
  aliases:
  - some alias
  - another alias

When aliases are found in the input, they are mapped to the canonical
name in the output. By default, ambiguities that result from non-unique
names/aliases are forbidden, as it can lead to conflicts and unexpected
behaviour. This can be overridden with --allow-ambiguous-schema, but you
do this at your own risk.

Keys/headers in the input, which are specified in the mapping schema,
can admit additional metadata that is used to create structured output.
(Unspecified keys/headers are not subject to this destructuring and will
be passed into the output verbatim.) The metadata is represented in the
input as a suffix to a common field name, taking the form:

  <FIELD> <DELIMITER> <PROPERTY> (<UNIT>)

A delimiter -- a pipe character, by default -- signals the presence of
metadata, with any interleaving whitespace ignored. The parenthetical
<UNIT> is optional, arbitrary text; whereas the <PROPERTY> drives the
transformation:

  []                    transform to a list
  <STRING>              transform to a key-value mapping, with the
                        property string as the key

Mixed property types, over common keys/headers, will result in failure.
"""

# NOTE This is one big script, rather than being nicely modularised, to
# satisfy the constraints imposed by Nextflow. (With apologies to future
# maintainers!)

from __future__ import annotations

import argparse
import csv
import dataclasses
import logging
import os
import re
import sys
import time
from abc import ABC, abstractmethod
from collections.abc import Iterator, Mapping
from enum import Enum
from io import StringIO
from typing import TextIO, TypedDict, cast

import schema  # type: ignore
import yaml

## CONSTANTS ###########################################################

# Pattern for Illumina Sample Sheet INI-style section markers
_SECTION_MARKER = re.compile(
    r"""
        (?<= ^ \[ )  # Anchor to opening square brace at start of string

        [A-Za-z\s]*  # Any number of alphabetic letters or whitespace

        (?= \] $ )   # Anchor to closing square brace at end of string
    """,
    re.VERBOSE,
)


# Pattern for destructured header names
def _DESTRUCTURE_HEADER(delimiter: str) -> re.Pattern:
    assert len(delimiter) > 0

    return re.compile(
        rf"""
            ^ \s*               # Anchor to start of string, ignoring whitespace

            (?P<field>          # 'field' capture group
                .+?             # At least one (non-greedy) character
            )

            (?:
                # Delimiter, ignoring whitespace
                \s* {re.escape(delimiter)} \s*

                (?P<property>   # 'property' capture group
                    .*?         # Any number (non-greedy) of characters
                )

                (?:
                    \s* \( \s*  # Opening paren, ignoring whitespace
                    (?P<unit>   # 'unit' capture group
                        .+?     # At least one (non-greedy) character
                    )

                    \s* \)      # Closing paren, ignoring whitespace
                )?              # Unit is optional
            )?                  # Destructuring metadata is optional

            \s* $               # Anchor to end of string, ignoring whitespace
        """,
        re.VERBOSE,
    )


## MAPPING SCHEMA HANDLING #############################################


class CSVType(Enum):
    """Known CSV data types"""

    # RHS is the expected serialisation (lowercase)
    TABULAR = "tabular"
    KEY_VALUE = "key-value"


# Schema definition for validating and deserialising CSV data types
_CSV_TYPE_SCHEMA = schema.And(str, schema.Use(str.lower), schema.Use(CSVType))


@dataclasses.dataclass
class Name:
    """Name with optional aliases"""

    name: str
    aliases: list[str] = dataclasses.field(default_factory=list)

    @property
    def all_names(self) -> list[str]:
        """
        @return  List of all defined names
        """
        return [self.name] + self.aliases

    @staticmethod
    def _has_duplicates(*names: Name) -> bool:
        # Flatten input names and aliases into a list
        everything = [name for each in names for name in each.all_names]
        return len(everything) != len(set(everything))


# Schema definition for validating and deserialising names with
# optional aliases
_NAME_SCHEMA = schema.Or(
    schema.And(str, schema.Use(Name)),
    schema.And(
        {"name": str, schema.Optional("aliases"): [str]},
        schema.Use(lambda value: Name(**value)),
    ),
)


@dataclasses.dataclass
class Section:
    """Illumina Sample Sheet section mapping schema"""

    # NOTE The Section class has a dual use, which covers the
    # deserialisation of the schema YAML and also the internal
    # representation when cross-referenced with the input Illumina
    # Sample Sheet. To that end, the `section` attribute can represent
    # either the Illumina Sample Sheet section name/aliases _or_ the
    # parsed output key for that section.

    section: Name | str
    type: CSVType
    case_sensitive: bool = False
    delimiter: str = "|"
    required: list[Name] = dataclasses.field(default_factory=list)
    optional: list[Name] = dataclasses.field(default_factory=list)


# Schema definition for validating and deserialising sections
_SECTION_SCHEMA = schema.And(
    {
        "section": _NAME_SCHEMA,
        "type": _CSV_TYPE_SCHEMA,
        schema.Optional(
            schema.And(
                "case-sensitive",
                schema.Use(lambda _: "case_sensitive"),
            )
        ): bool,
        schema.Optional("delimiter"): str,
        schema.Optional("required"): [_NAME_SCHEMA],
        schema.Optional("optional"): [_NAME_SCHEMA],
    },
    schema.Use(lambda value: Section(**value)),
)

# Schema definition for validating our YAML mapping schema
_MAPPING_SCHEMA = schema.Schema(
    {str: _SECTION_SCHEMA},
    ignore_extra_keys=True,
)


class MappingSchema(Mapping[str, Section]):
    """
    Illumina Sample Sheet mapping schema
    (i.e., over multiple sections)
    """

    _sections: list[str]
    _schema: dict[str, Section]
    _normalised_lookup: dict[str, tuple[str, bool]]
    _output: dict[str, Name]

    def __init__(self, schema_io: TextIO, *, strict: bool = False):
        try:
            # Validate and deserialise the YAML
            parsed = _MAPPING_SCHEMA.validate(yaml.safe_load(schema_io))

            # Check for duplicate sections and keys; warn, or bail out
            # in strict mode, if any are found
            sections = []
            for key, section in parsed.items():
                sections.append(section.section)
                if Name._has_duplicates(*section.required, *section.optional):
                    logging.warning(
                        "Duplicate field name/aliases mappings detected in "
                        "section [%s] (output key '%s')",
                        section.section.name,
                        key,
                    )

                    if strict:
                        raise RuntimeError("Strictness violation in field mappings")

            if Name._has_duplicates(*sections):
                logging.warning("Duplicate section names/aliases detected")
                if strict:
                    raise RuntimeError("Strictness violation in section names")

        except (schema.SchemaError, RuntimeError) as exc:
            logging.debug("Schema validation failure: %s", exc)
            logging.critical("Illumina Sample Sheet mapping schema is invalid")
            sys.exit(1)

        # The YAML defines mappings for the output, which is more user
        # friendly. However, we need to invert this to see it from the
        # input section's perspective. That is, rather than output keys
        # mapping to input section names/aliases; we map input section
        # names/aliases to output keys. Because of the potential for
        # aliases, we duplicate sections in the internal representation,
        # to allow for efficient lookups.
        self._schema = {
            name: dataclasses.replace(section, section=output_key)
            for output_key, section in parsed.items()
            for name in section.section.all_names
        }

        # Lookup table for case-normalised section names/aliases
        self._normalised_lookup = {
            name.casefold(): (name, section.case_sensitive)
            for name, section in self._schema.items()
        }

        # We only expose the canonical section names in the Mapping
        # interface (aliases are "hidden")
        self._sections = [section.section.name for section in parsed.values()]

        # Keep a record of the expected output mapping (just the section
        # name), so we don't have to recompute it from the inversion
        # when we sanity check the output.
        self._output = {
            output_key: section.section for output_key, section in parsed.items()
        }

    def __getitem__(self, item: str) -> Section:
        normalised_item, case_sensitive = self._normalised_lookup[item.casefold()]
        return self._schema[item if case_sensitive else normalised_item]

    def __iter__(self) -> Iterator[str]:
        yield from self._sections

    def __len__(self) -> int:
        return len(self._sections)

    @property
    def output(self) -> dict[str, Name]:
        """
        @return  Output to section mapping
        """
        return self._output


## ILLUMINA SAMPLE SHEET HANDLING ######################################


class _ValueUnitT(TypedDict):
    value: str
    unit: str | None


def _create_value(value: str, unit: str | None = None) -> _ValueUnitT:
    return {"value": value.strip(), "unit": unit}


def _csv_row(*columns: str) -> str:
    # Serialise columns back into a CSV row
    output = StringIO()
    writer = csv.writer(output)
    writer.writerow(columns)

    return output.getvalue().strip()


# Transformed value type:
# * None
# * Scalar string
# * List of {value, unit}
# * Mapping of string to {value, unit}
_TransformedValueT = None | str | list[_ValueUnitT] | dict[str, _ValueUnitT]


class AbstractTransformer(ABC):
    """ABC for munging CSV row field(s) into an appropriate data structure"""

    @abstractmethod
    def transform(self, *columns: str):
        """
        Transform the incoming row into its respective representation
        """


class KVTransformer(AbstractTransformer):
    """Transformer for dictionaries of scalars or {value, unit}"""

    _property_column_map: dict[str, tuple[int, str | None]]

    def __init__(self):
        self._property_column_map = {}

    def __str__(self) -> str:
        pretty = ", ".join(
            f"{key}: {column_id} ({unit})" if unit else f"{key}: {column_id}"
            for key, (column_id, unit) in self._property_column_map.items()
        )
        return f"Key/value from columns {{{pretty}}}"

    def add_key(self, column_id: int, key: str, unit: str | None = None):
        """
        @param  column_id  Column index
        @param  key        Key
        @param  unit       Unit, if any
        """
        if key in self._property_column_map:
            logging.warning(
                "Section header destructuring contains duplicate key '%s'; "
                "the former will be overwritten",
                key,
            )

        self._property_column_map[key] = (column_id, unit)

    def transform(self, *columns: str) -> dict[str, _ValueUnitT]:
        return {
            key: _create_value(columns[column_id], unit)
            for key, (column_id, unit) in self._property_column_map.items()
        }


class ListTransformer(AbstractTransformer):
    """Transformer for lists of scalars or {value, unit}"""

    _column_ids: list[tuple[int, str | None]]

    def __init__(self):
        self._column_ids = []

    def __str__(self) -> str:
        pretty = ", ".join(
            f"{column_id} ({unit})" if unit else str(column_id)
            for column_id, unit in self._column_ids
        )
        return f"List from columns [{pretty}]"

    def add_element(self, column_id: int, unit: str | None = None):
        """
        @param  column_id  Column index
        @param  unit       Unit, if any
        """
        self._column_ids.append((column_id, unit))

    def transform(self, *columns: str) -> list[_ValueUnitT]:
        return [
            _create_value(columns[column_id], unit)
            for column_id, unit in self._column_ids
        ]


class ScalarTransformer(AbstractTransformer):
    """Transformer for scalar values"""

    _column_id: int

    def __init__(self, column_id: int):
        self._column_id = column_id

    def __str__(self) -> str:
        return f"Scalar from column {self._column_id}"

    def transform(self, *columns: str) -> str:
        return columns[self._column_id].strip()


class NullTransformer(AbstractTransformer):
    """Transformer for null values"""

    def __str__(self) -> str:
        return "Null"

    def transform(self, *_columns: str) -> None:
        return None


class AbstractMunger(ABC):
    """ABC for CSV munger implementations"""

    @abstractmethod
    def consume(self, *columns: str):
        """
        Consume a new row from the CSV source

        @param  columns  CSV row to consume
        """

    @abstractmethod
    def discharge(self):
        """
        Munge the consumed rows and emit an appropriate data structure

        @return  Munged data structure
        """


class TabularMunger(AbstractMunger):
    """CSV munger for headered, tabular CSV data"""

    _header: None | dict[str, AbstractTransformer]
    _expected_columns: None | int

    _case_sensitive: bool
    _destructure_header: re.Pattern

    # Map required and optional fields (and their aliases),
    # case-normalised, to their case-preserved, output names
    # (i.e., per the mapping schema)
    _required_lookup: dict[str, str]
    _optional_lookup: dict[str, str]

    _contents: list[dict[str, _TransformedValueT]]

    def __init__(
        self,
        case_sensitive: bool = False,
        delimiter: str = "|",
        required: list[Name] | None = None,
        optional: list[Name] | None = None,
    ):
        assert required is not None or optional is not None

        self._header = None
        self._expected_columns = None

        self._case_sensitive = case_sensitive
        self._destructure_header = _DESTRUCTURE_HEADER(delimiter)

        # NOTE We can compute the actual output fields from this
        # dictionary (same for _optional_lookup) with something like:
        #
        # ```python
        # [
        #     output_key
        #     for normalised_key, output_key in self._required_lookup.items()
        #     if normalised_key == self._normalise_case(output_key)
        # ]
        # ```
        self._required_lookup = {
            self._normalise_case(name): each.name
            for each in (required or [])
            for name in each.all_names
        }

        self._optional_lookup = {
            self._normalise_case(name): each.name
            for each in (optional or [])
            for name in each.all_names
        }

        self._contents = []

    def _normalise_case(self, field: str) -> str:
        # Case-normalise header field
        return field if self._case_sensitive else field.casefold()

    def _output_field(self, destructured_field: str, verbatim: str) -> tuple[str, bool]:
        # Map the destructured input field to a required or optional
        # output field, taking case-sensitivity and aliases into
        # account; if there's no match, then fallback to the verbatim
        # (un-destructured) field. Return a tuple of this output field
        # and whether it's specified in the mapping schema.
        normalised_field = self._normalise_case(destructured_field)

        try:
            return self._required_lookup[normalised_field], True

        except KeyError:
            try:
                return self._optional_lookup[normalised_field], True

            except KeyError:
                return verbatim, False

    def _parse_header(self, *columns: str):
        # Trim trailing empty columns: CSV writers (can) make the number
        # of columns consistent across the whole file. That is, the
        # number of columns in the input will be equal to the widest
        # section's data. This trimming is so we can sanity check each
        # section's data, regardless of the widest section. (It also
        # prevents re.match from failing, below.)

        # We assume every row has at least the same number of columns as
        # headers (it could have more, but those are presumed empty and
        # will be ignored if not).
        columns = trim_empty(*columns)
        self._expected_columns = len(columns)

        fields: dict[str, AbstractTransformer] = {}

        for column_id, header_match in enumerate(
            map(self._destructure_header.match, columns)
        ):
            # header_match is an re.Match object, with the following
            # groups: 'field', 'property' and 'unit', per
            # _DESTRUCTURE_HEADER. This should always be the case,
            # unless empty strings exist within the header.
            try:
                assert header_match is not None

                # Verbatim (i.e., non-destructured) input from CSV
                verbatim = header_match.string

                # Destructured field name and metadata
                input_field = header_match["field"]
                input_property: str | None = header_match["property"]
                input_unit: str | None = header_match["unit"]

                # Output field name and whether it is specified in the
                # mapping schema (i.e., subject to destructuring)
                output_field, is_specified = self._output_field(input_field, verbatim)

            except AssertionError:
                logging.critical("Unparsable header found in column %s", column_id)
                sys.exit(1)

            try:
                transformer: AbstractTransformer
                match (is_specified, input_property):
                    # Unspecified and metadata-less specified fields
                    case (False, _) | (True, None):
                        if output_field in fields:
                            logging.warning(
                                "Section in input contains duplicate header '%s' "
                                "(mapped to '%s' in the output); "
                                "the former will be overwritten",
                                input_field,
                                output_field,
                            )

                        fields[output_field] = ScalarTransformer(column_id)

                    # Destructuring only happens for specified fields

                    case (True, "[]"):
                        transformer = fields.setdefault(output_field, ListTransformer())
                        assert isinstance(transformer, ListTransformer), "list"

                        transformer.add_element(column_id, input_unit)

                    case (True, key):
                        transformer = fields.setdefault(output_field, KVTransformer())
                        assert isinstance(transformer, KVTransformer), "key-value"

                        # 'key' has to be a string here, but the type
                        # checker doesn't seem to be able to infer that
                        transformer.add_key(column_id, cast(str, key), input_unit)

            except AssertionError as exc:
                logging.critical(
                    "Mismatched property type in destructuring metadata: "
                    "Cannot add a %s transform to existing header '%s' "
                    "(mapped to '%s' in the output) of a different type",
                    exc,
                    input_field,
                    output_field,
                )
                sys.exit(1)

        # Decant the header transformers parsed from the input such that
        # required keys come first (failing if they don't exist),
        # followed by optional keys (replaced with nulls if they don't
        # exist), and then everything else
        self._header = {}

        all_required = True
        for normalised_key, output_key in self._required_lookup.items():
            # NOTE We only want the canonical headers
            if self._normalise_case(output_key) == normalised_key:
                try:
                    self._header[output_key] = fields.pop(output_key)

                except KeyError:
                    logging.critical(
                        "Required header '%s' (or its aliases) not found in input",
                        output_key,
                    )

                    all_required = False

        if not all_required:
            sys.exit(1)

        for normalised_key, output_key in self._optional_lookup.items():
            # NOTE We only want the canonical headers
            if self._normalise_case(output_key) == normalised_key:
                self._header[output_key] = fields.pop(output_key, NullTransformer())

        self._header |= fields

    def consume(self, *columns: str):
        if self._header is None:
            self._parse_header(*columns)
            logging.debug("Following header transformations found:")
            for header, transform in cast(dict, self._header).items():
                logging.debug("* %s: %s", header, transform)

        else:
            if len(columns) < cast(int, self._expected_columns):
                # NOTE This should never happen as the number of columns
                # in each row will be equal to the widest section in the
                # input (presuming it's exported as a compliant CSV)
                logging.critical("Row with insufficient columns found in section")
                logging.debug(
                    "Expected %s columns, found: %s",
                    self._expected_columns,
                    _csv_row(*columns),
                )
                sys.exit(1)

            self._contents.append(
                {
                    key: transformer.transform(*columns)
                    for key, transformer in self._header.items()
                }
            )

    def discharge(self) -> list[dict[str, _TransformedValueT]]:
        return self._contents


class KVMunger(AbstractMunger):
    """CSV munger for key/value CSV data"""

    # NOTE Key/value CSV data is just tabular data, with a single row,
    # that is transposed. We exploit this fact and let TabularMunger do
    # all the work.

    _tabular: TabularMunger
    _rows: list[tuple[str, str]]

    def __init__(
        self,
        case_sensitive: bool = False,
        delimiter: str = "|",
        required: list[Name] | None = None,
        optional: list[Name] | None = None,
    ):
        self._tabular = TabularMunger(case_sensitive, delimiter, required, optional)
        self._rows = []

    def consume(self, *columns: str):
        key, value, *_ = columns
        self._rows.append((key, value))

    def discharge(self) -> dict[str, _TransformedValueT]:
        for row in zip(*self._rows):
            self._tabular.consume(*row)

        return self._tabular.discharge()[0]


def is_empty(*columns: str) -> bool:
    """
    Are all columns empty?

    @param   columns  CSV columns
    @return  Whether all columns are empty, after stripping
    """
    return all(not column.strip() for column in columns)


def trim_empty(*columns: str) -> tuple[str, ...]:
    """
    Trim empty columns from end of tuple

    @param   columns  CSV columns
    @return  Input tuple with trailing empty strings stripped
    """
    length = len(columns)
    for i in range(length):
        if columns[length - i - 1].strip() != "":
            return columns[: length - i]

    return ()


def get_section(*columns: str) -> str | None:
    """
    @param   columns  CSV columns
    @return  Illumina Sample Sheet section name, if found in the first
             column, otherwise None
    """
    head, *tail = columns
    if (match := _SECTION_MARKER.search(head.strip())) and is_empty(*tail):
        return match.group(0).strip()

    return None


## HELPERS #############################################################


@dataclasses.dataclass
class _ConsumeState:
    key: str
    section: str
    munger: AbstractMunger

    def discharge_into(self, output: dict[str, object]):
        """
        Discharge the consumption state into the given output

        @param  output  Dictionary in which to discharge
        """
        logging.info(
            "Discharging [%s] section into '%s'",
            self.section,
            self.key,
        )

        output[self.key] = self.munger.discharge()


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
    description, epilog = re.split(r"^---$", __doc__, flags=re.MULTILINE)

    parser = argparse.ArgumentParser(
        description=description.strip(),
        epilog=epilog.strip(),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "schema",
        metavar="SCHEMA",
        type=argparse.FileType(mode="rt"),
        help="Illumina Sample Sheet mapping schema (YAML)",
    )

    parser.add_argument(
        "--sample-sheet",
        type=argparse.FileType(mode="rt"),
        default="-",
        help="input Illumina Sample Sheet (CSV; defaults to stdin)",
    )

    parser.add_argument(
        "--output",
        type=argparse.FileType(mode="xt"),
        default="-",
        help="output file (YAML; defaults to stdout)",
    )

    parser.add_argument(
        "--allow-ambiguous-schema",
        action="store_true",
        help="allow duplicates in the schema section and field mappings",
    )

    return parser.parse_args(args)


## ENTRYPOINT ##########################################################


def main(*argv: str) -> int:
    """Entrypoint"""
    _configure_logging()
    args = _parse_arguments(*argv)

    if args.allow_ambiguous_schema:
        logging.warning(
            "Ambiguous schema section and field mappings have been allowed. "
            "This can lead to unexpected behaviour!"
        )

    logging.info(
        "Reading Illumina Sample Sheet mapping schema from %s", args.schema.name
    )
    mapping = MappingSchema(args.schema, strict=not args.allow_ambiguous_schema)
    logging.info("Illumina Sample Sheet mapping schema loaded and validated")

    logging.info("Reading Illumina Sample Sheet from %s", args.sample_sheet.name)
    reader = csv.reader(args.sample_sheet)

    output: dict[str, object] = {}
    state: _ConsumeState | None = None

    row_id = 0
    for row in reader:
        row_id += 1

        # Handle newly found sections in the input
        if new_section := get_section(*row):
            # Discharge any existing consumption and reset the state machine
            if state:
                state.discharge_into(output)
                state = None

            logging.debug("Found section [%s]", new_section)

            # Mapped section
            if section := mapping.get(new_section):
                logging.info("Section [%s] maps to '%s'", new_section, section.section)

                munger: AbstractMunger
                match section.type:
                    case CSVType.TABULAR:
                        munger = TabularMunger(
                            section.case_sensitive,
                            section.delimiter,
                            section.required,
                            section.optional,
                        )

                    case CSVType.KEY_VALUE:
                        munger = KVMunger(
                            section.case_sensitive,
                            section.delimiter,
                            section.required,
                            section.optional,
                        )

                state = _ConsumeState(cast(str, section.section), new_section, munger)

        # Otherwise, consume row (when appropriate)
        elif state and not is_empty(*row):
            logging.debug("Consuming [%s] data from CSV row %s", state.section, row_id)
            state.munger.consume(*row)

    # Discharge any remaining consumption
    if state:
        state.discharge_into(output)

    # Check we have all the desired output
    if output.keys() != mapping.output.keys():
        logging.critical(
            "The following sections (or their aliases) were not found in the input:"
        )

        # It should always be the case that output.keys() is a subset of
        # mapping.output.keys(); i.e., we don't have to cater for extra,
        # unspecified output, just incomplete output.
        for key in mapping.output.keys() - output.keys():
            logging.critical(
                "* [%s], which maps to '%s' in the output",
                mapping.output[key].name,
                key,
            )

        sys.exit(1)

    logging.info("Writing parsed YAML to %s", args.output.name)
    yaml.safe_dump(output, args.output, explicit_start=True, sort_keys=False)
    args.output.flush()

    logging.info("Parsing complete")
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
