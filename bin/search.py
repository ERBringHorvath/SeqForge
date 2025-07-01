# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Elijah Bring Horvath

import json
import os
import pandas as pd
from Bio import SeqIO
from datetime import datetime
import re


FIELD_ALIASES = {
    "accession": ["top_accession", "accession", "accessions"],
    "organism": ["top_organism", "organism", "organism_name"],
    "strain": ["strain"],
    "isolation_source": ["isolation_source", "source"],
    "host": ["host"],
    "region": ["geo_loc_name", "geo_location", "value"],
    "lat_lon": ["lat_lon"],
    "collection_date": ["collection_date"],
    "collected_by": ["collected_by"],
    "tax_id": ["top_tax_id", "tax_id"],
    "comment": ["comment"],
    "keywords": ["keywords"],
    "sequencing_tech": ["sequencing_tech", "sequencing_technology"],
    "release_date": ["top_release_date", "release_date", "date"]
}

def normalize_key(key):
    return re.sub(r'\W+', '_', key.strip().lower())


def flatten_biosample_attributes(attributes):
    flat = {}
    for attr in attributes:
        name = normalize_key(attr.get("name", "unknown"))
        value = attr.get("value", "not_specified")
        flat[name] = value
    return flat


def parse_json(json_file, extract_all, fields):
    with open(json_file) as f:
        data = json.load(f)

    if 'reports' not in data:
        print("\033[91mError: JSON format invalid. Expected 'reports' key.\033[0m")
        return []

    reports = data['reports']
    metadata = []
    for report in reports:
        try:
            biosample_attrs = report['assembly_info']['biosample']['attributes']
            flat_attrs = flatten_biosample_attributes(biosample_attrs)

            # Add top-level metadata with prefixed keys to prevent collisions
            flat_attrs["top_accession"] = report.get("accession", "not_specified")
            flat_attrs["top_organism"] = (
                report.get("organism", {}).get("organism_name") or
                report["assembly_info"]["biosample"].get("description", {}).get("organism", {}).get("organism_name") or
                "not_specified"
            )
            flat_attrs["top_tax_id"] = (
                report.get("organism", {}).get("tax_id") or
                report["assembly_info"]["biosample"].get("description", {}).get("organism", {}).get("tax_id") or
                "not_specified"
            )
            flat_attrs["top_release_date"] = report.get("annotation_info", {}).get("release_date", "not_specified")

            entry = {}

            if extract_all:
                entry.update(flat_attrs)
            else:
                for field in fields:
                    aliases = FIELD_ALIASES.get(field, [field])
                    found = False
                    for alias in aliases:
                        if alias in flat_attrs:
                            entry[field] = flat_attrs[alias]
                            found = True
                            break
                    if not found:
                        entry[field] = 'not_specified'

            metadata.append(entry)

        except KeyError as e:
            print(f"\033[93mKeyError for report: {report.get('accession', 'unknown')} - {e}\033[0m")
        except Exception as e:
            print(f"\033[91mUnexpected error: {e}\033[0m")

    return metadata

def parse_genbank(genbank_file, extract_all, fields):
    metadata = []
    for record in SeqIO.parse(genbank_file, "genbank"):
        entry = {'accession': record.annotations.get("accessions", [record.id])[0]}
        annotations = {normalize_key(k): v for k, v in record.annotations.items()}
        features = {normalize_key(k): v for k, v in record.features[0].qualifiers.items()} if record.features else {}

        if extract_all:
            entry.update(annotations)
            entry.update(features)
        else:
            for field in fields:
                aliases = FIELD_ALIASES.get(field, [field])
                found = False
                for alias in aliases:
                    if alias in annotations:
                        entry[field] = annotations[alias]
                        found = True
                        break
                    if alias in features:
                        entry[field] = features[alias][0] if isinstance(features[alias], list) else features[alias]
                        found = True
                        break
                if not found:
                    entry[field] = 'not_specified'

        metadata.append(entry)

    return metadata


def run_search(args):
    if not os.path.exists(args.input):
        print(f"\033[91mError: Input file not found at {args.input}\033[0m")
        return

    ext = os.path.splitext(args.input)[-1].lower()
    extract_all = args.all
    fields = [normalize_key(f) for f in args.fields] if args.fields else []

    if ext == ".json":
        metadata = parse_json(args.input, extract_all, fields)
    elif ext in [".gb", ".gbk", ".genbank"]:
        metadata = parse_genbank(args.input, extract_all, fields)
    else:
        print("\033[91mError: Unsupported file format. Must be .json or .gb/.gbk/.genbank\033[0m")
        return

    if not metadata:
        print("\033[93mNo metadata extracted.\033[0m")
        return

    df = pd.DataFrame(metadata)

    out_path = args.output
    ext = os.path.splitext(out_path)[-1].lower()
    if ext == '.tsv':
        df.to_csv(out_path, sep='\t', index=False)
    elif ext == '.json':
        df.to_json(out_path, orient='records', indent=2)
    else:
        df.to_csv(out_path, index=False)

    print(f"\n\033[92mMetadata extracted and saved to {out_path}\033[0m")
