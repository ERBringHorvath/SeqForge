# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Elijah Bring Horvath

import json
import os
import pandas as pd
from Bio import SeqIO
from datetime import datetime
import re
import logging
from datetime import datetime


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


def parse_json(json_file, extract_all, fields, logger):
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
            msg = f"KeyError for report: {report.get('accession', 'unknown')} - {e}"
            print(f"\033[93m{msg}\033[0m")
            logger.warning(msg)
        except Exception as e:
            msg = f"Unexpected error: {e}"
            print(f"\033[91m{msg}\033[0m")
            logger.warning(f"Unexpected error: {e}")

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
    logger = logging.getLogger("search")
    logger.setLevel(logging.INFO)
    logger.propagate = False

    if not any(isinstance(h, logging.FileHandler) and h.baseFilename.endswith("seqforge_search.log")
            for h in logger.handlers):
        file_handler = logging.FileHandler("seqforge_search.log", mode='a')
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    start_time = datetime.now()
    print(f"Search started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Search started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

    if not os.path.exists(args.input):
        msg = f"Error: Input file not found at {args.input}"
        print(f"\033[91m{msg}\033[0m")
        logger.warning(msg)
        return

    ext = os.path.splitext(args.input)[-1].lower()
    extract_all = args.all
    fields = [normalize_key(f) for f in args.fields] if args.fields else []

    if ext == ".json":
        metadata = parse_json(args.input, extract_all, fields)
    elif ext in [".gb", ".gbk", ".genbank", ".gp"]:
        metadata = parse_genbank(args.input, extract_all, fields)
    else:
        msg = "Error: Unsupported file format. Must be .json or .gb/.gbk/.genbank/.gp (NCBI)"
        print(f"\033[91m{msg}\033[0m")
        logger.warning(msg)
        return
    
    if args.all and args.fields:
        print("\n\033[91mError: --all and --fields are incompatible arguments. Please choose one or the other\033[0m")
        return

    if not metadata:
        msg = "No metadata extracted"
        print(f"\033[93m{msg}\033[0m")
        logger.info(msg)
        return

    df = pd.DataFrame(metadata)

    out_path = args.output
    ext = os.path.splitext(out_path)[-1].lower()
    if ext == '.tsv':
        df.to_csv(out_path, sep='\t', index=False)
    elif ext == '.json':
        df.to_json(out_path, orient='records', indent=2)
    else:
        df.to_csv(out_path, index=False, sep=",", encoding="utf-8")

    print(f"\n\033[92mMetadata extracted and saved to {out_path}\033[0m\n")
    logger.info(f"Metadata extracted and saved to {out_path}")

    end_time = datetime.now()
    print(f"Search completed at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Search completed at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    msg = f"Total runtime: {str(end_time - start_time)}"
    logger.info(msg)
    print(msg)
