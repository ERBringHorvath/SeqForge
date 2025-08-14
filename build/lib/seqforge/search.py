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

from .shared.constants import FIELD_ALIASES

#replace whitespace with underscore
def normalize_key(key):
    return re.sub(r'\W+', '_', key.strip().lower())

def flatten_biosample_attributes(attributes):
    flat = {}
    for attr in attributes:
        name = normalize_key(attr.get("name", "unknown"))
        value = attr.get("value", "not_specified")
        flat[name] = value
    return flat

def parse_json(json_input, extract_all, fields, logger):
    json_files = []

    if os.path.isdir(json_input):
        json_files = [
            os.path.join(json_input, f)
            for f in os.listdir(json_input)
            if f.lower().endswith('.json')
        ]
        if not json_files:
            msg = f"No JSON filoes found in directory: {json_input}"
            print(f"\n\033[91m{msg}\033[0m")
            logger.warning(msg)
            return []
    else:
        json_files = [json_input]

    metadata = []
    for json_file in json_files:
        try:
            with open(json_file) as f:
                data = json.load(f)

            if 'reports' not in data:
                msg = f"File skipped (no 'reports' key): {json_file}"
                print(f"\n\033[91m{msg}\033[0m")
                logger.warning(msg)
                continue

            reports = data['reports']
            for report in reports:
                try:
                    biosample_attrs = report['assembly_info']['biosample']['attributes']
                    flat_attrs = flatten_biosample_attributes(biosample_attrs)

                    flat_attrs["top_accession"] = report.get("accession", "not_specified")
                    flat_attrs["top_organism"] = (
                        report.get("organism", {}).get("organism_name") or
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
                    msg = f"Keyerror in {json_file}: {e}"
                    print(f"\n\033[91m{msg}\033[0m")
                    logger.warning(msg)
                except Exception as e:
                    msg = f"Unexpected error in {json_file}: {e}"
                    print(f"\n\033[91m{msg}\033[0m")
                    logger.warning(msg)

        except Exception as e:
            msg = f"Failed to load JSON file {json_file}"
            print(f"\n\033[91m{msg}\033[0m")
            logger.warning(msg)

    return metadata

def parse_genbank(genbank_input, extract_all, fields, logger):
    metadata = []
    gb_extensions = ('.gb', '.gbk', '.genbank', '.gp')

    input_files = []
    if os.path.isdir(genbank_input):
        input_files = [
            os.path.join(genbank_input, f)
            for f in os.listdir(genbank_input)
            if f.lower().endswith(gb_extensions)
        ]
        if not input_files:
            msg = f"No GenBank files found in directory: {genbank_input}"
            print(f"\n\033[91m{msg}\033[0m")
            logger.warning(msg)
            return []
    else:
        input_files = [genbank_input]
    
    for gb_file in input_files:
        try:
            for record in SeqIO.parse(gb_file, "genbank"):
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
                            entry['field'] = 'not_specified'
                    
                metadata.append(entry)
        except Exception as e:
            msg = f"Failed to parse GenBank file '{gb_file}': {e}"
            print(f"\n\033[91m{msg}\033[0m")
            logger.warning(msg)
            continue
    
    return metadata

def run_search(args):
    logger = logging.getLogger("search")
    logger.setLevel(logging.INFO)
    logger.propagate = False

    if not args.all and not args.fields:
        print("\n\033[93mPlease use either '--all' to extract all available metadata "
              "or '--fields <fields>' to extract specific items\n"
              "Use 'seqforge search -h' to print a list of available fields\033[0m")
        return

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

    ext = None
    if os.path.isdir(args.input):
        ext = 'dir'
    else:
        ext = os.path.splitext(args.input)[-1].lower()

    extract_all = args.all
    fields = [normalize_key(f) for f in args.fields] if args.fields else []

    metadata = []

    if os.path.isdir(args.input):
        if args.json:
            metadata.extend(parse_json(args.input, extract_all, fields, logger))
        elif args.gb: 
            metadata.extend(parse_genbank(args.input, extract_all, fields, logger))
        else:
            metadata.extend(parse_json(args.input, extract_all, fields, logger))
            metadata.extend(parse_genbank(args.input, extract_all, fields, logger))
    else:
        if ext == ".json":
            metadata = parse_json(args.input, extract_all, fields, logger)
        elif ext in [".gb", ".gbk", ".genbank", ".gp"]:
            metadata = parse_genbank(args.input, extract_all, fields, logger)
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
