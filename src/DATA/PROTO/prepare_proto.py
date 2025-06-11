import json
import math
from collections import defaultdict
from pathlib import Path
import re


def gcd(numbers):
    if not numbers:
        return 0
    result = numbers[0]
    for i in range(1, len(numbers)):
        result = math.gcd(result, numbers[i])
    return result
def stoichiometry_key(label):
    stoichiometry = label.split("_", maxsplit=1)[0]
    stoichiometry = re.split(r'[A-Z]', stoichiometry)[1:]
    key = []
    for s in stoichiometry:
        if s == '':
            key.append(1)
        else:
            key.append(int(s))
    gcd_value = gcd(key)
    if gcd_value != 1:
        key = [k/gcd_value for k in key]
    return ":".join([str(k) for k in sorted(key, reverse=True)])

base_folder = Path(__file__).parent
collection_folder = base_folder / "aflow_prototype_encyclopedia" / "data"
legacy_file = base_folder / "legacy_data.json"
legacy_data = json.loads(legacy_file.read_text())
expected_keys = ["uid",
                 'label',
                 'alias',
                 'title',
                 'icsd',
                 'ccdc',
                 'prototype',
                 'prototype_tex',
                 'number_of_species',
                 'pearson_symbol',
                 'space_group_number',
                 'strukturbericht',
                 'strukturbericht_tex',
                 'part',
                 'parameter_symbols',
                 'parameter_values',
                 'compounds',
                 'legacy',
                 'notes',
                 'comments']

def prepare():
    lookup = dict(icsd={}, alias={}, label={}, ccdc={},
                  prototype=defaultdict(list), pearson_symbol=defaultdict(list),
                  space_group_number=defaultdict(list), part=defaultdict(list),
                  number_of_species=defaultdict(list), stoichiometry=defaultdict(list),
                  all=[], legacy=[], encyclopedia=[])
    combined = dict()
    for lib_entry_path in collection_folder.glob("*/info.json"):
        lib_entry = json.loads(lib_entry_path.read_text())
        for key in expected_keys:
            if not key in lib_entry:
                lib_entry[key] = None
        combined[lib_entry["uid"]] = lib_entry
        if lib_entry["label"] is None:
            print(f"problem: {key} {lib_entry[key]}")
            exit()
        lookup["label"][lib_entry["label"]] = lib_entry["uid"]
        base_label, num_id = lib_entry["label"].rsplit("-", maxsplit=1)
        for key in ["icsd", "ccdc"]:
            try:
                if lib_entry[key] is not None:
                    if "," in lib_entry[key]:
                        for se in lib_entry[key].split(","):
                            lookup[key][f"{key.upper()}_{se.strip()}"] = lib_entry["uid"]
                    else:
                        lookup[key][f"{key.upper()}_{lib_entry[key]}"] = lib_entry["uid"]

            except KeyError:
                print(lib_entry)
                raise

        if lib_entry["alias"] is not None:
            if isinstance(lib_entry["alias"], list):
                for entry in lib_entry["alias"]:
                    lookup["alias"][entry] = lib_entry["uid"]
            else:
                lookup["alias"][lib_entry["alias"]] = lib_entry["uid"]

        for key in ["prototype", "pearson_symbol", "space_group_number", "part", "number_of_species"]:
            if lib_entry[key] is not None:
                lookup[key][lib_entry[key]].append(lib_entry["uid"])

        stoichiometry = stoichiometry_key(lib_entry["label"])
        lookup["stoichiometry"][stoichiometry].append(lib_entry["uid"])
        for key in ["all", "encyclopedia"]:
            lookup[key].append(lib_entry["uid"])

        if int(num_id) == 1:
            lookup["alias"][base_label] = lib_entry["uid"]

    for uid, lib_entry in legacy_data.items():
        found=None
        for key in ["icsd", 'ccdc', 'label']:
            if lib_entry[key] in lookup[key]:
                found=lookup[key][lib_entry[key]]
                break
        if lib_entry["label"] in lookup["alias"]:
            found=lookup["alias"][lib_entry["label"]]

        if found is not None:
            if lib_entry["alias"] is not None:
                if combined[found]["alias"] is not None:
                    combined[found]["alias"] += lib_entry["alias"]
                    combined[found]["alias"] = list(set(combined[found]["alias"]))
                else:
                    combined[found]["alias"] = list(set(lib_entry["alias"]))
                for a in lib_entry["alias"]:
                    lookup["alias"][a] = found
            if combined[found]["prototype"] is None and lib_entry["prototype"] is not None:
                combined[found]["prototype"] = lib_entry["prototype"]
                lookup["prototype"][lib_entry["prototype"]].append(found)
            continue
        combined[uid] = lib_entry
        lookup["label"][lib_entry["label"]] = lib_entry["uid"]
        for key in ["icsd", "ccdc"]:
            try:
                if lib_entry[key] is not None:
                    if "," in lib_entry[key]:
                        for se in lib_entry[key].split(","):
                            lookup[key][f"{key.upper()}_{se.strip()}"] = lib_entry["uid"]
                    else:
                        lookup[key][f"{key.upper()}_{lib_entry[key]}"] = lib_entry["uid"]

            except KeyError:
                print(lib_entry)
                raise
        if lib_entry["alias"] is not None:
            if isinstance(lib_entry["alias"], list):
                for entry in lib_entry["alias"]:
                    lookup["alias"][entry] = lib_entry["uid"]
            else:
                lookup["alias"][lib_entry["alias"]] = lib_entry["uid"]

        for key in ["prototype", "pearson_symbol", "space_group_number", "number_of_species"]:
            if lib_entry[key] is not None:
                lookup[key][lib_entry[key]].append(lib_entry["uid"])

        stoichiometry = stoichiometry_key(lib_entry["label"])
        lookup["stoichiometry"][stoichiometry].append(lib_entry["uid"])
        lookup["legacy"].append(lib_entry["uid"])

    (base_folder / "lookup.json").write_text(json.dumps(lookup))
    (base_folder / "data.json").write_text(json.dumps(combined, indent=4))



if __name__ == '__main__':
    prepare()