
import random
import requests
import json

import pandas as pd
from pandas import MultiIndex, CategoricalDtype
from pandas import IndexSlice as idx


FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"


def check_value(value, parent_key=None):
    # print(f"Checking value: '{value}'")
    if isinstance(value, dict):
        return walk_dict(value, parent_key)
    elif isinstance(value, list):
        value = check_list(value)
        if isinstance(value, list):
            return value
            # for subvalue in value:
            #     check_value(subvalue, parent_key)
        else:
            return check_value(value, parent_key)
    elif isinstance(value, (str, int, float, tuple, type(None))):
        return value


def check_list(input_list):
    # print(f"Checking list: '{input_list}'")
    list_len = len(input_list)
    if list_len == 0:
        return None
    if list_len == 1:
        return input_list[0]
    return input_list


def walk_dict(input_dict, parent_key=None):
    # print(f"Checking dict (parent_key='{parent_key}': '{input_dict}'")
    new_dict = dict()
    for key, value in input_dict.items():
        # print(f"Current key: '{key}'")
        if parent_key is None:
            new_key = key
        elif isinstance(parent_key, tuple):
            new_key = parent_key[0], f"{parent_key[-1]}.{key}"
        else:
            new_key = (parent_key, key)
        results = check_value(value, new_key)
        if isinstance(results, dict):
            new_dict |= results
        else:
            new_dict[new_key] = results
    new_dict = {x if isinstance(x, tuple) else ("file", x): y
                for x, y in new_dict.items()}
    # print(new_dict)
    return new_dict


def retrieve_file_information(sequencing_type):
    sequencing_type = sequencing_type.upper()

    file_fields = ",".join([
        "file_name",
        "file_size",
        "data_type",
        # "cases.case_id",
        "experimental_strategy",
        "cases.samples.sample_type",
        "cases.project.program.name",
        "md5sum"
        # "cases.diagnoses.tissue_or_organ_of_origin",
        # "cases.samples.tissue_type",
        ])

    file_expand = ",".join([
        "cases",
        "cases.project",
        "analysis.metadata.read_groups"
        # "cases.project.program",
        # "cases.samples"
        # "cases.tissue_source_site",
        # "associated_entities"
        ])

    file_filter = dict(
        op="and",
        content=[
            dict(op="=", content=dict(field="experimental_strategy", value=sequencing_type)),
            dict(op="=", content=dict(field="data_type", value="Aligned Reads")),
            dict(op="in", content=dict(field="cases.samples.sample_type", value="Primary Tumor")),
            dict(op="=", content=dict(field="cases.project.program.name", value="TCGA"))
            ]
        )

    file_params = dict(
        filters=json.dumps(file_filter),
        fields=file_fields,
        expand=file_expand,
        size="100000"
        )

    response = requests.post(FILES_ENDPOINT,
                             headers={"Content-Type": "application/json"},
                             json=file_params)
    if not response.ok:
        raise Exception

    file_table = json.loads(response.content.decode("utf-8"))["data"]["hits"]
    file_table = [walk_dict(entry) for entry in file_table]
    mdex = MultiIndex.from_tuples(file_table[0].keys())
    file_table = pd.DataFrame(file_table, columns=mdex)
    file_table[("file", "experimental_strategy")] = file_table[("file", "experimental_strategy")].astype(CategoricalDtype(["WGS", "WXS"]))
    file_table.set_index([("cases", "case_id"), ("file", "id"), ("file", "experimental_strategy")],
                         inplace=True, drop=False)
    file_table.index.rename(["case_id", "file_id", "sequencing"], inplace=True)
    file_table = file_table[file_table.columns.sort_values(ascending=False)]
    return file_table


def retrieve_all_data():
    wgs = retrieve_file_information("WGS")
    wxs = retrieve_file_information("WXS")
    data = pd.concat([wgs, wxs])
    data = data.sort_index()
    return data


def filter_for_paired_sequences(table: pd.DataFrame, true_pairs: bool = False):
    seq_set = {"WGS", "WXS"}
    bool_filter = [set(vals := [x[1] for x in table.loc[case].index]) == seq_set
                   and (not true_pairs or len(vals) == 2)
                   for (case, _, _), _ in table.iterrows()]
    return table.loc[bool_filter]


def filter_for_samples_with_known_capture_kits(table: pd.DataFrame,
                                               include_wgs: bool = True):
    has_kit = table.loc[~ pd.isna(
        table.analysis["metadata.read_groups.target_capture_kit_name"]
        )]
    if not include_wgs:
        return has_kit
    has_kit.index = has_kit.index.remove_unused_levels()
    wgs = table.loc[idx[:, :, ["WGS"]]]
    # wgs.index = wgs.index.remove_unused_levels()
    wgs = wgs.loc[[x in has_kit.index.levels[0] for x in wgs.reset_index()["case_id"]]]
    has_kit = pd.concat([has_kit, wgs]).sort_index()
    return has_kit


def sort_by_project_size(table):
    proj_names = ("cases", "project.project_id")
    counts = table[proj_names].value_counts()
    table = table.sort_values(by=[proj_names, ("cases", "case_id"), ("file", "experimental_strategy")],
                              key=lambda col: [counts[x] if col.name == proj_names
                                               else x for x in col],
                              ascending=False)
    return table


def randomly_select_cases_from_project(table: pd.DataFrame, project: str, k: int):
    subset = table.loc[table.cases["project.project_id"] == project]
    choices = random.sample(range(subset.shape[0]), k=k)
    choices = [subset.iloc[choice].cases.case_id for choice in choices]
    results = table.loc[choices]
    return results


def export_manifests(table, output_prefix, seq_type):
    # manifest = table.file
    cols = ["id", "file_name", "md5sum", "file_size"]
    cols = list(zip(["file"]*len(cols), cols))
    cols += [("cases", "case_id"),
             ("cases", "project.project_id"),
             ("analysis", "metadata.read_groups.target_capture_kit_name")]
    manifest = table.loc[idx[:, cols]]
    manifest.columns = [x[1].split(".")[-1] for x in manifest.columns]
    if seq_type == "paired":
        wxs_manifest = manifest.loc[idx[:, :, "WXS"]]
        wgs_manifest = manifest.loc[idx[:, :, "WGS"]]

        for seq, manifest in zip(["wxs", "wgs"], [wxs_manifest, wgs_manifest]):
            output = f"{output_prefix}.manifest.{seq}.tsv"
            manifest.to_csv(output, sep="\t", index=False)
    else:
        output = f"{output_prefix}.manifest.{seq_type}.tsv"
        manifest.to_csv(output, sep="\t", index=False)
