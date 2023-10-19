
import shutil
import subprocess

import pandas as pd


def read_manifest(manifest_file):
    manifest = pd.read_csv(manifest_file, sep="\t")
    manifest["downloaded"] = False
    return manifest


def download_gdc_file(file_id, token, download_directory):
    print(f"Downloading file: {file_id}...", end="")
    subprocess.run(["gdc-client", "download",
                    "--token-file", token,
                    "--dir", download_directory,
                    file_id,
                    "1>", f"{download_directory}/download.log",
                    "2>&1"],
                   check=True, capture_output=True, text=True)
    shutil.move(f"{download_directory}/{file_id}")
    print(f"\rDownloading file: {file_id} complete.")








