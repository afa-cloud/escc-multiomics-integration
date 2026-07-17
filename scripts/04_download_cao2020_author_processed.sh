#!/usr/bin/env bash

set -euo pipefail

DATA_ROOT="${RESEARCH_DATA_ROOT:-${HOME}/ResearchDataHub}"
INCOMING_ROOT="${DATA_ROOT}/_incoming"
ARIA2C="${ARIA2C:-aria2c}"

download_file() {
  local url="$1"
  local destination="$2"
  local expected_size="$3"
  local destination_dir
  local destination_name
  local actual_size

  destination_dir="$(dirname "${destination}")"
  destination_name="$(basename "${destination}")"
  mkdir -p "${destination_dir}"

  if command -v "${ARIA2C}" >/dev/null 2>&1; then
    "${ARIA2C}" \
      --allow-overwrite=true \
      --auto-file-renaming=false \
      --continue=true \
      --file-allocation=none \
      --max-connection-per-server=16 \
      --min-split-size=8M \
      --retry-wait=5 \
      --max-tries=10 \
      --split=16 \
      --dir="${destination_dir}" \
      --out="${destination_name}" \
      "${url}"
  else
    curl \
      --fail \
      --location \
      --retry 8 \
      --retry-all-errors \
      --continue-at - \
      --output "${destination}" \
      "${url}"
  fi

  actual_size="$(stat -f '%z' "${destination}")"
  if [[ "${expected_size}" != "unknown" && "${actual_size}" != "${expected_size}" ]]; then
    echo "文件大小不一致：${destination} expected=${expected_size} observed=${actual_size}" >&2
    exit 1
  fi
}

prepare_dataset_dir() {
  local root="$1"
  mkdir -p "${root}/00_source" "${root}/10_metadata" "${root}/90_manifests"
}

GSE149608_ROOT="${INCOMING_ROOT}/GEO/GSE149608/retrieved_20260711_author_processed"
GSE149609_ROOT="${INCOMING_ROOT}/GEO/GSE149609/retrieved_20260711_author_processed"
GSE151838_ROOT="${INCOMING_ROOT}/GEO/GSE151838/retrieved_20260711_author_processed"
NATURE_ROOT="${INCOMING_ROOT}/Nature_Communications/doi_10.1038_s41467-020-17227-z/supplement_snapshot_20231204"

prepare_dataset_dir "${GSE149608_ROOT}"
prepare_dataset_dir "${GSE149609_ROOT}"
prepare_dataset_dir "${GSE151838_ROOT}"
prepare_dataset_dir "${NATURE_ROOT}"

download_file \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149608/suppl/GSE149608_RAW.tar" \
  "${GSE149608_ROOT}/00_source/GSE149608_RAW.tar" \
  "2639984640"
download_file \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149608/soft/GSE149608_family.soft.gz" \
  "${GSE149608_ROOT}/10_metadata/GSE149608_family.soft.gz" \
  "4445"
download_file \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149608/suppl/filelist.txt" \
  "${GSE149608_ROOT}/10_metadata/filelist.txt" \
  "1546"

download_file \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149609/suppl/GSE149609_RAW.tar" \
  "${GSE149609_ROOT}/00_source/GSE149609_RAW.tar" \
  "22517760"
download_file \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149609/soft/GSE149609_family.soft.gz" \
  "${GSE149609_ROOT}/10_metadata/GSE149609_family.soft.gz" \
  "4108"
download_file \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149609/suppl/filelist.txt" \
  "${GSE149609_ROOT}/10_metadata/filelist.txt" \
  "2476"

download_file \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151838/suppl/GSE151838_RAW.tar" \
  "${GSE151838_ROOT}/00_source/GSE151838_RAW.tar" \
  "30720"
download_file \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151838/soft/GSE151838_family.soft.gz" \
  "${GSE151838_ROOT}/10_metadata/GSE151838_family.soft.gz" \
  "2219"
download_file \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151838/suppl/filelist.txt" \
  "${GSE151838_ROOT}/10_metadata/filelist.txt" \
  "230"

download_file \
  "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-17227-z/MediaObjects/41467_2020_17227_MOESM9_ESM.xlsx" \
  "${NATURE_ROOT}/00_source/41467_2020_17227_MOESM9_ESM.xlsx" \
  "6783218"

tar -tf "${GSE149608_ROOT}/00_source/GSE149608_RAW.tar" >/dev/null
tar -tf "${GSE149609_ROOT}/00_source/GSE149609_RAW.tar" >/dev/null
tar -tf "${GSE151838_ROOT}/00_source/GSE151838_RAW.tar" >/dev/null
gzip -t "${GSE149608_ROOT}/10_metadata/GSE149608_family.soft.gz"
gzip -t "${GSE149609_ROOT}/10_metadata/GSE149609_family.soft.gz"
gzip -t "${GSE151838_ROOT}/10_metadata/GSE151838_family.soft.gz"
unzip -tq "${NATURE_ROOT}/00_source/41467_2020_17227_MOESM9_ESM.xlsx" >/dev/null

shasum -a 256 \
  "${GSE149608_ROOT}/00_source/GSE149608_RAW.tar" \
  "${GSE149609_ROOT}/00_source/GSE149609_RAW.tar" \
  "${GSE151838_ROOT}/00_source/GSE151838_RAW.tar" \
  "${NATURE_ROOT}/00_source/41467_2020_17227_MOESM9_ESM.xlsx"
