#!/usr/bin/env bash

set -euo pipefail

DATA_ROOT="${RESEARCH_DATA_ROOT:-${HOME}/ResearchDataHub}"
ROOT="${DATA_ROOT}/_incoming/Metabolomics_Workbench/PR001876/retrieved_20260711"
ARIA2C="${ARIA2C:-aria2c}"

[[ -d "${DATA_ROOT}" ]] || {
  echo "ResearchDataHub 不存在或未挂载：${DATA_ROOT}" >&2
  exit 1
}
[[ -r "${DATA_ROOT}/CATALOG.tsv" ]] || {
  echo "CATALOG.tsv 不可读：${DATA_ROOT}/CATALOG.tsv" >&2
  exit 1
}

if awk -F '\t' '
  NR > 1 && $13 == "verified" &&
  ($5 ~ /(^|;)PR001876(;|$)/ || $5 ~ /ST003013|ST003014|ST003015|ST003025|ST003027/) {
    found = 1
  }
  END { exit(found ? 0 : 1) }
' "${DATA_ROOT}/CATALOG.tsv"; then
  echo "CATALOG.tsv 已存在 PR001876 或子研究的 verified 记录；请直接复用并核对版本。"
  exit 0
fi

mkdir -p \
  "${ROOT}/00_source" \
  "${ROOT}/10_metadata" \
  "${ROOT}/20_reusable" \
  "${ROOT}/90_manifests"

verify_md5() {
  local expected="$1"
  local file="$2"
  local observed

  observed="$(md5 -q "${file}")"
  if [[ "${observed}" != "${expected}" ]]; then
    echo "MD5 不一致：${file} expected=${expected} observed=${observed}" >&2
    return 1
  fi
}

size_matches() {
  local expected="$1"
  local file="$2"
  local observed

  observed="$(stat -f '%z' "${file}")"
  if [[ "${expected}" == minimum:* ]]; then
    [[ "${observed}" -ge "${expected#minimum:}" ]]
  elif [[ "${expected}" == "nonempty" ]]; then
    [[ "${observed}" -gt 0 ]]
  else
    [[ "${observed}" == "${expected}" ]]
  fi
}

download_file() {
  local url="$1"
  local destination="$2"
  local expected_size="$3"
  local expected_md5="${4:-}"
  local destination_dir
  local destination_name
  local partial
  local partial_name
  local actual_size

  destination_dir="$(dirname "${destination}")"
  destination_name="$(basename "${destination}")"
  partial="${destination}.part"
  partial_name="$(basename "${partial}")"
  mkdir -p "${destination_dir}"

  if [[ -f "${destination}" ]]; then
    actual_size="$(stat -f '%z' "${destination}")"
    if size_matches "${expected_size}" "${destination}" && \
      { [[ -z "${expected_md5}" ]] || verify_md5 "${expected_md5}" "${destination}"; }; then
      echo "已通过已有文件校验，跳过下载：${destination}"
      return 0
    fi

    mv "${destination}" "${destination}.invalid.$(date '+%Y%m%d%H%M%S')"
  fi

  if [[ -f "${partial}" ]] && \
    size_matches "${expected_size}" "${partial}" && \
    { [[ -z "${expected_md5}" ]] || verify_md5 "${expected_md5}" "${partial}"; }; then
    mv "${partial}" "${destination}"
    echo "已通过断点文件校验：${destination}"
    return 0
  fi

  if command -v "${ARIA2C}" >/dev/null 2>&1; then
    "${ARIA2C}" \
      --allow-overwrite=false \
      --auto-file-renaming=false \
      --continue=true \
      --file-allocation=none \
      --max-connection-per-server=8 \
      --min-split-size=4M \
      --retry-wait=5 \
      --max-tries=10 \
      --split=8 \
      --dir="${destination_dir}" \
      --out="${partial_name}" \
      "${url}"
  else
    curl --fail --location --retry 8 --retry-all-errors --continue-at - \
      --output "${partial}" "${url}"
  fi

  actual_size="$(stat -f '%z' "${partial}")"
  if ! size_matches "${expected_size}" "${partial}"; then
    echo "文件大小不一致：${partial} expected=${expected_size} observed=${actual_size}" >&2
    exit 1
  fi

  if [[ -n "${expected_md5}" ]]; then
    verify_md5 "${expected_md5}" "${partial}"
  fi

  mv "${partial}" "${destination}"
}

download_file \
  "https://www.metabolomicsworkbench.org/studydownload/ST003013_Rawfiles.zip" \
  "${ROOT}/00_source/ST003013_Rawfiles.zip" \
  "111464557" \
  "044e1a26d894560a825c616773e8f580"
download_file \
  "https://www.metabolomicsworkbench.org/studydownload/ST003014_Rawfiles.zip" \
  "${ROOT}/00_source/ST003014_Rawfiles.zip" \
  "88718718" \
  "cf0a3c8ba0a5fa6dfcfaca877f6da5cf"
download_file \
  "https://www.metabolomicsworkbench.org/studydownload/ST003015_Rawfiles.zip" \
  "${ROOT}/00_source/ST003015_Rawfiles.zip" \
  "58144183" \
  "73c2626b9c28385c008d7dc9f3195d8f"
download_file \
  "https://www.metabolomicsworkbench.org/studydownload/ST003025_Rawfiles.zip" \
  "${ROOT}/00_source/ST003025_Rawfiles.zip" \
  "637059" \
  "d444ca5d7b5dcf72a37b2f232596aa4f"
download_file \
  "https://www.metabolomicsworkbench.org/studydownload/ST003027_Rawfiles.zip" \
  "${ROOT}/00_source/ST003027_Rawfiles.zip" \
  "46783597" \
  "6b11b4d20538ecff1b808c41785b8b24"

download_file \
  "https://www.metabolomicsworkbench.org/data/study_textformat_view.php?ANALYSIS_ID=AN004946&MODE=d&STUDY_ID=ST003013" \
  "${ROOT}/10_metadata/ST003013_AN004946_mwtab.txt" \
  "minimum:10000"
download_file \
  "https://www.metabolomicsworkbench.org/data/study_textformat_view.php?ANALYSIS_ID=AN004947&MODE=d&STUDY_ID=ST003014" \
  "${ROOT}/10_metadata/ST003014_AN004947_mwtab.txt" \
  "minimum:10000"
download_file \
  "https://www.metabolomicsworkbench.org/data/study_textformat_view.php?ANALYSIS_ID=AN004948&MODE=d&STUDY_ID=ST003015" \
  "${ROOT}/10_metadata/ST003015_AN004948_mwtab.txt" \
  "minimum:10000"
download_file \
  "https://www.metabolomicsworkbench.org/data/study_textformat_view.php?ANALYSIS_ID=AN004960&MODE=d&STUDY_ID=ST003025" \
  "${ROOT}/10_metadata/ST003025_AN004960_mwtab.txt" \
  "minimum:10000"
download_file \
  "https://www.metabolomicsworkbench.org/data/study_textformat_view.php?ANALYSIS_ID=AN004962&MODE=d&STUDY_ID=ST003027" \
  "${ROOT}/10_metadata/ST003027_AN004962_mwtab.txt" \
  "minimum:10000"
download_file \
  "https://www.metabolomicsworkbench.org/data/study_textformat_view.php?ANALYSIS_ID=AN004963&MODE=d&STUDY_ID=ST003027" \
  "${ROOT}/10_metadata/ST003027_AN004963_mwtab.txt" \
  "minimum:10000"

verify_mwtab() {
  local file="$1"
  local study_id="$2"
  local analysis_id="$3"

  grep -q '^#METABOLOMICS WORKBENCH ' "${file}"
  grep -q 'PROJECT_ID:PR001876' "${file}"
  grep -q "STUDY_ID:${study_id}" "${file}"
  grep -q "ANALYSIS_ID:${analysis_id}" "${file}"
}

verify_mwtab "${ROOT}/10_metadata/ST003013_AN004946_mwtab.txt" "ST003013" "AN004946"
verify_mwtab "${ROOT}/10_metadata/ST003014_AN004947_mwtab.txt" "ST003014" "AN004947"
verify_mwtab "${ROOT}/10_metadata/ST003015_AN004948_mwtab.txt" "ST003015" "AN004948"
verify_mwtab "${ROOT}/10_metadata/ST003025_AN004960_mwtab.txt" "ST003025" "AN004960"
verify_mwtab "${ROOT}/10_metadata/ST003027_AN004962_mwtab.txt" "ST003027" "AN004962"
verify_mwtab "${ROOT}/10_metadata/ST003027_AN004963_mwtab.txt" "ST003027" "AN004963"

for archive in "${ROOT}"/00_source/*.zip; do
  unzip -tq "${archive}" >/dev/null
done

(
  cd "${ROOT}"
  shasum -a 256 00_source/*.zip 10_metadata/*.txt > 90_manifests/SHA256SUMS
  shasum -a 256 -c 90_manifests/SHA256SUMS
)
