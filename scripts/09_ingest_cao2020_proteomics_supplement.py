#!/usr/bin/env python3

import csv
import hashlib
import math
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile
import xml.etree.ElementTree as ET
import zipfile


PROJECT_ROOT = Path.cwd().resolve()
DATA_ROOT = Path(
    os.environ.get("RESEARCH_DATA_ROOT", str(Path.home() / "ResearchDataHub"))
)
STAGE_ROOT = (
    DATA_ROOT
    / "_incoming"
    / "Nature_Communications"
    / "doi_10.1038_s41467-020-17227-z"
    / "supplement_snapshot_20231204"
)
CANONICAL_ROOT = (
    DATA_ROOT
    / "datasets"
    / "public"
    / "Nature_Communications"
    / "doi_10.1038_s41467-020-17227-z"
    / "supplement_data6_snapshot_20231204"
)
SOURCE_NAME = "41467_2020_17227_MOESM9_ESM.xlsx"
SOURCE_PATH = STAGE_ROOT / "00_source" / SOURCE_NAME
DATASET_KEY = (
    "NATURE_COMMUNICATIONS_10.1038-s41467-020-17227-z_"
    "supplement_data6_snapshot_20231204"
)
DOWNLOAD_DATE = "2026-07-11"
REPAIR_EXISTING = "--repair-existing" in sys.argv
SOURCE_URL = (
    "https://static-content.springer.com/esm/art%3A10.1038%2F"
    "s41467-020-17227-z/MediaObjects/41467_2020_17227_MOESM9_ESM.xlsx"
)

NS_MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
NS_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
NS_PKG_REL = "http://schemas.openxmlformats.org/package/2006/relationships"

SHEET_PATIENTS = {
    "T_N_Proteomic_Dataset 1": {"M": "P1", "P": "P2", "S": "P6", "V": "P7", "Y": "Z1"},
    "T_N_Proteomic_Dataset 2": {"M": "P8", "P": "P9", "S": "P11", "V": "P7", "Y": "Z2"},
    "T_N_Proteomic_Dataset 3": {"M": "P12", "P": "P13", "S": "P15", "V": "P7", "Y": "Z3"},
}

RATIO_TRIPLETS = {
    "M": ("M", "N", "O"),
    "P": ("P", "Q", "R"),
    "S": ("S", "T", "U"),
    "V": ("V", "W", "X"),
    "Y": ("Y", "Z", "AA"),
}

BASE_COLUMNS = {
    "A": "group_id",
    "B": "hit_number",
    "C": "description",
    "D": "score",
    "E": "mass",
    "F": "coverage",
    "G": "sequence",
    "H": "same_sets",
    "I": "spectrum_count",
    "J": "unique_spectrum_count",
    "K": "peptide_count",
    "L": "unique_peptide_count",
    "AB": "accession",
    "AI": "uniprot_swissprot_accession",
    "AJ": "uniprot_swissprot_description",
    "AM": "uniprot_trembl_accession",
    "AO": "protein_or_domain",
    "AV": "kegg_annotation",
    "AW": "go_biological_process",
    "AX": "go_molecular_function",
    "AY": "go_cellular_component",
}


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def column_letters(reference):
    match = re.match(r"([A-Z]+)", reference)
    return match.group(1) if match else reference


def load_shared_strings(archive):
    root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    return [
        "".join(node.text or "" for node in item.iter(f"{{{NS_MAIN}}}t"))
        for item in root.findall(f"{{{NS_MAIN}}}si")
    ]


def workbook_sheets(archive):
    workbook = ET.fromstring(archive.read("xl/workbook.xml"))
    relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    targets = {
        rel.attrib["Id"]: rel.attrib["Target"]
        for rel in relationships.findall(f"{{{NS_PKG_REL}}}Relationship")
    }
    sheets = []
    for sheet in workbook.find(f"{{{NS_MAIN}}}sheets"):
        rel_id = sheet.attrib[f"{{{NS_REL}}}id"]
        target = targets[rel_id]
        member = target.lstrip("/") if target.startswith("/") else f"xl/{target}"
        sheets.append(
            {
                "name": sheet.attrib["name"],
                "state": sheet.attrib.get("state", "visible"),
                "member": member,
            }
        )
    return sheets


def cell_value(cell, strings):
    cell_type = cell.attrib.get("t", "")
    value_node = cell.find(f"{{{NS_MAIN}}}v")
    if cell_type == "inlineStr":
        return "".join(
            node.text or "" for node in cell.iter(f"{{{NS_MAIN}}}t")
        )
    if value_node is None or value_node.text is None:
        return ""
    value = value_node.text
    if cell_type == "s":
        return strings[int(value)]
    return value


def parse_sheet(archive, sheet, strings):
    root = ET.fromstring(archive.read(sheet["member"]))
    dimension_node = root.find(f"{{{NS_MAIN}}}dimension")
    dimension = dimension_node.attrib.get("ref", "") if dimension_node is not None else ""
    formula_count = len(root.findall(f".//{{{NS_MAIN}}}f"))
    merged = root.find(f"{{{NS_MAIN}}}mergeCells")
    merged_count = int(merged.attrib.get("count", "0")) if merged is not None else 0
    table_parts = root.find(f"{{{NS_MAIN}}}tableParts")
    table_count = int(table_parts.attrib.get("count", "0")) if table_parts is not None else 0

    dimension_match = re.match(r"A1:([A-Z]+)([0-9]+)", dimension)
    if not dimension_match:
        raise RuntimeError(f"无法解析工作表范围：{sheet['name']} {dimension}")
    last_column, last_row = dimension_match.groups()

    rows = root.findall(f".//{{{NS_MAIN}}}row")
    headers = {}
    protein_rows = []
    for row in rows:
        row_number = int(row.attrib.get("r", "0"))
        values = {
            column_letters(cell.attrib.get("r", "")): cell_value(cell, strings)
            for cell in row.findall(f"{{{NS_MAIN}}}c")
        }
        if row_number == 1:
            headers = values
        elif row_number > 1:
            protein_rows.append((row_number, values))

    inventory = {
        "sheet_name": sheet["name"],
        "sheet_state": sheet["state"],
        "dimension": dimension,
        "row_count_including_header": last_row,
        "column_count": 51 if last_column == "AY" else "",
        "formula_count": formula_count,
        "merged_cell_count": merged_count,
        "table_count": table_count,
    }
    return inventory, headers, protein_rows


def parse_channel_header(header):
    match = re.match(r"^(\d+)N_(\d+)-VS-(\d+)T_(\d+)$", header)
    if not match:
        return "", "", "", ""
    patient_n, normal_channel, patient_t, tumor_channel = match.groups()
    return patient_n, normal_channel, patient_t, tumor_channel


def write_tsv(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)


def write_tsv_atomic(path, fieldnames, rows):
    path = Path(path)
    with tempfile.NamedTemporaryFile(
        mode="w",
        newline="",
        encoding="utf-8",
        dir=path.parent,
        prefix=f"{path.name}.",
        delete=False,
    ) as handle:
        temp_path = Path(handle.name)
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temp_path, path)


def main():
    if not SOURCE_PATH.is_file():
        raise RuntimeError(f"缺少来源工作簿：{SOURCE_PATH}")
    if CANONICAL_ROOT.exists() and not REPAIR_EXISTING:
        raise RuntimeError(f"规范目录已存在，拒绝覆盖：{CANONICAL_ROOT}")

    metadata_dir = STAGE_ROOT / "10_metadata"
    reusable_dir = STAGE_ROOT / "20_reusable"
    manifest_dir = STAGE_ROOT / "90_manifests"
    for directory in (metadata_dir, reusable_dir, manifest_dir):
        directory.mkdir(parents=True, exist_ok=True)

    inventory_rows = []
    header_rows = []
    mapping_rows = []
    patient_map_rows = []
    anomaly_rows = []
    protein_annotation_rows = []
    ratio_rows = []

    with zipfile.ZipFile(SOURCE_PATH) as archive:
        bad_member = archive.testzip()
        if bad_member is not None:
            raise RuntimeError(f"xlsx ZIP 成员损坏：{bad_member}")
        strings = load_shared_strings(archive)
        sheets = workbook_sheets(archive)
        if [sheet["name"] for sheet in sheets] != list(SHEET_PATIENTS):
            raise RuntimeError("工作表名称或顺序与预期不一致")

        for batch_index, sheet in enumerate(sheets, start=1):
            inventory, headers, protein_rows = parse_sheet(archive, sheet, strings)
            inventory["artifact_tool_verified"] = "yes"
            for column, header in headers.items():
                header_rows.append(
                    {
                        "sheet_name": sheet["name"],
                        "column": column,
                        "header": header,
                    }
                )

            for ratio_column, patient_id in SHEET_PATIENTS[sheet["name"]].items():
                ratio_col, quant_col, sig_col = RATIO_TRIPLETS[ratio_column]
                header = headers.get(ratio_col, "")
                patient_n, normal_channel, patient_t, tumor_channel = parse_channel_header(header)
                is_z = patient_id.startswith("Z")
                if is_z:
                    mapping_status = "excluded_unknown_z_ratio"
                    mapping_confidence = "low"
                    analysis_eligibility = "excluded_pending_identity"
                    note = "Z1/Z2/Z3 不是论文患者 ID；保留审计但不进入患者级矩阵"
                else:
                    expected_number = patient_id.removeprefix("P")
                    if patient_n != expected_number or patient_t != expected_number:
                        raise RuntimeError(
                            f"{sheet['name']} {ratio_column} 的患者标题与映射不一致：{header}"
                        )
                    mapping_status = "patient_ratio_mapped"
                    mapping_confidence = "high"
                    analysis_eligibility = "patient_ratio_available"
                    note = (
                        "P7 为三个 iTRAQ 批次共同技术参考，三次比值保留但不得计作三位患者"
                        if patient_id == "P7"
                        else "工作簿原始标题直接给出患者和 N/T 通道"
                    )

                mapping_row = {
                    "batch": batch_index,
                    "sheet_name": sheet["name"],
                    "ratio_column": ratio_col,
                    "quant_number_column": quant_col,
                    "significance_column": sig_col,
                    "original_header": header,
                    "paper_patient_id": "" if is_z else patient_id,
                    "source_label": patient_id,
                    "normal_channel": normal_channel,
                    "tumor_channel": tumor_channel,
                    "mapping_status": mapping_status,
                    "mapping_confidence": mapping_confidence,
                    "analysis_eligibility": analysis_eligibility,
                    "mapping_evidence": "Supplementary Data 6 workbook header + Cao 2020 patient list",
                    "notes": note,
                }
                mapping_rows.append(mapping_row)
                if not is_z:
                    patient_map_rows.append(mapping_row)

            sheet_anomaly_count = 0
            for row_number, values in protein_rows:
                group_id = values.get("A", "")
                if not re.fullmatch(r"[0-9]+(?:\.0+)?", group_id):
                    sheet_anomaly_count += 1
                    anomaly_rows.append(
                        {
                            "batch": batch_index,
                            "sheet_name": sheet["name"],
                            "source_row": row_number,
                            "populated_cell_count": len(values),
                            "populated_columns": ";".join(values),
                            "raw_a_preview": group_id[:500],
                            "exclusion_reason": "row does not begin with numeric GroupID; source cells are shifted or continuation text",
                            "analysis_eligibility": "excluded_structural_anomaly",
                        }
                    )
                    continue
                base = {
                    "batch": batch_index,
                    "sheet_name": sheet["name"],
                    "source_row": row_number,
                }
                for column, field in BASE_COLUMNS.items():
                    base[field] = values.get(column, "")
                protein_annotation_rows.append(base)

                for ratio_column in ("M", "P", "S", "V"):
                    patient_id = SHEET_PATIENTS[sheet["name"]][ratio_column]
                    ratio_col, quant_col, sig_col = RATIO_TRIPLETS[ratio_column]
                    ratio_header = headers.get(ratio_col, "")
                    _, normal_channel, _, tumor_channel = parse_channel_header(ratio_header)
                    source_ratio = values.get(ratio_col, "")
                    log2_tumor_vs_normal = ""
                    if source_ratio not in ("", "NA"):
                        try:
                            numeric_ratio = float(source_ratio)
                            if numeric_ratio > 0:
                                log2_tumor_vs_normal = f"{-math.log2(numeric_ratio):.12g}"
                        except ValueError:
                            pass
                    ratio_rows.append(
                        {
                            "batch": batch_index,
                            "sheet_name": sheet["name"],
                            "source_row": row_number,
                            "paper_patient_id": patient_id,
                            "ratio_header": ratio_header,
                            "normal_channel": normal_channel,
                            "tumor_channel": tumor_channel,
                            "normal_vs_tumor_ratio": source_ratio,
                            "log2_tumor_vs_normal": log2_tumor_vs_normal,
                            "quant_number": values.get(quant_col, ""),
                            "source_significance": values.get(sig_col, ""),
                            "cross_batch_reference": "yes" if patient_id == "P7" else "no",
                        }
                    )

            inventory["structural_anomaly_rows"] = sheet_anomaly_count
            inventory_rows.append(inventory)

    inventory_fields = [
        "sheet_name",
        "sheet_state",
        "dimension",
        "row_count_including_header",
        "column_count",
        "formula_count",
        "merged_cell_count",
        "table_count",
        "artifact_tool_verified",
        "structural_anomaly_rows",
    ]
    write_tsv(metadata_dir / "workbook_inventory.tsv", inventory_fields, inventory_rows)
    write_tsv(metadata_dir / "workbook_headers.tsv", ["sheet_name", "column", "header"], header_rows)
    mapping_fields = [
        "batch",
        "sheet_name",
        "ratio_column",
        "quant_number_column",
        "significance_column",
        "original_header",
        "paper_patient_id",
        "source_label",
        "normal_channel",
        "tumor_channel",
        "mapping_status",
        "mapping_confidence",
        "analysis_eligibility",
        "mapping_evidence",
        "notes",
    ]
    write_tsv(metadata_dir / "proteomics_mapping_audit.tsv", mapping_fields, mapping_rows)
    anomaly_fields = [
        "batch",
        "sheet_name",
        "source_row",
        "populated_cell_count",
        "populated_columns",
        "raw_a_preview",
        "exclusion_reason",
        "analysis_eligibility",
    ]
    write_tsv(
        metadata_dir / "source_structure_anomalies.tsv",
        anomaly_fields,
        anomaly_rows,
    )
    write_tsv(reusable_dir / "proteomics_patient_sample_map.tsv", mapping_fields, patient_map_rows)

    protein_fields = [
        "batch",
        "sheet_name",
        "source_row",
        *BASE_COLUMNS.values(),
    ]
    write_tsv(
        reusable_dir / "proteomics_protein_identification.tsv",
        protein_fields,
        protein_annotation_rows,
    )

    ratio_fields = [
        "batch",
        "sheet_name",
        "source_row",
        "paper_patient_id",
        "ratio_header",
        "normal_channel",
        "tumor_channel",
        "normal_vs_tumor_ratio",
        "log2_tumor_vs_normal",
        "quant_number",
        "source_significance",
        "cross_batch_reference",
    ]
    write_tsv(reusable_dir / "proteomics_ratio_long.tsv", ratio_fields, ratio_rows)

    provenance_fields = [
        "output_scope",
        "corresponding_source_file",
        "input_sha256",
        "generating_script",
        "key_parameters",
        "software",
        "generated_date",
        "regenerable",
    ]
    provenance_rows = [
        {
            "output_scope": "10_metadata/workbook_inventory.tsv;10_metadata/workbook_headers.tsv;10_metadata/source_structure_anomalies.tsv",
            "corresponding_source_file": f"00_source/{SOURCE_NAME}",
            "input_sha256": sha256_file(SOURCE_PATH),
            "generating_script": "scripts/09_ingest_cao2020_proteomics_supplement.py",
            "key_parameters": "artifact-tool workbook inspection followed by OOXML header and structure audit; nonnumeric GroupID rows retained as source anomalies and excluded from quantitation",
            "software": "@oai/artifact-tool 2.8.6+; Python standard library OOXML parser",
            "generated_date": DOWNLOAD_DATE,
            "regenerable": "yes",
        },
        {
            "output_scope": "10_metadata/proteomics_mapping_audit.tsv;20_reusable/proteomics_patient_sample_map.tsv",
            "corresponding_source_file": f"00_source/{SOURCE_NAME}",
            "input_sha256": sha256_file(SOURCE_PATH),
            "generating_script": "scripts/09_ingest_cao2020_proteomics_supplement.py",
            "key_parameters": "map explicit N-vs-T headers; retain P7 three-batch technical ratios; exclude Z1/Z2/Z3",
            "software": "Python standard library",
            "generated_date": DOWNLOAD_DATE,
            "regenerable": "yes",
        },
        {
            "output_scope": "20_reusable/proteomics_protein_identification.tsv;20_reusable/proteomics_ratio_long.tsv",
            "corresponding_source_file": f"00_source/{SOURCE_NAME}",
            "input_sha256": sha256_file(SOURCE_PATH),
            "generating_script": "scripts/09_ingest_cao2020_proteomics_supplement.py",
            "key_parameters": "separate one-row-per-source-protein annotations from patient ratios; source ratio=N/T; derived log2(T/N)=-log2(N/T); no protein-level deduplication at intake",
            "software": "Python standard library",
            "generated_date": DOWNLOAD_DATE,
            "regenerable": "yes",
        },
    ]
    write_tsv(reusable_dir / "PROVENANCE.tsv", provenance_fields, provenance_rows)

    dataset_lines = [
        "# Cao 2020 Supplementary Data 6 蛋白组比值说明",
        "",
        "## 基本信息",
        "",
        f"- `dataset_key`：`{DATASET_KEY}`。",
        "- 来源：Cao 等 2020，Nature Communications，Supplementary Data 6。",
        "- DOI：`10.1038/s41467-020-17227-z`。",
        "- 来源对象快照标签：`supplement_data6_snapshot_20231204`；实际获取日期：2026-07-11。",
        "- 物种：人（Homo sapiens）；疾病：食管鳞状细胞癌（ESCC）。",
        "- 内容：3 个 iTRAQ 批次工作表，分别 8,344、8,666、7,070 行和 51 列；无公式、无合并单元格。",
        "- 患者：P1、P2、P6、P7、P8、P9、P11、P12、P13、P15 共 10 对肿瘤/癌旁组织。",
        "- 当前完整性状态：`verified`；原工作簿 ZIP 完整，结构经电子表格工具与 OOXML 双重核查，生成文件已通过 SHA256。",
        "",
        "## 通道与患者映射",
        "",
        "- 批次 1：P1、P2、P6、P7；批次 2：P8、P9、P11、P7；批次 3：P12、P13、P15、P7。",
        "- P7 在三个批次均作为公共技术参考；保留 3 个技术比值，但不能计作 3 位独立患者。",
        "- 来源比值列为 `N-VS-T`，即癌旁/肿瘤；`proteomics_ratio_long.tsv` 同时保存来源 N/T 比值和派生 `log2(T/N)=-log2(N/T)`。",
        "- Z1、Z2、Z3 不是论文患者 ID，身份未明；只保留在映射审计表，不进入患者级蛋白矩阵。",
        "",
        "## 目录与主要文件",
        "",
        f"- `00_source/{SOURCE_NAME}`：出版商原始补充工作簿，不原位改写。",
        "- `10_metadata/workbook_inventory.tsv`：工作表可见性、维度、公式、合并单元格和表对象审计。",
        "- `10_metadata/workbook_headers.tsv`：3 个工作表的全部 51 列原始标题。",
        "- `10_metadata/proteomics_mapping_audit.tsv`：15 个比值组的患者、通道、证据、置信度和排除状态。",
        f"- `10_metadata/source_structure_anomalies.tsv`：{len(anomaly_rows)} 个未以数值 GroupID 开始的错列或续行来源行；保留审计但排除定量。",
        "- `20_reusable/proteomics_patient_sample_map.tsv`：12 个可用患者比值映射，其中 P7 跨 3 批次重复。",
        "- `20_reusable/proteomics_protein_identification.tsv`：每个来源蛋白行只保存一次的鉴定与注释。",
        "- `20_reusable/proteomics_ratio_long.tsv`：来源蛋白行展开为患者级长表；未在入库阶段强行去重蛋白。",
        "- `20_reusable/PROVENANCE.tsv`：生成规则和输入 SHA256。",
        "- `90_manifests/MANIFEST.tsv`：逐文件大小、SHA256、来源和完整性状态。",
        "",
        "## 推荐使用场景",
        "",
        "- Cao 2020 同患者甲基化—转录—蛋白候选轴的蛋白方向校准。",
        "- 先按批次和蛋白标识完成映射/缺失处理，再以患者为推断单位分析。",
        "- P7 可用于跨批次重复一致性检查，但不增加生物学样本量。",
        "",
        "## 不适用场景与限制",
        "",
        "- 本表是作者补充蛋白比值/鉴定表，不是原始质谱强度或完整 PRIDE 原始档案。",
        "- 三个批次的蛋白行、同源蛋白和缺失模式不同，不能直接把行号当作统一蛋白主键。",
        "- `Sig` 为来源字段，正式统计不得只依赖该标签；需保留效应量、患者级重复和批次敏感性。",
        "- Z1/Z2/Z3 身份未解析前不得用于患者级整合。",
        f"- 来源工作簿含 {len(anomaly_rows)} 个错列或长文本续行，不能按表头解释；正式矩阵已排除并保留审计。",
        "- 公共数据只支持跨层级候选轴，不证明因果机制。",
        "",
        "## 来源、引用与许可",
        "",
        "- 论文：https://www.nature.com/articles/s41467-020-17227-z",
        f"- 来源工作簿：{SOURCE_URL}",
        "- 使用时应引用原始论文并遵守出版商和 PRIDE 条款。",
        f"- 使用项目：`{PROJECT_ROOT}`。",
    ]
    (STAGE_ROOT / "DATASET.md").write_text("\n".join(dataset_lines) + "\n", encoding="utf-8")

    generated_files = [
        metadata_dir / "workbook_inventory.tsv",
        metadata_dir / "workbook_headers.tsv",
        metadata_dir / "proteomics_mapping_audit.tsv",
        metadata_dir / "source_structure_anomalies.tsv",
        reusable_dir / "proteomics_patient_sample_map.tsv",
        reusable_dir / "proteomics_protein_identification.tsv",
        reusable_dir / "proteomics_ratio_long.tsv",
        reusable_dir / "PROVENANCE.tsv",
    ]
    manifest_fields = [
        "relative_path",
        "file_level",
        "size_bytes",
        "sha256",
        "source_url",
        "download_date",
        "file_status",
        "corresponding_source_file",
        "generation_method",
        "notes",
    ]
    manifest_rows = [
        {
            "relative_path": f"00_source/{SOURCE_NAME}",
            "file_level": "00_source",
            "size_bytes": SOURCE_PATH.stat().st_size,
            "sha256": sha256_file(SOURCE_PATH),
            "source_url": SOURCE_URL,
            "download_date": DOWNLOAD_DATE,
            "file_status": "verified",
            "corresponding_source_file": "",
            "generation_method": "source bytes retrieved without modification",
            "notes": "xlsx ZIP integrity and workbook structure verified",
        }
    ]
    for path in generated_files:
        relative = path.relative_to(STAGE_ROOT).as_posix()
        manifest_rows.append(
            {
                "relative_path": relative,
                "file_level": relative.split("/", 1)[0],
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "source_url": "generated from Supplementary Data 6 workbook",
                "download_date": DOWNLOAD_DATE,
                "file_status": "generated_verified",
                "corresponding_source_file": f"00_source/{SOURCE_NAME}",
                "generation_method": "scripts/09_ingest_cao2020_proteomics_supplement.py",
                "notes": "",
            }
        )
    manifest_rows.sort(key=lambda row: row["relative_path"])
    write_tsv(manifest_dir / "MANIFEST.tsv", manifest_fields, manifest_rows)

    CANONICAL_ROOT.parent.mkdir(parents=True, exist_ok=True)
    if CANONICAL_ROOT.exists():
        if not REPAIR_EXISTING:
            raise RuntimeError(f"规范目录已存在，拒绝覆盖：{CANONICAL_ROOT}")
        repair_files = [
            STAGE_ROOT / "DATASET.md",
            manifest_dir / "MANIFEST.tsv",
            *generated_files,
        ]
        for source in repair_files:
            destination = CANONICAL_ROOT / source.relative_to(STAGE_ROOT)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
    else:
        shutil.copytree(STAGE_ROOT, CANONICAL_ROOT)
    canonical_manifest_fields, canonical_manifest = read_tsv(
        CANONICAL_ROOT / "90_manifests" / "MANIFEST.tsv"
    )
    if canonical_manifest_fields != manifest_fields:
        raise RuntimeError("规范路径 MANIFEST 字段变化")
    for row in canonical_manifest:
        path = CANONICAL_ROOT / row["relative_path"]
        if not path.is_file():
            raise RuntimeError(f"规范路径缺失文件：{path}")
        if path.stat().st_size != int(row["size_bytes"]):
            raise RuntimeError(f"规范路径大小不符：{path}")
        if sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"规范路径 SHA256 不符：{path}")

    catalog_path = DATA_ROOT / "CATALOG.tsv"
    catalog_fields, catalog_rows = read_tsv(catalog_path)
    catalog_record = {
        "dataset_key": DATASET_KEY,
        "record_type": "publication_supplement_proteomics_ratio_dataset",
        "access_level": "public",
        "source": "Nature Communications",
        "accession": "10.1038/s41467-020-17227-z;Supplementary Data 6;PXD019834",
        "version": "supplement_data6_snapshot_20231204",
        "species": "Homo sapiens",
        "disease": "esophageal squamous cell carcinoma",
        "tissue": "primary tumor; adjacent normal esophagus",
        "assay": "iTRAQ proteomics; author supplementary N-vs-T ratio tables",
        "sample_summary": "3 iTRAQ batches; 10 patient pairs; P7 repeated in all 3 batches as technical reference; Z1/Z2/Z3 excluded",
        "available_levels": "00_source;10_metadata;20_reusable",
        "status": "verified",
        "local_path": str(CANONICAL_ROOT),
        "manifest_path": str(CANONICAL_ROOT / "90_manifests" / "MANIFEST.tsv"),
        "source_url": SOURCE_URL,
        "download_date": DOWNLOAD_DATE,
        "last_verified": DOWNLOAD_DATE,
        "license_or_access": "public publication supplement; cite Cao et al. 2020 and follow source terms",
        "projects_using": str(PROJECT_ROOT),
        "recommended_use": "matched-patient protein direction for Cao 2020 methylation-RNA-protein candidate axes",
        "limitations": (
            "ratio table not raw intensities; batch-specific protein rows; "
            f"P7 technical repeats; Z ratios unresolved; {len(anomaly_rows)} "
            "structurally shifted or continuation rows excluded; PXD archive not included"
        ),
        "notes": (
            "Workbook headers resolve all 10 patient ratios; Z1/Z2/Z3 and "
            f"{len(anomaly_rows)} structurally anomalous source rows remain excluded."
        ),
    }
    catalog_matches = [
        index
        for index, row in enumerate(catalog_rows)
        if row["dataset_key"] == DATASET_KEY
    ]
    if len(catalog_matches) > 1:
        raise RuntimeError("CATALOG.tsv 存在重复 dataset_key")
    if catalog_matches:
        if not REPAIR_EXISTING:
            raise RuntimeError("CATALOG.tsv 已存在 dataset_key")
        catalog_rows[catalog_matches[0]].update(catalog_record)
    else:
        catalog_rows.append(catalog_record)
    write_tsv_atomic(catalog_path, catalog_fields, catalog_rows)

    project_dataset_path = PROJECT_ROOT / "data" / "datasets.tsv"
    project_fields, project_rows = read_tsv(project_dataset_path)
    project_record = {
        "logical_name": "cao2020_proteomics_supplement_data6",
        "dataset_key": DATASET_KEY,
        "accession": "10.1038/s41467-020-17227-z;Supplementary Data 6;PXD019834",
        "version": "supplement_data6_snapshot_20231204",
        "data_level": "00_source;10_metadata;20_reusable",
        "central_path": str(CANONICAL_ROOT),
        "project_purpose": "Cao 2020 同患者甲基化-RNA-蛋白候选轴的蛋白方向校准",
        "inclusion_exclusion": (
            "纳入 P1/P2/P6/P7/P8/P9/P11/P12/P13/P15；P7 三批次按技术重复处理；"
            f"排除 Z1/Z2/Z3 和 {len(anomaly_rows)} 个结构异常来源行"
        ),
        "limitations": (
            "作者比值表而非原始强度；需统一蛋白标识和批次缺失；"
            f"{len(anomaly_rows)} 个错列或续行已排除；PXD 原始档案尚未纳入"
        ),
        "local_link": "",
        "last_verified": DOWNLOAD_DATE,
    }
    project_matches = [
        index
        for index, row in enumerate(project_rows)
        if row["dataset_key"] == DATASET_KEY
    ]
    if len(project_matches) > 1:
        raise RuntimeError("data/datasets.tsv 存在重复 dataset_key")
    if project_matches:
        if not REPAIR_EXISTING:
            raise RuntimeError("data/datasets.tsv 已存在 dataset_key")
        project_rows[project_matches[0]].update(project_record)
    else:
        project_rows.append(project_record)
    write_tsv_atomic(project_dataset_path, project_fields, project_rows)

    print(
        "已完成 Cao 2020 Supplementary Data 6 的结构审计、10 位患者通道映射、"
        "蛋白比值长表、规范路径复制、全 SHA 校验和索引登记。"
    )


if __name__ == "__main__":
    main()
