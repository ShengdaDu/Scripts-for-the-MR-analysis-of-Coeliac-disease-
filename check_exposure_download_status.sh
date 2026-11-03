#!/usr/bin/env bash

# ==============================================================================
# 根据metadata核对GWAS暴露文件下载状态
#
# 功能：
#   1. 以metadata中include=True的记录为预期下载清单
#   2. 检查最终文件、.part文件和.validated标记
#   3. 区分下载成功、文件缺失、仅有部分文件、空文件和gzip损坏
#   4. 生成完整状态表、成功表、失败表和仅含失败项的重下载metadata
#   5. 统计exposure目录中不属于本metadata的额外文件
#
# 默认使用快速模式：
#   只要最终*.h.tsv.gz文件存在且大小大于0，即计为成功。
#
# 严格模式：
#   CHECK_GZIP=1 ./check_exposure_download_status.sh
#   会对没有.validated标记的最终文件执行gzip -t。
#   由于文件较大，严格模式可能耗时较长。
# ==============================================================================

set -u
set -o pipefail

# ------------------------------- 用户设置 --------------------------------------
METADATA="${METADATA:-./config/exposure_metadata_remaining_565_extended.tsv}"
OUTDIR="${OUTDIR:-./exposure}"
REPORT_DIR="${REPORT_DIR:-./download_logs/exposure_status_check}"

# 0：快速检查文件是否存在且非空
# 1：对无.validated标记的最终文件执行gzip -t
CHECK_GZIP="${CHECK_GZIP:-0}"
# ------------------------------------------------------------------------------

mkdir -p "${REPORT_DIR}"

ALL_STATUS="${REPORT_DIR}/download_status_all.tsv"
SUCCESS_STATUS="${REPORT_DIR}/download_success.tsv"
FAILED_STATUS="${REPORT_DIR}/download_failed.tsv"
SUMMARY_FILE="${REPORT_DIR}/download_summary.tsv"
RETRY_METADATA="${REPORT_DIR}/exposure_metadata_retry_failed.tsv"
EXTRA_FILES="${REPORT_DIR}/extra_files_not_in_metadata.tsv"

TASK_FILE="${REPORT_DIR}/expected_download_tasks.tsv"
SUCCESS_IDS="${REPORT_DIR}/success_exposure_ids.txt"
FAILED_IDS="${REPORT_DIR}/failed_exposure_ids.txt"
EXPECTED_NAMES="${REPORT_DIR}/expected_file_names.txt"
ACTUAL_NAMES="${REPORT_DIR}/actual_final_file_names.txt"

if [[ ! -f "${METADATA}" ]]; then
    echo "错误：找不到metadata文件：${METADATA}" >&2
    exit 1
fi

if [[ ! -d "${OUTDIR}" ]]; then
    echo "错误：找不到暴露文件目录：${OUTDIR}" >&2
    exit 1
fi

for cmd in awk gzip find sort comm wc stat cut; do
    command -v "${cmd}" >/dev/null 2>&1 || {
        echo "错误：缺少必要命令：${cmd}" >&2
        exit 1
    }
done

if [[ "${CHECK_GZIP}" != "0" && "${CHECK_GZIP}" != "1" ]]; then
    echo "错误：CHECK_GZIP只能设置为0或1。" >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# 根据列名读取metadata，避免依赖固定列顺序
# ------------------------------------------------------------------------------
awk -F '\t' '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        gsub(/\r$/, "", $i)
        col[$i] = i
    }

    required[1] = "include"
    required[2] = "exposure_id"
    required[3] = "exposure_name"
    required[4] = "file_name"

    for (j = 1; j <= 4; j++) {
        if (!(required[j] in col)) {
            print "错误：metadata缺少列：" required[j] > "/dev/stderr"
            exit 2
        }
    }
    next
}
{
    include = tolower($col["include"])
    gsub(/\r$/, "", include)

    if (include == "true" || include == "1" ||
        include == "yes"  || include == "y") {

        exposure_id = $col["exposure_id"]
        exposure_name = $col["exposure_name"]
        file_name = $col["file_name"]

        gsub(/\r$/, "", exposure_id)
        gsub(/\r$/, "", exposure_name)
        gsub(/\r$/, "", file_name)

        if (exposure_id != "" && file_name != "") {
            print exposure_id "\t" exposure_name "\t" file_name
        }
    }
}
' "${METADATA}" > "${TASK_FILE}"

if [[ ! -s "${TASK_FILE}" ]]; then
    echo "错误：metadata中没有include=True的有效任务。" >&2
    exit 1
fi

cut -f3 "${TASK_FILE}" | sort -u > "${EXPECTED_NAMES}"
: > "${SUCCESS_IDS}"
: > "${FAILED_IDS}"

printf "exposure_id\texposure_name\tfile_name\tlocal_path\tfile_size_bytes\tpart_size_bytes\tvalidated_marker\tstatus\treason\n" \
    > "${ALL_STATUS}"
printf "exposure_id\texposure_name\tfile_name\tlocal_path\tfile_size_bytes\tstatus\n" \
    > "${SUCCESS_STATUS}"
printf "exposure_id\texposure_name\tfile_name\tlocal_path\tpart_path\tfile_size_bytes\tpart_size_bytes\tstatus\treason\n" \
    > "${FAILED_STATUS}"

expected_count=0
success_count=0
failed_count=0
marker_success_count=0
gzip_success_count=0
unchecked_success_count=0
missing_count=0
partial_count=0
empty_count=0
corrupt_count=0
orphan_marker_count=0

echo "开始核对下载状态……"
echo "Metadata：${METADATA}"
echo "文件目录：${OUTDIR}"
echo "gzip严格校验：${CHECK_GZIP}"
echo

while IFS=$'\t' read -r exposure_id exposure_name file_name; do
    expected_count=$((expected_count + 1))

    local_file="${OUTDIR}/${file_name}"
    part_file="${local_file}.part"
    marker_file="${local_file}.validated"

    file_size=0
    part_size=0
    marker="no"
    status=""
    reason=""
    success=0

    if [[ -e "${local_file}" ]]; then
        file_size="$(stat -c '%s' "${local_file}" 2>/dev/null || printf '0')"
    fi

    if [[ -e "${part_file}" ]]; then
        part_size="$(stat -c '%s' "${part_file}" 2>/dev/null || printf '0')"
    fi

    if [[ -f "${marker_file}" ]]; then
        marker="yes"
    fi

    if [[ -s "${local_file}" ]]; then
        if [[ -f "${marker_file}" ]]; then
            status="success_validated_marker"
            reason="最终文件存在且已有validated标记"
            marker_success_count=$((marker_success_count + 1))
            success=1

        elif [[ "${CHECK_GZIP}" == "1" ]]; then
            if gzip -t "${local_file}" 2>/dev/null; then
                status="success_gzip_valid"
                reason="最终文件存在且gzip完整性校验通过"
                gzip_success_count=$((gzip_success_count + 1))
                success=1
            else
                status="failed_gzip_invalid"
                reason="最终文件存在，但gzip完整性校验失败"
                corrupt_count=$((corrupt_count + 1))
            fi

        else
            status="success_exists_unchecked"
            reason="最终文件存在且非空；未执行gzip完整性校验"
            unchecked_success_count=$((unchecked_success_count + 1))
            success=1
        fi

    elif [[ -e "${local_file}" ]]; then
        status="failed_empty_final"
        reason="最终文件存在但大小为0"
        empty_count=$((empty_count + 1))

    elif [[ -s "${part_file}" ]]; then
        status="failed_partial_only"
        reason="最终文件缺失，仅存在未完成的.part文件"
        partial_count=$((partial_count + 1))

    elif [[ -e "${part_file}" ]]; then
        status="failed_empty_partial"
        reason=".part文件存在但大小为0"
        empty_count=$((empty_count + 1))

    elif [[ -f "${marker_file}" ]]; then
        status="failed_orphan_marker"
        reason="存在validated标记，但最终文件缺失"
        orphan_marker_count=$((orphan_marker_count + 1))

    else
        status="failed_missing"
        reason="最终文件和.part文件均不存在"
        missing_count=$((missing_count + 1))
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${exposure_id}" "${exposure_name}" "${file_name}" "${local_file}" \
        "${file_size}" "${part_size}" "${marker}" "${status}" "${reason}" \
        >> "${ALL_STATUS}"

    if [[ "${success}" -eq 1 ]]; then
        success_count=$((success_count + 1))
        printf "%s\n" "${exposure_id}" >> "${SUCCESS_IDS}"

        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${exposure_id}" "${exposure_name}" "${file_name}" \
            "${local_file}" "${file_size}" "${status}" \
            >> "${SUCCESS_STATUS}"
    else
        failed_count=$((failed_count + 1))
        printf "%s\n" "${exposure_id}" >> "${FAILED_IDS}"

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${exposure_id}" "${exposure_name}" "${file_name}" \
            "${local_file}" "${part_file}" "${file_size}" "${part_size}" \
            "${status}" "${reason}" \
            >> "${FAILED_STATUS}"
    fi

    if (( expected_count % 50 == 0 )); then
        echo "已检查 ${expected_count} 个表型……"
    fi
done < "${TASK_FILE}"

sort -u -o "${SUCCESS_IDS}" "${SUCCESS_IDS}"
sort -u -o "${FAILED_IDS}" "${FAILED_IDS}"

# ------------------------------------------------------------------------------
# 生成只包含失败表型的重下载metadata
# 保留原metadata全部列，可直接交给下载脚本。
# ------------------------------------------------------------------------------
awk -F '\t' '
NR == FNR {
    failed[$1] = 1
    next
}
FNR == 1 {
    for (i = 1; i <= NF; i++) {
        header = $i
        gsub(/\r$/, "", header)
        col[header] = i
    }

    if (!("exposure_id" in col)) {
        print "错误：metadata缺少exposure_id列" > "/dev/stderr"
        exit 2
    }

    print
    next
}
{
    exposure_id = $col["exposure_id"]
    gsub(/\r$/, "", exposure_id)

    if (exposure_id in failed) {
        print
    }
}
' "${FAILED_IDS}" "${METADATA}" > "${RETRY_METADATA}"

# ------------------------------------------------------------------------------
# 查找exposure目录中存在、但不属于本metadata的额外最终文件
# ------------------------------------------------------------------------------
find "${OUTDIR}" -maxdepth 1 -type f -name '*.h.tsv.gz' -printf '%f\n' |
    sort -u > "${ACTUAL_NAMES}"

comm -13 "${EXPECTED_NAMES}" "${ACTUAL_NAMES}" > "${EXTRA_FILES}"
extra_count="$(wc -l < "${EXTRA_FILES}")"

# ------------------------------------------------------------------------------
# 输出汇总统计
# ------------------------------------------------------------------------------
{
    printf "metric\tcount\n"
    printf "metadata_expected\t%s\n" "${expected_count}"
    printf "download_success\t%s\n" "${success_count}"
    printf "download_failed\t%s\n" "${failed_count}"
    printf "success_validated_marker\t%s\n" "${marker_success_count}"
    printf "success_gzip_valid\t%s\n" "${gzip_success_count}"
    printf "success_exists_unchecked\t%s\n" "${unchecked_success_count}"
    printf "failed_missing\t%s\n" "${missing_count}"
    printf "failed_partial_only\t%s\n" "${partial_count}"
    printf "failed_empty_file\t%s\n" "${empty_count}"
    printf "failed_gzip_invalid\t%s\n" "${corrupt_count}"
    printf "failed_orphan_marker\t%s\n" "${orphan_marker_count}"
    printf "extra_final_files_not_in_metadata\t%s\n" "${extra_count}"
} > "${SUMMARY_FILE}"

echo
echo "================================================================"
echo "核对完成。"
echo "Metadata预期文件数 ：${expected_count}"
echo "下载成功数         ：${success_count}"
echo "下载失败/未完成数  ：${failed_count}"
echo "额外最终文件数     ：${extra_count}"
echo
echo "完整状态表         ：${ALL_STATUS}"
echo "成功文件表         ：${SUCCESS_STATUS}"
echo "失败文件表         ：${FAILED_STATUS}"
echo "重下载metadata     ：${RETRY_METADATA}"
echo "汇总统计           ：${SUMMARY_FILE}"
echo "额外文件表         ：${EXTRA_FILES}"
echo "================================================================"
