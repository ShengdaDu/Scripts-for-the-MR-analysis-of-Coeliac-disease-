#!/usr/bin/env bash

# ==============================================================================
# GWAS Catalog harmonised 暴露数据并行下载脚本（安全停止版）
#
# 主要功能：
#   1. 从 metadata 中读取 exposure_id、file_name 和 harmonised_directory_url
#   2. 自动解析每个 harmonised 目录中唯一的 *.h.tsv.gz 文件
#   3. 同时下载多个暴露文件
#   4. 使用 *.part 文件支持断点续传
#   5. 下载完成后使用 gzip -t 进行完整性校验
#   6. 使用 *.validated 标记，避免重复扫描已经验证的完整 gzip 文件
#   7. 使用 flock 保护并发日志写入
#   8. 将 xargs、子 shell、wget 和 gzip 放入同一个独立进程组
#   9. 主脚本收到 Ctrl+C、TERM 或 HUP 时，自动终止整个任务进程组
#  10. 即使主脚本被 kill -9，独立 watchdog 仍会清理残留任务
#  11. 使用单实例锁，避免重复启动多个相同下载任务
#
# 推荐 metadata：
#   ./config/exposure_metadata_remaining_565_extended.tsv
# ==============================================================================

set -u
set -o pipefail

# ------------------------------- 用户设置 --------------------------------------
METADATA="./config/exposure_metadata_remaining_565_extended.tsv"
OUTDIR="./exposure"
LOGDIR="./download_logs/exposure_download_remaining_565"

# 同时运行的下载任务数。建议从 4 开始；网络和磁盘稳定后可改为 6。
MAX_JOBS=4

# wget 重试参数
TRIES=5
TIMEOUT=60
WAIT_RETRY=10
# ------------------------------------------------------------------------------

mkdir -p "${OUTDIR}" "${LOGDIR}"

SUCCESS_LOG="${LOGDIR}/download_success.tsv"
FAILED_LOG="${LOGDIR}/download_failed.tsv"
URL_LOG="${LOGDIR}/resolved_download_urls.tsv"
TASK_FILE="${LOGDIR}/download_tasks.tsv"
TASK_NUL_FILE="${LOGDIR}/download_tasks.nul"

# 单实例锁和 PID 信息文件
LOCK_FILE="${LOGDIR}/downloader.lock"
PID_FILE="${LOGDIR}/downloader.pid"

# xargs 下载进程组和 watchdog 的 PID
# 使用 setsid 启动后，XARGS_PID 同时也是整个下载任务组的 PGID。
XARGS_PID=""
WATCHDOG_PID=""

# ------------------------------------------------------------------------------
# 检查依赖命令
# ------------------------------------------------------------------------------
for cmd in wget gzip awk grep sed sort xargs flock setsid; do
    command -v "${cmd}" >/dev/null 2>&1 || {
        echo "错误：未安装必要命令：${cmd}" >&2
        exit 1
    }
done

if [[ ! -f "${METADATA}" ]]; then
    echo "错误：找不到 metadata 文件：${METADATA}" >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# 单实例锁
#
# 同一时间只允许一个脚本实例运行，避免两个实例同时下载同一批文件。
# 主脚本退出后，文件描述符 200 会自动关闭并释放锁。
# ------------------------------------------------------------------------------
exec 200>"${LOCK_FILE}"

if ! flock -n 200; then
    echo "错误：检测到另一个下载脚本实例正在运行。" >&2
    echo "锁文件：${LOCK_FILE}" >&2
    echo "PID 文件：${PID_FILE}" >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# 初始化本次运行日志
# ------------------------------------------------------------------------------
printf "exposure_id\tlocal_file\tremote_file\tdownload_url\tstatus\n" \
    > "${SUCCESS_LOG}"
printf "exposure_id\tlocal_file\tharmonised_directory_url\treason\n" \
    > "${FAILED_LOG}"
printf "exposure_id\tremote_file\tdownload_url\n" \
    > "${URL_LOG}"

# ------------------------------------------------------------------------------
# 从 metadata 生成任务表
# ------------------------------------------------------------------------------
awk -F '\t' '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        gsub(/\r$/, "", $i)
        col[$i] = i
    }

    required[1] = "include"
    required[2] = "exposure_id"
    required[3] = "file_name"
    required[4] = "harmonised_directory_url"

    for (j = 1; j <= 4; j++) {
        if (!(required[j] in col)) {
            print "错误：metadata 缺少列：" required[j] > "/dev/stderr"
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
        local_name = $col["file_name"]
        harmonised_url = $col["harmonised_directory_url"]

        gsub(/\r$/, "", exposure_id)
        gsub(/\r$/, "", local_name)
        gsub(/\r$/, "", harmonised_url)

        if (exposure_id != "" && local_name != "" && harmonised_url != "") {
            print exposure_id "\t" local_name "\t" harmonised_url
        }
    }
}
' "${METADATA}" > "${TASK_FILE}"

if [[ ! -s "${TASK_FILE}" ]]; then
    echo "错误：metadata 中没有 include=True 的有效任务。" >&2
    exit 1
fi

# 预先生成 NUL 分隔的任务输入文件。
# 不再使用“while | xargs”管道，减少额外的管道子进程和进程树复杂度。
: > "${TASK_NUL_FILE}"
while IFS=$'\t' read -r exposure_id local_name harmonised_url; do
    printf '%s\0%s\0%s\0' \
        "${exposure_id}" "${local_name}" "${harmonised_url}" \
        >> "${TASK_NUL_FILE}"
done < "${TASK_FILE}"

echo "Metadata 文件：${METADATA}"
echo "任务数量     ：$(wc -l < "${TASK_FILE}")"
echo "并发任务数   ：${MAX_JOBS}"
echo "输出目录     ：${OUTDIR}"
echo

# ------------------------------------------------------------------------------
# 下载单个暴露文件
# ------------------------------------------------------------------------------
download_one() {
    local exposure_id="$1"
    local local_name="$2"
    local harmonised_url="$3"

    local local_file="${OUTDIR}/${local_name}"
    local partial_file="${local_file}.part"
    local marker_file="${local_file}.validated"

    local listing=""
    local remote_file=""
    local download_url=""
    local wget_status=0
    local -a candidates=()

    echo "[${exposure_id}] 开始处理"

    # 已下载且已经验证的文件直接跳过，不再重复运行 gzip -t。
    if [[ -s "${local_file}" && -f "${marker_file}" ]]; then
        echo "[${exposure_id}] 文件已有 validated 标记，跳过"

        {
            flock 9
            printf "%s\t%s\tNA\tNA\talready_present\n" \
                "${exposure_id}" "${local_name}"
        } 9>>"${SUCCESS_LOG}"

        return 0
    fi

    # 最终文件存在但没有 validated 标记时，只进行一次完整校验。
    if [[ -s "${local_file}" ]]; then
        if gzip -t "${local_file}" 2>/dev/null; then
            touch "${marker_file}"
            echo "[${exposure_id}] 现有文件校验通过，创建 validated 标记"

            {
                flock 9
                printf "%s\t%s\tNA\tNA\texisting_file_validated\n" \
                    "${exposure_id}" "${local_name}"
            } 9>>"${SUCCESS_LOG}"

            return 0
        else
            echo "[${exposure_id}] 现有最终文件损坏，删除后重新下载"
            rm -f "${local_file}" "${marker_file}"
        fi
    fi

    # 获取 harmonised 目录页面。
    if ! listing="$(
        wget -qO- \
            --tries="${TRIES}" \
            --timeout="${TIMEOUT}" \
            --waitretry="${WAIT_RETRY}" \
            --retry-connrefused \
            "${harmonised_url}"
    )"; then
        listing=""
    fi

    if [[ -z "${listing}" ]]; then
        echo "[${exposure_id}] 无法读取 harmonised 目录"

        {
            flock 9
            printf "%s\t%s\t%s\tdirectory_listing_failed\n" \
                "${exposure_id}" "${local_name}" "${harmonised_url}"
        } 9>>"${FAILED_LOG}"

        return 1
    fi

    # 精确提取以 .h.tsv.gz 结尾的链接，不匹配 meta.yaml 等附属文件。
    mapfile -t candidates < <(
        printf "%s" "${listing}" |
            grep -oE 'href="[^"]+\.h\.tsv\.gz"' |
            sed -E 's/^href="//; s/"$//' |
            sort -u
    )

    if [[ "${#candidates[@]}" -eq 0 ]]; then
        echo "[${exposure_id}] 目录中没有找到 *.h.tsv.gz 文件"

        {
            flock 9
            printf "%s\t%s\t%s\tno_harmonised_h_file\n" \
                "${exposure_id}" "${local_name}" "${harmonised_url}"
        } 9>>"${FAILED_LOG}"

        return 1
    fi

    if [[ "${#candidates[@]}" -gt 1 ]]; then
        echo "[${exposure_id}] 找到多个 *.h.tsv.gz 文件，拒绝模糊下载："
        printf "  %s\n" "${candidates[@]}"

        {
            flock 9
            printf "%s\t%s\t%s\tmultiple_harmonised_h_files\n" \
                "${exposure_id}" "${local_name}" "${harmonised_url}"
        } 9>>"${FAILED_LOG}"

        return 1
    fi

    remote_file="${candidates[0]}"

    if [[ "${remote_file}" =~ ^https?:// ]]; then
        download_url="${remote_file}"
    else
        download_url="${harmonised_url%/}/${remote_file}"
    fi

    {
        flock 9
        printf "%s\t%s\t%s\n" \
            "${exposure_id}" "${remote_file}" "${download_url}"
    } 9>>"${URL_LOG}"

    echo "[${exposure_id}] 开始下载：${remote_file}"

    # 下载到 .part 文件，脚本重启后可继续断点续传。
    wget \
        --continue \
        --tries="${TRIES}" \
        --timeout="${TIMEOUT}" \
        --waitretry="${WAIT_RETRY}" \
        --retry-connrefused \
        --no-verbose \
        --output-document="${partial_file}" \
        "${download_url}"

    wget_status=$?

    if [[ "${wget_status}" -ne 0 ]]; then
        echo "[${exposure_id}] 下载失败，保留 .part 文件供下次续传"

        {
            flock 9
            printf "%s\t%s\t%s\twget_failed_%s\n" \
                "${exposure_id}" "${local_name}" \
                "${harmonised_url}" "${wget_status}"
        } 9>>"${FAILED_LOG}"

        return 1
    fi

    # 下载完成后校验压缩文件；通过后才改为最终文件名。
    if gzip -t "${partial_file}" 2>/dev/null; then
        mv -f "${partial_file}" "${local_file}"
        touch "${marker_file}"

        echo "[${exposure_id}] 下载完成且 gzip 校验通过"

        {
            flock 9
            printf "%s\t%s\t%s\t%s\tdownloaded\n" \
                "${exposure_id}" "${local_name}" \
                "${remote_file}" "${download_url}"
        } 9>>"${SUCCESS_LOG}"

        return 0
    else
        echo "[${exposure_id}] gzip 校验失败，保留 .part 文件"

        {
            flock 9
            printf "%s\t%s\t%s\tgzip_validation_failed\n" \
                "${exposure_id}" "${local_name}" "${harmonised_url}"
        } 9>>"${FAILED_LOG}"

        return 1
    fi
}

export -f download_one
export OUTDIR SUCCESS_LOG FAILED_LOG URL_LOG
export TRIES TIMEOUT WAIT_RETRY

# ------------------------------------------------------------------------------
# 判断完整下载任务进程组是否仍然存在
# ------------------------------------------------------------------------------
worker_group_alive() {
    [[ -n "${XARGS_PID}" ]] &&
        kill -0 -- "-${XARGS_PID}" 2>/dev/null
}

# ------------------------------------------------------------------------------
# 终止整个下载任务进程组
#
# xargs 使用 setsid 启动，因此 XARGS_PID 同时是 PGID。
# 向负 PGID 发送信号会同时终止：
#   xargs、bash -c download_one、wget、gzip -t
# ------------------------------------------------------------------------------
stop_worker_group() {
    if ! worker_group_alive; then
        return 0
    fi

    echo
    echo "收到停止信号，正在终止全部下载和校验任务……"

    kill -TERM -- "-${XARGS_PID}" 2>/dev/null || true

    # 最多等待 5 秒，让 wget 和 gzip 正常退出。
    local i
    for i in {1..10}; do
        if ! worker_group_alive; then
            break
        fi
        sleep 0.5
    done

    # 仍有残留时，强制终止整个任务进程组。
    if worker_group_alive; then
        echo "部分任务没有正常退出，执行强制终止……"
        kill -KILL -- "-${XARGS_PID}" 2>/dev/null || true
    fi

    wait "${XARGS_PID}" 2>/dev/null || true
}

# ------------------------------------------------------------------------------
# 主脚本退出清理函数
# ------------------------------------------------------------------------------
cleanup() {
    local exit_code="${1:-$?}"

    # 防止 cleanup 内部再次触发相同 trap。
    trap - INT TERM HUP EXIT

    stop_worker_group

    # 正常终止主脚本时，同时关闭 watchdog。
    if [[ -n "${WATCHDOG_PID}" ]] &&
       kill -0 "${WATCHDOG_PID}" 2>/dev/null; then
        kill -TERM "${WATCHDOG_PID}" 2>/dev/null || true
        wait "${WATCHDOG_PID}" 2>/dev/null || true
    fi

    rm -f "${PID_FILE}"
    exit "${exit_code}"
}

# Ctrl+C、普通 kill、终端断开及普通退出均执行清理。
trap 'cleanup 130' INT
trap 'cleanup 143' TERM
trap 'cleanup 129' HUP
trap 'cleanup $?' EXIT

# ------------------------------------------------------------------------------
# 在独立 session/进程组中启动 xargs 并行调度器
#
# 这是防止“主脚本被 kill 后，wget/gzip 仍继续运行”的核心修改。
# ------------------------------------------------------------------------------
setsid xargs \
    -0 \
    -r \
    -n 3 \
    -P "${MAX_JOBS}" \
    bash -c 'download_one "$1" "$2" "$3"' _ \
    < "${TASK_NUL_FILE}" &

XARGS_PID=$!

# ------------------------------------------------------------------------------
# 启动独立 watchdog
#
# SIGKILL（kill -9）无法被主脚本 trap 捕获。
# watchdog 每 2 秒检查主脚本是否存在；主脚本消失后，它会终止整个下载任务组。
# watchdog 本身通过 setsid 启动，不属于下载任务组，因此不会被提前一起杀掉。
# ------------------------------------------------------------------------------
setsid bash -c '
    main_pid="$1"
    worker_pgid="$2"
    pid_file="$3"

    while kill -0 "${main_pid}" 2>/dev/null; do
        sleep 2
    done

    if kill -0 -- "-${worker_pgid}" 2>/dev/null; then
        kill -TERM -- "-${worker_pgid}" 2>/dev/null || true
        sleep 3

        if kill -0 -- "-${worker_pgid}" 2>/dev/null; then
            kill -KILL -- "-${worker_pgid}" 2>/dev/null || true
        fi
    fi

    rm -f "${pid_file}"
' _ "$$" "${XARGS_PID}" "${PID_FILE}" >/dev/null 2>&1 &

WATCHDOG_PID=$!

# 记录主脚本、工作进程组和 watchdog 的 PID，便于查看和手动停止。
printf "main_pid\t%s\nworker_pgid\t%s\nwatchdog_pid\t%s\n" \
    "$$" "${XARGS_PID}" "${WATCHDOG_PID}" > "${PID_FILE}"

echo "主脚本 PID       ：$$"
echo "任务进程组 PGID  ：${XARGS_PID}"
echo "Watchdog PID     ：${WATCHDOG_PID}"
echo "PID 信息文件     ：${PID_FILE}"
echo
echo "推荐停止方式："
echo "  kill -TERM $$"
echo "脚本会自动停止 xargs、子 shell、wget 和 gzip -t。"
echo

# 等待全部并行任务结束。
wait "${XARGS_PID}"
xargs_status=$?

# 正常完成后清空变量，避免 EXIT trap 重复终止已经结束的进程组。
XARGS_PID=""

if [[ -n "${WATCHDOG_PID}" ]] &&
   kill -0 "${WATCHDOG_PID}" 2>/dev/null; then
    kill -TERM "${WATCHDOG_PID}" 2>/dev/null || true
    wait "${WATCHDOG_PID}" 2>/dev/null || true
fi
WATCHDOG_PID=""

rm -f "${PID_FILE}"

# 正常完成后关闭 trap，再输出汇总。
trap - INT TERM HUP EXIT

echo
echo "================================================================"
echo "下载任务结束。"
echo "成功记录数：$(( $(wc -l < "${SUCCESS_LOG}") - 1 ))"
echo "失败记录数：$(( $(wc -l < "${FAILED_LOG}") - 1 ))"
echo "成功日志  ：${SUCCESS_LOG}"
echo "失败日志  ：${FAILED_LOG}"
echo "URL 日志  ：${URL_LOG}"
echo "================================================================"

exit "${xargs_status}"
