#負責把三個benchmark工具、對每種資料大小，依序自動跑一輪
#
#每種細胞數量子樣本，把三個 benchmark 腳本各自當成獨立、全新的程序執行，
#這樣峰值記憶體的量測才不會被前一次執行的殘留污染。三個工具各自用自己
#原生的方法選高變異基因(不再有共用的precompute前置步驟)，比較貼近
#真實使用情境，也不會因為precompute失敗而連累三個工具都沒有資料
#
#每次執行只跑「一輪」，不會自動連跑 N_ROUNDS 次——因為一次連續跑滿3輪
#耗時太長(尤其大尺寸的Seurat)，電腦要連續負荷太久不實際。改成每次執行
#自動判斷「這是第幾輪」:檢查 benchmark_results_round1.csv、round2.csv、
#round3.csv 裡哪一個還不存在，就把這次的結果寫進那一個。也就是說，你
#只要照樣執行同一個指令(bash run_all.sh)三次(可以分開幾天跑，中間電腦
#正常使用也沒關係)，它會自動接續成round1、round2、round3。全部3輪都
#跑完後，再執行會直接印出提示、不會覆蓋掉任何資料。
#跑完3輪後用 average_results.py 把幾輪的結果平均起來，讓數字比較有公信力
#(單次量測容易受系統當下的背景負載影響，多輪平均可以把這種雜訊降低)
#
#不管是前置步驟(precompute)失敗，還是三個工具裡任何一個失敗，都不會讓
#整個腳本中斷，也不會讓失敗「完全不留痕跡」——每一筆結果(不管成功或
#失敗)都會寫進 CSV，用 status 欄位(OK/FAILED)標記，失敗時 elapsed_sec、
#peak_memory_MB 留空。這樣之後可以直接在表格裡看到「這個工具在這個規模
#下失敗了」，不用去終端機翻紀錄
#
#每個工具呼叫都用 timeout 包住(見 TOOL_TIMEOUT_SEC)，超過上限就強制判定
#失敗、繼續跑下一個，不用整晚坐在電腦前盯著有沒有卡死。這個上限只是抓
#「明顯卡住」的安全網，不是效能門檻，所以設得比正常執行時間寬鬆很多；
#如果某個尺寸的正常執行時間真的比預設上限長，可以用環境變數覆寫，
#例如 TOOL_TIMEOUT_SEC=7200 bash run_all.sh。
#
#子樣本資料夾是照數字大小(由小到大)依序執行，不是字母順序——如果照字母
#順序排，"n_100000" 會排在 "n_20000" 前面(字串比較"1"<"2")，導致最大的
#那筆(n_2047431，記憶體風險最高)提早出現在中段，而不是像預期一樣留到
#最後。改成數字排序後，風險最高的尺寸永遠是最後一個，前面所有尺寸的結果
#都已經安全寫進CSV，不會被它拖累

#管線中有一環節出錯，終止腳本
set -euo pipefail

#算出路徑
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
BENCHMARK_DATA_DIR="$PROJECT_ROOT/results/benchmark_data"
RESULTS_DIR="$PROJECT_ROOT/results"

N_ROUNDS=3

#每個工具呼叫的安全逾時上限(秒)，可用環境變數覆寫。60分鐘對目前已知
#的尺寸(到n_1500000為止)都綽綽有餘，204萬筆那個尺寸如果需要更寬鬆，
#執行時另外指定 TOOL_TIMEOUT_SEC 即可，不用改這裡
TOOL_TIMEOUT_SEC="${TOOL_TIMEOUT_SEC:-3600}"

#防呆:確認03輸出資料夾存在
if [ ! -d "$BENCHMARK_DATA_DIR" ]; then
  echo "ERROR: $BENCHMARK_DATA_DIR not found. Run 03_pool_and_subsample.py first."
  exit 1
fi

#從barcodes.tsv/features.tsv的行數，便宜地估出n_cells/n_genes，
#不用真的載入矩陣——失敗時要在CSV裡寫這兩個欄位，但這時候工具本身
#已經當掉了，沒辦法從工具內部拿到這兩個數字，只好從檔案行數推回來
get_n_cells() { wc -l < "$1/barcodes.tsv" 2>/dev/null | tr -d ' ' || echo 0; }
get_n_genes() { wc -l < "$1/features.tsv" 2>/dev/null | tr -d ' ' || echo 0; }

#寫一筆失敗的紀錄進CSV，欄位順序跟三個benchmark腳本自己寫的OK紀錄一致
#(tool,label,n_cells,n_genes,elapsed_sec,peak_memory_MB,status)，
#elapsed_sec、peak_memory_MB留空，status寫FAILED
write_failed_row() {
  local tool="$1" label="$2" n_cells="$3" n_genes="$4" csv="$5"
  if [ ! -f "$csv" ]; then
    echo "tool,label,n_cells,n_genes,elapsed_sec,peak_memory_MB,status" > "$csv"
  fi
  echo "${tool},${label},${n_cells},${n_genes},,,FAILED" >> "$csv"
}

#用 timeout 包住一次工具呼叫(-k 60:逾時後先送SIGTERM，60秒內還沒停
#就強制SIGKILL)，逾時或執行失敗都寫一筆FAILED紀錄，不讓腳本卡死或中斷。
#exit code 124是GNU timeout的慣例，代表「逾時被砍」，訊息裡特別標出來，
#方便之後在終端機紀錄裡分辨是「真的卡住」還是「單純crash」
run_tool() {
  local tool_name="$1" label="$2" n_cells="$3" n_genes="$4"
  shift 4
  #注意:不能寫成 "if ! timeout ...; then local rc=$?"——"!" 本身會翻轉
  #$?，所以 then 區塊裡的 $? 其實是被"!"改過的值，不是timeout原本的
  #結束碼。要用 if/else(不加"!")才能拿到timeout真正的結束碼
  local rc
  if timeout -k 60 "$TOOL_TIMEOUT_SEC" "$@"; then
    return 0
  else
    rc=$?
  fi
  if [ "$rc" -eq 124 ]; then
    echo "$tool_name run TIMED OUT (>${TOOL_TIMEOUT_SEC}s) for $label"
  else
    echo "$tool_name run FAILED (exit $rc) for $label"
  fi
  write_failed_row "$tool_name" "$label" "$n_cells" "$n_genes" "$RESULTS_CSV"
}

#自動判斷這次要跑第幾輪:找「還沒有對應CSV檔案」的第一個輪數
round=""
for i in $(seq 1 "$N_ROUNDS"); do
  candidate="$RESULTS_DIR/benchmark_results_round${i}.csv"
  if [ ! -f "$candidate" ]; then
    round=$i
    break
  fi
done

if [ -z "$round" ]; then
  echo "All $N_ROUNDS rounds already have results:"
  for i in $(seq 1 "$N_ROUNDS"); do
    echo "  $RESULTS_DIR/benchmark_results_round${i}.csv"
  done
  echo "Run average_results.py to combine them, or delete one of the files above to re-run that round."
  exit 0
fi

RESULTS_CSV="$RESULTS_DIR/benchmark_results_round${round}.csv"
echo ""
echo "##################################################"
echo "ROUND $round / $N_ROUNDS -- writing to $RESULTS_CSV"
echo "##################################################"

#依數字大小(由小到大)列出所有子樣本資料夾名稱，取代直接glob "n_*"
#(那樣是字母順序，n_100000會排在n_20000前面，見上方註解)
mapfile -t subset_labels < <(cd "$BENCHMARK_DATA_DIR" && ls -d n_* | sort -t_ -k2 -n)

#主迴圈:對每個子樣本大小(由小到大)，依序跑三個工具
for label in "${subset_labels[@]}"; do
  subset_dir="$BENCHMARK_DATA_DIR/$label"
  n_cells=$(get_n_cells "$subset_dir")
  n_genes=$(get_n_genes "$subset_dir")
  echo ""
  echo "=================================================="
  echo "Round $round -- Subset: $label"
  echo "=================================================="

  echo "--- Scanpy ---"
  run_tool "Scanpy" "$label" "$n_cells" "$n_genes" \
    python "$SCRIPT_DIR/run_scanpy_benchmark.py" "$subset_dir" "$label" "$RESULTS_CSV"

  echo "--- Seurat ---"
  run_tool "Seurat" "$label" "$n_cells" "$n_genes" \
    Rscript "$SCRIPT_DIR/run_seurat_benchmark.R" "$subset_dir" "$label" "$RESULTS_CSV"

  echo "--- BPCells ---"
  run_tool "BPCells" "$label" "$n_cells" "$n_genes" \
    Rscript "$SCRIPT_DIR/run_bpcells_benchmark.R" "$subset_dir" "$label" "$RESULTS_CSV"
done

echo ""
echo "Round $round / $N_ROUNDS complete. Results in: $RESULTS_CSV"
if [ "$round" -lt "$N_ROUNDS" ]; then
  echo "Run 'bash run_all.sh' again to run round $((round + 1))."
else
  echo "All $N_ROUNDS rounds done. Run average_results.py to combine them, then summarize_results.py to plot."
fi
