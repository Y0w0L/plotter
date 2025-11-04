import subprocess
import datetime
from tqdm import tqdm
import os
import argparse

parser = argparse.ArgumentParser(description="Run analysis jobs.")
parser.add_argument(
    "-l",
    "--log",
    action="store_true",
    help="Enable writing a detailed log file."
)
args = parser.parse_args()
WRITE_LOG = args.log

# ======================================================================
# 1. パラメータ空間の定義
# ======================================================================
model_ = ["masetti"]
threshold_ = [0, 10, 25, 50, 100, 200, 300, 400]
voltage_ = [10, 7, 4]
pixel_pitch_ = [15, 22.5]
chip_type_ = ["std", "gap"]
beam_info_ = ["e3GeV"]

seed_threshold_ = [0]
neighbor_threshold_ = [0, 50, 100, 200, 300, 400, 500, 1000, 1500]

# ======================================================================
# 2. 実行する全ジョブのリストを生成
# ======================================================================
param_list = []
for model in model_:
    for beam_info in beam_info_:
        for chip_type in chip_type_:
            for voltage in voltage_:
                for pixel_pitch in pixel_pitch_:
                    # enumerateを使って、thresholdのインデックス(i)を取得
                    for i, threshold in enumerate(threshold_):
                        for seed_threshold in seed_threshold_:
                            for neighbor_threshold in neighbor_threshold_:
                                param_list.append(
                                    (model, beam_info, chip_type, voltage, pixel_pitch, i, threshold, seed_threshold, neighbor_threshold)
                                )

# ======================================================================
# 3. ログファイルとジョブの実行
# ======================================================================
# ログを保存するディレクトリを作成
# log_dir = "log"
# if not os.path.exists(log_dir):
#     os.makedirs(log_dir)

#timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
log_filename = None
if WRITE_LOG:
    log_dir = "log"
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = os.path.join(log_dir, f"analysis_log_{timestamp}.txt")

print(f"Starting all jobs... Total jobs: {len(param_list)}")

try:
    log_file_context = open(log_filename, 'a', encoding='utf-8') if WRITE_LOG else open(os.devnull, 'w')

    with log_file_context as log_file:
        # tqdmを使って進捗バーを表示
        for params in tqdm(param_list, desc="Overall Progress"):
            # パラメータをアンパック
            model, beam_info, chip_type, voltage, pixel_pitch, i, threshold, seed_threshold, neighbor_threshold = params

            # --- ファイル名を生成 ---
            # pixel_pitchを文字列に変換
            pitch_name = "15" if pixel_pitch == 15 else "22p5"

            input_filename = f"/home/towa/alice3/data/ce65_sim_202505/n{voltage}v/ce65_p{pitch_name}_{chip_type}_Thd{threshold}e_{beam_info}_{model}.root"
            output_filename = f"/home/towa/alice3/hist/ce65driftTime/ce65driftTime_{pitch_name}_{chip_type}_{voltage}V_SeedThd{seed_threshold}e_NeighborThd{neighbor_threshold}e_{i}.root"

            if os.path.exists(output_filename):
                tqdm.write(f"Skipping, file exists: {os.path.basename(output_filename)}")
                if WRITE_LOG:
                    log_file.write(f"--- Skipping, file exists: {os.path.basename(output_filename)}")
                    log_file.flush()
                continue

            # --- 実行するコマンドをリストとして作成 ---
            command = [
                "python3", "AnalysisCarrier.py",
                f"-i={input_filename}",
                f"-o={output_filename}",
                f"-st={seed_threshold}",
                f"-nt={neighbor_threshold}",
                f"-p={pixel_pitch}"
            ]
            
            # --- ジョブヘッダーをログに書き込む ---
            if WRITE_LOG:
                header = f"""
========================================================================
>>> Starting Job: 
    Input: {os.path.basename(input_filename)}
    Output: {os.path.basename(output_filename)}
    Params: pitch={pixel_pitch}, seed={seed_threshold}, neighbor={neighbor_threshold}
========================================================================
"""
                log_file.write(header)
                log_file.flush() # バッファを即座にファイルに書き込む

            stdout_dest = log_file if WRITE_LOG else subprocess.DEVNULL
            stderr_dest = log_file if WRITE_LOG else subprocess.DEVNULL

            # --- サブプロセスとして解析スクリプトを実行 ---
            # stdoutとstderrを両方ともログファイルにリダイレクトする
            subprocess.run(
                command,
                check=True,       # エラーが発生した場合に例外を発生させる
                stdout=stdout_dest,
                stderr=stderr_dest
            )

except Exception as e:
    print(f"\nAn error occurred: {e}")
    print(f"Please check '{log_filename}' for details.")
finally:
    print(f"\nAll jobs finished. You can check the detailed logs in '{log_filename}'.")