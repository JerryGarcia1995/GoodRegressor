import os

SubDirPrefix = "005"
StartCh = 21
EndCh = 60
RegressionSummaryFile = "output.txt"
TargetY = "Fracture"

BASE_DIR = os.getcwd()   # run this from the parent directory

def parse_hdat_file(hdat_path):
    """
    Returns:
      all_chain_files : list[str]
      recommended_file : str
    """
    all_chain_files = []
    recommended_file = None

    with open(hdat_path, "r") as f:
        lines = f.readlines()

    # 1) collect filenames from the main table (4th column)
    for line in lines:
        if line.strip().startswith(TargetY):
            parts = line.split()
            if len(parts) >= 4:
                all_chain_files.append(parts[3])

    # 2) find RECOMMEND section
    for i, line in enumerate(lines):
        if "==========RECOMMEND==========" in line:
            # "next next line"
            if i + 2 < len(lines):
                parts = lines[i + 2].split()
                if len(parts) >= 4:
                    recommended_file = parts[3]
            break

    if recommended_file is None:
        print(f"Warning: No RECOMMEND section in {hdat_path}")
        return all_chain_files, None

    return all_chain_files, recommended_file

for d in sorted(os.listdir(BASE_DIR)):
    dir_path = os.path.join(BASE_DIR, d)
    
    if not os.path.isdir(dir_path):
        continue

    if not d.startswith(SubDirPrefix):
        continue
    
    #try:
        #ch = int(d.split("_ch")[-1])
    #except ValueError:
        #continue

    #if not (StartCh <= ch <= EndCh):
        #continue

    try:
        hdat_path = os.path.join(dir_path, RegressionSummaryFile)
        basic_path = os.path.join(dir_path, "basic_properties_GoodRegressor.txt")

        if not (os.path.isfile(hdat_path) and os.path.isfile(basic_path)):
            print(f"Skipping {d}: missing files")
            continue

        all_chain_files, recommended_file = parse_hdat_file(hdat_path)
        if recommended_file is None:
            print(f"Skipping {d}: no recommended file")
            continue
            
        if recommended_file not in all_chain_files:
            print(f"Skipping {d}: recommended file not in list")
            continue

        idx = all_chain_files.index(recommended_file)
        selected_files = all_chain_files[:idx + 1]

        out_path = os.path.join(BASE_DIR, f"{d}.txt")

        with open(out_path, "w") as out:
            out.write("name " + TargetY + "\n")
            out.write("beta_target_type(e.g.,none-logit-probit-cloglog-cauchit-nloglog) none\n")
            out.write(f"basic_properties_GoodRegressor_filename {os.path.abspath(basic_path)}\n")
            out.write("chain_filenames\n")

            for fname in selected_files:
                out.write(os.path.abspath(os.path.join(dir_path, fname)) + "\n")

        print(f"Written: {out_path}")

    except Exception as e:
        print(f"Error in {d}: {e}")
        continue