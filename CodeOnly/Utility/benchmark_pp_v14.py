#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Postprocess benchmark results (no pandas/numpy):

1) Parse model_logs/{Model}_fold{*}_{sbatch#}_{*}.log.txt
   log format expected:
     R2=0.7967 | RMSE=185.9119 | MAE=123.8197
   -> compute mean/std per model and rank by mean R2 (desc)

2) Merge out/*_predictions_{Model}_fold*.csv (5 folds per model)
   -> output wide TSV with columns:
      {Model}_Pred_<metric>    {Model}_True_<metric>   ...
   ordered by ranking
"""

import os, glob, re, csv, math, argparse

RE_MODEL_FOLD = re.compile(r"^(?P<model>.+?)_fold(?P<fold>\d+)_", re.IGNORECASE)

# log line example: R2=0.7967 | RMSE=185.9119 | MAE=123.8197
RE_METRICS = re.compile(
    r"R2\s*=\s*([-+0-9.eE]+)\s*\|\s*RMSE\s*=\s*([-+0-9.eE]+)\s*\|\s*MAE\s*=\s*([-+0-9.eE]+)"
)

def mean(xs):
    return sum(xs) / len(xs) if xs else float("nan")

def std(xs):
    # sample std (ddof=1)
    n = len(xs)
    if n <= 1:
        return 0.0
    mu = mean(xs)
    return math.sqrt(sum((x - mu) ** 2 for x in xs) / (n - 1))

def parse_log_metrics(path):
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        txt = f.read()
    ms = list(RE_METRICS.finditer(txt))
    if not ms:
        return None
    m = ms[-1]  # take last occurrence
    r2   = float(m.group(1))
    rmse = float(m.group(2))
    mae  = float(m.group(3))
    return rmse, mae, r2

def collect_metrics(model_logs_dir):
    paths = sorted(glob.glob(os.path.join(model_logs_dir, "*.log.txt")))
    data = {}  # model -> dict(metric -> list)
    for p in paths:
        base = os.path.basename(p)
        mm = RE_MODEL_FOLD.match(base)
        if not mm:
            continue
        model = mm.group("model")
        triple = parse_log_metrics(p)
        if triple is None:
            continue
        rmse, mae, r2 = triple
        data.setdefault(model, {"RMSE": [], "MAE": [], "R2": []})
        data[model]["RMSE"].append(rmse)
        data[model]["MAE"].append(mae)
        data[model]["R2"].append(r2)
    return data

def summarize(data):
    rows = []
    for model, d in data.items():
        rows.append({
            "Model": model,
            "R2_mean": mean(d["R2"]),
            "R2_std": std(d["R2"]),
            "RMSE_mean": mean(d["RMSE"]),
            "RMSE_std": std(d["RMSE"]),
            "MAE_mean": mean(d["MAE"]),
            "MAE_std": std(d["MAE"]),
            "n_folds_found": len(d["R2"]),
        })
    rows.sort(key=lambda r: (-r["R2_mean"], r["RMSE_mean"]))
    return rows

def save_metrics_tsv(rows, out_path):
    if not rows:
        return
    keys = list(rows[0].keys())
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

def print_table(rows):
    print("\nModel\t<R2>\t<RMSE>\t<MAE>\t(n_folds)\n")
    for r in rows:
        print(f"{r['Model']}\t{r['R2_mean']:.6f}\t{r['RMSE_mean']:.6f}\t{r['MAE_mean']:.6f}\t{r['n_folds_found']}")
    print("\n(Extra) mean ± std")
    for r in rows:
        print(
            f"{r['Model']}: "
            f"R2 {r['R2_mean']:.6f} ± {r['R2_std']:.6f}, "
            f"RMSE {r['RMSE_mean']:.6f} ± {r['RMSE_std']:.6f}, "
            f"MAE {r['MAE_mean']:.6f} ± {r['MAE_std']:.6f}"
        )

def find_fold_csvs(out_dir, model):
    paths = glob.glob(os.path.join(out_dir, f"*predictions_{model}_fold*.csv"))
    def fold_key(p):
        m = re.search(r"_fold(\d+)\.csv$", os.path.basename(p))
        return int(m.group(1)) if m else 999
    return sorted(paths, key=fold_key)

def read_two_cols_csv(path, metric):
    true_col = f"True_{metric}"
    pred_col = f"Pred_{metric}"
    true_vals, pred_vals = [], []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        r = csv.reader(f)
        header = next(r, None)
        if header is None:
            return [], []
        try:
            i_true = header.index(true_col)
            i_pred = header.index(pred_col)
        except ValueError:
            # fallback: first two columns
            i_true, i_pred = 0, 1
        for row in r:
            if not row:
                continue
            try:
                true_vals.append(float(row[i_true]))
                pred_vals.append(float(row[i_pred]))
            except Exception:
                continue
    return true_vals, pred_vals

def concat_model_folds(out_dir, model, metric):
    paths = find_fold_csvs(out_dir, model)
    if not paths:
        return None
    true_all, pred_all = [], []
    for p in paths:
        t, pr = read_two_cols_csv(p, metric)
        true_all.extend(t)
        pred_all.extend(pr)
    return true_all, pred_all

def build_wide_tsv(out_dir, ordered_models, metric, out_path):
    arrays = []
    maxlen = 0
    for model in ordered_models:
        res = concat_model_folds(out_dir, model, metric)
        if res is None:
            continue
        true_vals, pred_vals = res
        maxlen = max(maxlen, len(true_vals))
        arrays.append((model, true_vals, pred_vals))

    if not arrays:
        print("No fold prediction CSVs found.")
        return

    def pad(xs):
        if len(xs) < maxlen:
            return xs + [""] * (maxlen - len(xs))
        return xs

    header = []
    cols = []
    for model, true_vals, pred_vals in arrays:
        header.append(f"{model}_Pred_{metric}")
        header.append(f"{model}_True_{metric}")
        cols.append(pad([str(x) for x in pred_vals]))
        cols.append(pad([str(x) for x in true_vals]))

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for i in range(maxlen):
            w.writerow([col[i] for col in cols])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-logs", default="model_logs")
    ap.add_argument("--out-dir", default="out")
    ap.add_argument("--metric", default="Ea")
    ap.add_argument("--metrics-tsv", default="metrics_summary.tsv")
    ap.add_argument("--wide-tsv", default="predictions_wide.tsv")
    args = ap.parse_args()

    data = collect_metrics(args.model_logs)
    rows = summarize(data)

    if not rows:
        print("No metrics parsed. Check your log format includes: R2=... | RMSE=... | MAE=...")
        return

    print_table(rows)
    save_metrics_tsv(rows, args.metrics_tsv)
    print(f"\nSaved metrics summary: {args.metrics_tsv}")

    ordered_models = [r["Model"] for r in rows]
    build_wide_tsv(args.out_dir, ordered_models, args.metric, args.wide_tsv)
    print(f"Saved wide predictions TSV: {args.wide_tsv}")

if __name__ == "__main__":
    main()
