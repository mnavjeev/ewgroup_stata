#!/usr/bin/env python3
"""Build and test isolated sibling checkouts; optionally benchmark exact fits."""
import argparse
import csv
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time


def run(command, cwd, logfile, env, timed=False):
    if timed:
        # A fresh helper gives each process its own getrusage high-water mark.
        # Unlike macOS `time -l`, this does not need restricted sysctl access.
        meter = ("import resource,subprocess,sys; "
                 "p=subprocess.run(sys.argv[1:]); "
                 "rss=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss; "
                 "print('__PEAK_RSS_BYTES__='+str(int(rss)*(1 if sys.platform=='darwin' else 1024))); "
                 "sys.exit(p.returncode)")
        command = [sys.executable, "-c", meter] + command
    started = time.monotonic()
    with logfile.open("w") as out:
        result = subprocess.run(command, cwd=cwd, env=env, stdout=out, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"Command failed ({result.returncode}); see {logfile}\n" +
                           logfile.read_text(errors="replace")[-5000:])
    return time.monotonic() - started


def stata_run(executable, script, cwd, work, env, success, timed=False):
    log = work / (script.stem + "-process.log")
    stata_log = cwd / (script.stem + ".log")
    if stata_log.exists():
        stata_log.unlink()
    elapsed = run([executable, "-b", "do", str(script)], cwd, log, env, timed)
    contents = stata_log.read_text(errors="replace") if stata_log.exists() else ""
    # Stata batch can exit zero after an error; require a completion marker.
    if success not in contents or re.search(r"^r\([0-9]+\);", contents, re.M):
        body = contents[contents.find(". do "):] if ". do " in contents else contents
        raise RuntimeError(f"Stata check failed; see {stata_log}\n" + body[-6000:])
    return elapsed, log, contents


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stata", default=os.environ.get("STATA_BIN"))
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--benchmarks", action="store_true",
                        help="also run exact 10k/50k scalar and 10k full 3D fits")
    parser.add_argument("--benchmark-only", action="store_true")
    parser.add_argument("--library", type=Path, help="reuse an installed R library for benchmarks")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    work = args.work_dir or Path(tempfile.mkdtemp(prefix="ewgroup-checks-"))
    work = work.resolve(); work.mkdir(parents=True, exist_ok=True)
    print(f"Logs and isolated checkouts: {work}", flush=True)
    stata = args.stata
    if not stata:
        stata = next((shutil.which(s) for s in ("stata-mp", "stata-se", "stata")
                      if shutil.which(s)), None)
    if not stata and Path("/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp").exists():
        stata = "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp"
    if not stata:
        parser.error("Specify a licensed Stata executable with --stata or STATA_BIN.")
    r = shutil.which("R"); rscript = shutil.which("Rscript")
    if not r or not rscript:
        parser.error("R and Rscript must be on PATH.")
    env = dict(os.environ, LC_ALL="C")
    for package in ("ewgroup", "ewgroup_stata"):
        shutil.copytree(root / package, work / package, dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns(".git", "*.o", "*.so", "*.dll", "*.log", "*.Rcheck"))
    lib = args.library.resolve() if args.library else work / "library"
    lib.mkdir(exist_ok=True)
    if not args.library:
        print("Installing R package", flush=True)
        run([r, "CMD", "INSTALL", "--no-multiarch", "--library=" + str(lib), "ewgroup"],
            work, work / "r-install.log", env)
    env["R_LIBS"] = str(lib) + os.pathsep + env.get("R_LIBS", "")

    if not args.benchmark_only:
        print("Building and checking R source package", flush=True)
        run([r, "CMD", "build", "ewgroup"], work, work / "r-build.log", env)
        archive = max(work.glob("ewgroup_*.tar.gz"), key=lambda p: p.stat().st_mtime)
        run([r, "CMD", "check", "--no-manual", "--no-build-vignettes", archive.name],
            work, work / "r-check.log", env)
        print("Checking compiled R, fallback R, and paper-reference fixtures", flush=True)
        run([rscript, "tests/make_reference.R", "--library=" + str(lib)],
            work / "ewgroup_stata", work / "fixtures.log", env)
        for name, marker in (("certify", "ewgroup certification checks passed"),
                             ("regressions", "ewgroup regression checks passed"),
                             ("check_net_install", "local net install check passed")):
            print(f"Stata: {name}", flush=True)
            stata_run(stata, work / "ewgroup_stata" / "tests" / (name + ".do"),
                      work / "ewgroup_stata", work, env, marker)
        print("R and Stata correctness/installation checks passed", flush=True)

    if args.benchmarks or args.benchmark_only:
        results = []
        for n, d in ((10000, 1), (50000, 1), (10000, 3)):
            tag = f"benchmark_{n}_{d}"
            script = work / (tag + ".R")
            script.write_text(f'''library(ewgroup)
set.seed(907)
n <- {n}; d <- {d}
B <- matrix(rnorm(n*d), n, d)
S <- lapply(seq_len(n), function(j) {{
  if (d == 1) matrix(runif(1, .2, 1)) else {{
    A <- matrix(rnorm(d*d), d); crossprod(A)/d + diag(.2, d)
  }}
}})
data <- as.data.frame(B); names(data) <- paste0("b", seq_len(d))
for (a in seq_len(d)) for (b in seq_len(d)) {{
  data[[paste0("s", a, b)]] <- vapply(S, function(s) s[a,b], numeric(1))
}}
elapsed <- system.time(fit <- ewgroup(B, S, sigma2=.15))[["elapsed"]]
stopifnot(all(is.finite(fit$theta)), is.finite(fit$alpha))
for (k in seq_len(d)) data[[paste0("expected",k)]] <- fit$theta[,k]
write.csv(data, "{tag}.csv", row.names=FALSE)
cat("ESTIMATOR_SECONDS=", elapsed, "\\n", sep="")
cat("ALPHA=", fit$alpha, "\\n", sep="")
''')
            print(f"Benchmark R: J={n}, d={d}", flush=True)
            rlog = work / (tag + "-R.log")
            run([rscript, str(script)], work, rlog, env, timed=True)
            bvars = " ".join(f"b{k}" for k in range(1, d + 1))
            svars = " ".join(f"s{a}{b}" for a in range(1, d + 1) for b in range(1, d + 1))
            doscript = work / (tag + ".do")
            doscript.write_text(f'''clear all
set more off
set type double
adopath ++ "{work / 'ewgroup_stata'}"
import delimited using "{tag}.csv", clear
timer clear
timer on 1
quietly ewgroup {bvars}, sigma({svars}) sigma2(.15) prefix(theta_) nomatrices
timer off 1
timer list 1
display "ESTIMATOR_SECONDS=" r(t1)
forvalues k=1/{d} {{
    assert abs(theta_`k' - expected`k') < 1e-8
}}
assert e(sample) == 1
display "BENCHMARK_PASSED"
''')
            print(f"Benchmark Stata: J={n}, d={d}", flush=True)
            _, slog, contents = stata_run(stata, doscript, work, work, env, "BENCHMARK_PASSED", timed=True)
            for language, logfile, timingtext in (("R", rlog, rlog.read_text()),
                                                  ("Stata", slog, contents)):
                text = logfile.read_text(errors="replace")
                seconds = re.findall(r"^ESTIMATOR_SECONDS=\s*([\d.eE+-]+)", timingtext, re.M)
                peak = re.search(r"__PEAK_RSS_BYTES__=(\d+)", text)
                results.append(dict(language=language, cells=n, dimension=d,
                                    seconds=seconds[-1] if seconds else "",
                                    peak_mib=round(int(peak.group(1))/1048576, 2) if peak else ""))
            with (work / "benchmarks.csv").open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=results[0].keys())
                writer.writeheader(); writer.writerows(results)
            print(results[-2:], flush=True)
    print(f"Completed. Results: {work}", flush=True)


if __name__ == "__main__":
    main()
