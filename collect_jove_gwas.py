import argparse
import csv
import gzip
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

OUTCOMES = dict(APPENDACUT="K11_APPENDACUT", APPENDACUT_COMPLIC="APPENDACUT_COMPLIC", APPENDOTH="K11_APPENDOTH")
EXPOSURES = {
    "saturated": "saturated_fatty_acids.gz",
    "Polyunsaturated": "polyunsaturated_fatty_acids.gz",
    "Monounsaturated": "monounsaturated_fatty_acids.gz",
    "OMEGA6": "omega6_fatty_acids.gz",
    "OMEGA3": "omega3_fatty_acids.gz",
    "Monoradio": "monounsaturated_fatty_acid_ratio.gz",
    "Polyradio": "polyunsaturated_fatty_acid_ratio.gz",
    "OMEGAradio": "omega6_to_omega3_ratio.gz",
}
FIELDS = "SNP chrom pos ea oa beta se p eaf source".split()

def save(path, rows):
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(FIELDS)
        w.writerows(rows)

def scan_finngen(task):
    loc, dest, short, endpoint, target = task
    hits, lookup = [], []
    with gzip.open(loc / f"finngen_R11_{endpoint}.gz", "rt") as f:
        assert next(f).split()[:11] == "#chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt".split()
        for line in f:
            a = line.split("\t")
            p = float(a[6])
            rsids = a[4].split(",")
            gws = p < 5e-8 and a[0].isdigit() and 1 <= int(a[0]) <= 22
            for snp in rsids:
                if not snp.startswith("rs") or not (gws or snp in target):
                    continue
                rec = (snp, a[0], a[1], a[3], a[2], a[8], a[9], a[6], a[10], endpoint)
                if snp in target:
                    lookup.append(rec)
                if gws and a[3] in ("A","C","G","T") and a[2] in ("A","C","G","T"):
                    hits.append(rec)
    save(dest / f"{short}_gws.tsv", hits)
    save(dest / f"{short}_forward.tsv", lookup)
    print(endpoint, "genome-wide significant SNP records", len(hits), "forward records", len(lookup), flush=True)
    return {h[0] for h in hits}

def scan_exposure(task):
    loc, dest, key, filename, target = task
    found = []
    with gzip.open(loc / filename, "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                sid = line.rstrip().split("\t")[9]
                break
        for line in f:
            a = line.split("\t", 4)
            if a[2] not in target:
                continue
            tail = a[4].rstrip().split("\t")
            fmt = dict(zip(tail[4].split(":"), tail[5].split(":")))
            try:
                b, se, lp, af = [float(fmt[n]) for n in ("ES","SE","LP","AF")]
            except (KeyError, ValueError):
                continue
            if all(map(math.isfinite, (b,se,lp,af))) and se > 0 and 0 < af < 1:
                found.append((a[2],a[0],a[1],tail[0],a[3],b,se,10**(-min(lp,300)),af,sid))
    assert found, key
    save(dest / f"{key}_lookup.tsv", found)
    print(key, sid, "lookup records", len(found), flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("gwas", type=Path, metavar="data_directory")
    p.add_argument("working",type=Path)
    p.add_argument("forward_snps",type=Path)
    a=p.parse_args()
    a.working.mkdir(parents=True,exist_ok=True)
    s=set(a.forward_snps.read_text().splitlines())
    with ProcessPoolExecutor(3) as pool:
        sets=list(pool.map(scan_finngen,[(a.gwas / "outcome",a.working,o,e,s) for o,e in OUTCOMES.items()]))
    all_s=s.union(*sets)
    with ProcessPoolExecutor(4) as pool:
        list(pool.map(scan_exposure,[(a.gwas / "exposure",a.working,k,f,all_s) for k,f in EXPOSURES.items()]))

if __name__=="__main__":
    main()
