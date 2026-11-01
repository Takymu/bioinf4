#!/usr/bin/env python3
"""Download exact-allele GRCh38 annotations; cache unmodified JSON responses.

Python standard library only. Run from any working directory.
"""
import concurrent.futures
import datetime as dt
import json
from pathlib import Path
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "raw"
RAW.mkdir(exist_ok=True)
ROWS = [line.split("\t") for line in (ROOT / "input/student_02.vcf").read_text().splitlines() if not line.startswith("#")]
LOG = []


def request(name, url, body=None, timeout=100):
    path = RAW / (name + ".json")
    if path.exists():
        return json.loads(path.read_text())
    headers = {"Accept": "application/json", "Content-Type": "application/json", "User-Agent": "BioinformaticsAssignment4/1.0"}
    payload = None if body is None else json.dumps(body).encode()
    for attempt in range(3):
        started = dt.datetime.now(dt.timezone.utc).isoformat()
        try:
            req = urllib.request.Request(url, data=payload, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout) as response:
                data = response.read()
                obj = json.loads(data)
                info = {"name": name, "url": url, "method": "GET" if body is None else "POST", "body": body, "retrieved_utc": started, "http_status": response.status, "headers": dict(response.headers)}
            path.write_bytes(data)
            (RAW / (name + ".request.json")).write_text(json.dumps(info, ensure_ascii=False, indent=2))
            LOG.append(info)
            print("saved", name, len(data), flush=True)
            return obj
        except Exception as exc:
            print("retry", name, attempt + 1, type(exc).__name__, str(exc)[:200], flush=True)
            if attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))


def fetch_vep():
    params = {"canonical": 1, "mane": 1, "hgvs": 1, "numbers": 1, "protein": 1, "variant_class": 1, "transcript_version": 1, "tsl": 1, "check_existing": 1}
    url = "https://rest.ensembl.org/vep/human/region?" + urllib.parse.urlencode(params)
    body = {"variants": [" ".join(r) for r in ROWS]}
    result = request("vep", url, body, timeout=180)
    assert len(result) == len(ROWS), (len(result), len(ROWS))
    request("ensembl_release", "https://rest.ensembl.org/info/data?content-type=application/json")
    request("ensembl_software", "https://rest.ensembl.org/info/software?content-type=application/json")
    request("ensembl_assembly", "https://rest.ensembl.org/info/assembly/homo_sapiens?content-type=application/json")


def fetch_clinvar():
    # Search both rs identifiers and exact GRCh38 positions. Match REF/ALT in esummary later.
    clauses = [r[2] for r in ROWS]
    for chrom in sorted({r[0] for r in ROWS}):
        positions = " OR ".join(r[1] + "[chrpos38]" for r in ROWS if r[0] == chrom)
        clauses.append("(" + chrom + "[Chromosome] AND (" + positions + "))")
    term = " OR ".join(clauses)
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode({"db": "clinvar", "term": term, "retmode": "json", "retmax": 10000})
    result = request("clinvar_search", url)
    sr = result["esearchresult"]
    assert int(sr["count"]) == len(sr["idlist"]), "Incomplete ClinVar search export"
    ids = sr["idlist"]
    if ids:
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + urllib.parse.urlencode({"db": "clinvar", "id": ",".join(ids), "retmode": "json"})
        request("clinvar_summary", url)
    else:
        (RAW / "clinvar_summary.json").write_text(json.dumps({"result": {"uids": []}}))
    xml_path = RAW / "clinvar_selected_vcv.xml"
    if not xml_path.exists():
        xml_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode({"db": "clinvar", "rettype": "vcv", "is_variationid": "", "id": "544554,3305447", "from_esearch": "true"})
        started = dt.datetime.now(dt.timezone.utc).isoformat()
        with urllib.request.urlopen(xml_url, timeout=100) as response:
            data = response.read()
            ET.fromstring(data)
            metadata = {"url": xml_url, "retrieved_utc": started, "method": "GET", "http_status": response.status, "headers": dict(response.headers)}
        xml_path.write_bytes(data)
        (RAW / "clinvar_selected_vcv.request.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2))


def fetch_gnomad_one(row):
    variant = "-".join([row[0], row[1], row[3], row[4]])
    query = """query Variant($variantId: String!) {
      variant(variantId: $variantId, dataset: gnomad_r4) {
        variant_id rsids
        exome { ac an af filters }
        genome { ac an af filters }
      }
    }"""
    result = request("gnomad/" + variant, "https://gnomad.broadinstitute.org/api", {"query": query, "variables": {"variantId": variant}})
    # Retain errors explicitly, never interpret a failed API call as zero AF.
    if result.get("errors"):
        print("gnomAD graphql errors", variant, result["errors"], flush=True)
    time.sleep(0.6)
    return variant, result


def fetch_gnomad():
    (RAW / "gnomad").mkdir(exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        result = dict(pool.map(fetch_gnomad_one, ROWS))
    (RAW / "gnomad.json").write_text(json.dumps(result, ensure_ascii=False, indent=2))


def main():
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        jobs = [pool.submit(fn) for fn in (fetch_vep, fetch_clinvar, fetch_gnomad)]
        for job in jobs:
            job.result()
    print("All annotations downloaded.")


if __name__ == "__main__":
    main()
