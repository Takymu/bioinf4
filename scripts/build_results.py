#!/usr/bin/env python3
"""Create one-row-per-input-variant results from cached source responses."""
import collections
import csv
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "raw"
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)


def load(name):
    return json.loads((RAW / (name + ".json")).read_text())


def dump_table(name, rows):
    assert rows, name
    with (OUT / name).open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def text(value):
    if isinstance(value, list):
        return "; ".join(str(x) for x in value)
    return "" if value is None else value


def priority(tc):
    return (not bool(tc.get("mane_select")), not bool(tc.get("canonical")), tc.get("biotype") != "protein_coding", tc["transcript_id"])


def main():
    lines = (ROOT / "input/student_02.vcf").read_text().splitlines()
    inputs = [x.split("\t") for x in lines if not x.startswith("#")]
    assert len(inputs) == 50 and len({r[2] for r in inputs}) == 50
    assert all(len(r) == 8 and len(r[3]) == len(r[4]) == 1 and r[3] in "ACGT" and r[4] in "ACGT" and r[3] != r[4] for r in inputs)
    vep = {x["id"]: x for x in load("vep")}
    assert set(vep) == {r[2] for r in inputs}
    clin = load("clinvar_summary")["result"]
    gnom = load("gnomad")
    output, transcripts, clin_output, gn_output, excluded = [], [], [], [], []
    matched_uids = set()
    ts = {"AG", "GA", "CT", "TC"}
    af_checks = 0
    for chrom, pos, ident, ref, alt, qual, filt, info in inputs:
        variant = "-".join([chrom, pos, ref, alt])
        vx = vep[ident]
        assert vx["assembly_name"] == "GRCh38"
        assert str(vx["start"]) == pos and str(vx["seq_region_name"]) == chrom
        assert vx["allele_string"] == ref + "/" + alt, (ident, vx["allele_string"])
        tc = [t for t in vx["transcript_consequences"] if t["variant_allele"] == alt]
        severe = vx["most_severe_consequence"]
        candidates = [t for t in tc if severe in t["consequence_terms"]]
        assert candidates, (ident, severe)
        selected = min(candidates, key=priority)
        mane_candidates = [t for t in tc if t.get("mane_select") and t["gene_id"] == selected["gene_id"]]
        mane = min(mane_candidates, key=priority) if mane_candidates else {}
        for t in tc:
            transcripts.append({"original_id": ident, "variant_grch38": variant, "gene": t.get("gene_symbol", t["gene_id"]), "gene_id": t["gene_id"], "transcript_id": t["transcript_id"], "consequence_terms": text(t["consequence_terms"]), "impact": t["impact"], "biotype": t["biotype"], "hgvsc": t.get("hgvsc", ""), "hgvsp": t.get("hgvsp", ""), "mane_select": t.get("mane_select", ""), "canonical": t.get("canonical", 0), "flags": text(t.get("flags", [])), "sift_prediction": t.get("sift_prediction", ""), "sift_score": t.get("sift_score", ""), "polyphen_prediction": t.get("polyphen_prediction", ""), "polyphen_score": t.get("polyphen_score", ""), "selected": int(t is selected)})

        accession = {"21": "NC_000021.9", "22": "NC_000022.11"}[chrom]
        target_spdi = f"{accession}:{int(pos)-1}:{ref}:{alt}"
        exact, overlapping = [], []
        for uid in clin["uids"]:
            record = clin[uid]
            for variation in record["variation_set"]:
                locations = [loc for loc in variation.get("variation_loc", []) if loc.get("assembly_name") == "GRCh38" and loc.get("chr") == chrom and loc.get("start") == pos and loc.get("stop") == pos]
                if not locations:
                    continue
                if variation.get("canonical_spdi") == target_spdi:
                    exact.append((uid, record, variation))
                    matched_uids.add(uid)
                else:
                    overlapping.append(uid)
                    excluded.append({"original_id": ident, "input_spdi": target_spdi, "clinvar_variation_id": uid, "clinvar_spdi": variation.get("canonical_spdi", ""), "title": record["title"], "reason": "different REF/ALT at the same GRCh38 position"})
        classifications = [r.get("germline_classification", {}).get("description", "") or "not provided" for _, r, _ in exact]
        conditions = sorted({trait["trait_name"] for _, r, _ in exact for trait in r.get("germline_classification", {}).get("trait_set", [])})
        evaluated = [r.get("germline_classification", {}).get("last_evaluated", "") for _, r, _ in exact]
        evaluated = ["" if value.startswith("1/01/01") else value for value in evaluated]
        clinical = {"original_id": ident, "variant_grch38": variant, "clinvar_present": "yes" if exact else "no", "clinvar_interpretation": text(classifications) if exact else "no exact-allele record found", "clinvar_variation_id": text([uid for uid, _, _ in exact]), "clinvar_vcv": text([r.get("accession", "") for _, r, _ in exact]), "clinvar_title": text([r["title"] for _, r, _ in exact]), "clinvar_review_status": text([r.get("germline_classification", {}).get("review_status", "") for _, r, _ in exact]), "clinvar_last_evaluated": text(evaluated), "clinvar_conditions": text(conditions), "same_position_other_allele_ids": text(overlapping), "match_rule": "GRCh38 chromosome and position plus exact canonical SPDI REF/ALT", "clinvar_url": text(["https://www.ncbi.nlm.nih.gov/clinvar/variation/" + uid + "/" for uid, _, _ in exact])}
        clin_output.append(clinical)

        response = gnom[variant]
        assert not response.get("errors"), (ident, response.get("errors"))
        gv = response["data"]["variant"]
        gn = {"original_id": ident, "variant_grch38": variant, "gnomad_dataset": "gnomad_r4", "gnomad_status": "found" if gv else "variant not found", "gnomad_url": "https://gnomad.broadinstitute.org/variant/" + variant + "?dataset=gnomad_r4"}
        if gv:
            assert gv["variant_id"] == variant
        for callset in ["exome", "genome"]:
            call = gv.get(callset) if gv else None
            for key in ["ac", "an", "af"]:
                gn["gnomad_" + callset + "_" + key] = call.get(key) if call else None
            gn["gnomad_" + callset + "_filters"] = text(call["filters"]) if call else "not available"
            if call and call["an"]:
                assert 0 <= call["ac"] <= call["an"]
                assert abs(call["af"] - call["ac"] / call["an"]) < 1e-12
                af_checks += 1
        gn_output.append(gn)

        coloc = [c for c in vx.get("colocated_variants", []) if c.get("id", "").startswith("rs") and c.get("allele_string", "").split("/")[0] == ref and alt in c.get("allele_string", "").split("/")[1:]]
        rsids = sorted({c["id"] for c in coloc})
        assert ident in rsids, (ident, rsids)
        result = {"original_id": ident, "chrom": chrom, "pos_grch38": int(pos), "ref": ref, "alt": alt, "variant_grch38": variant, "gene": selected.get("gene_symbol", selected["gene_id"]), "gene_id": selected["gene_id"], "most_severe_consequence": severe, "selected_transcript": selected["transcript_id"], "selected_biotype": selected["biotype"], "selected_consequences": text(selected["consequence_terms"]), "impact": selected["impact"], "hgvsc": selected.get("hgvsc", ""), "hgvsp": selected.get("hgvsp", ""), "selected_mane_select": selected.get("mane_select", ""), "selected_canonical": selected.get("canonical", 0), "selected_flags": text(selected.get("flags", [])), "rsid": text(rsids), "mane_transcript": mane.get("transcript_id", ""), "mane_consequence": text(mane.get("consequence_terms", [])), "mane_hgvsc": mane.get("hgvsc", ""), "mane_hgvsp": mane.get("hgvsp", ""), "mane_refseq": mane.get("mane_select", ""), "severity_differs_from_mane": "yes" if mane and severe not in mane["consequence_terms"] else "no" if mane else "MANE not available", "transcript_count": len(tc), "substitution_type": "transition" if ref + alt in ts else "transversion", "selection_rule": "VEP most_severe_consequence, then MANE Select, canonical, protein_coding, transcript ID lexicographically", "ensembl_url": "https://www.ensembl.org/Homo_sapiens/Variation/Explore?v=" + ident}
        result.update({k: v for k, v in gn.items() if k not in ["original_id", "variant_grch38"]})
        result.update({k: v for k, v in clinical.items() if k not in ["original_id", "variant_grch38", "match_rule"]})
        output.append(result)
    assert len(output) == 50 and len({r["variant_grch38"] for r in output}) == 50
    assert sum(t["selected"] for t in transcripts) == 50
    (OUT / "annotated_variants.json").write_text(json.dumps(output, ensure_ascii=False, indent=2))
    dump_table("annotated_variants.tsv", output)
    dump_table("vep_all_transcripts.tsv", transcripts)
    dump_table("clinvar_annotation.tsv", clin_output)
    dump_table("gnomad_annotation.tsv", gn_output)
    if (RAW / "clinvar_selected_vcv.xml").exists():
        xml_root = ET.parse(RAW / "clinvar_selected_vcv.xml").getroot()
        submissions = []
        for archive in xml_root:
            for assertion in archive.findall(".//ClinicalAssertion"):
                acc = assertion.find("ClinVarAccession")
                submissions.append({"variation_id": archive.get("VariationID"), "vcv": archive.get("Accession") + "." + archive.get("Version"), "record_last_updated": archive.get("DateLastUpdated"), "scv": acc.get("Accession") + "." + acc.get("Version"), "submitter": acc.get("SubmitterName"), "classification": assertion.findtext(".//GermlineClassification"), "review_status": assertion.findtext(".//ReviewStatus")})
        (OUT / "clinvar_selected_submissions.json").write_text(json.dumps(submissions, ensure_ascii=False, indent=2))
        dump_table("clinvar_selected_submissions.tsv", submissions)
    if excluded:
        dump_table("clinvar_excluded_other_alleles.tsv", excluded)
    summary = {
        "input_file": "input/student_02.vcf", "reference": "GRCh38", "ensembl_release": load("ensembl_release"), "ensembl_rest_software": load("ensembl_software"),
        "input_sha256": hashlib.sha256((ROOT / "input/student_02.vcf").read_bytes()).hexdigest(),
        "variant_count": len(output), "by_chromosome": dict(collections.Counter(r[0] for r in inputs)),
        "substitution_counts": dict(collections.Counter(r[3] + ">" + r[4] for r in inputs)),
        "transitions": sum(r[3] + r[4] in ts for r in inputs), "transversions": sum(r[3] + r[4] not in ts for r in inputs),
        "ti_tv_ratio": sum(r[3] + r[4] in ts for r in inputs) / sum(r[3] + r[4] not in ts for r in inputs),
        "most_severe_consequence_counts": dict(collections.Counter(r["most_severe_consequence"] for r in output)),
        "clinvar_exact_variant_count": sum(r["clinvar_present"] == "yes" for r in output),
        "clinvar_classification_counts": dict(collections.Counter(r["clinvar_interpretation"] for r in output)),
        "clinvar_candidate_records": len(clin["uids"]), "clinvar_matched_records": len(matched_uids), "clinvar_excluded_other_allele_records": len(excluded),
        "gnomad_found_count": sum(r["gnomad_status"] == "found" for r in output), "gnomad_ac_an_af_checks": af_checks,
        "total_transcript_annotations": len(transcripts), "severity_differs_from_mane_count": sum(r["severity_differs_from_mane"] == "yes" for r in output),
        "all_input_ref_alt_match_vep": True, "no_api_errors_in_final_gnomad": True,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2))
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
