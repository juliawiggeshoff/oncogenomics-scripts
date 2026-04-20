#!/usr/bin/env python3
"""
Fetch [HUGO_SYMBOL] mutations + clinical data from cBioPortal, then enrich with
Genome Nexus annotations (HGVSg, HGVSc, Exon, gnomAD, ClinVar, dbSNP,
SIGNAL, and optionally OncoKB Annotation)

Requires session.json from a browser cBioPortal session
Generated like this:
curl -s "https://www.cbioportal.org/api/session/main_session/[session_id]" | python3 -m json.tool > session.json

session id information is retrieve after "session_id=" in the https

Usage:
    python3 fetch_gene_cbioportal.annotate_genomenexus.py --session session.json --gene [GENE_SYMBOL] --token [ONCOKB_TOKEN]

    ONCOKB_TOKEN is optional. If omitted, the Annotation column stays empty.
    You can also set the environment variable ONCOKB_TOKEN instead.

Outputs:
    [GENE_SYMBOL]_mutations_final.tsv              --> cBioPortal data only
    [GENE_SYMBOL]_mutations_final.tsv              --> + all Genome Nexus annotations
"""

import json, sys, os, time, csv, requests
from collections import defaultdict

BASE_CBP = "https://www.cbioportal.org/api"
BASE_GN  = "https://www.genomenexus.org"

## Target column information to match the browser-downloaded TSV
TARGET_COLS = [
    "Gene","Study of Origin","Sample ID","Cancer Type","Cancer Type Detailed",
    "Protein Change","Annotation","Custom Driver","Custom Driver Tiers",
    "Functional Impact","Mutation Type","Variant Type","Copy #","MS","VS",
    "Center","Chromosome","Start Pos","End Pos","Ref","Var","HGVSg","HGVSc",
    "Allele Freq (T)","Allele Freq (N)","Variant Reads","Ref Reads",
    "Variant Reads (Normal)","Ref Reads (Normal)","# Mut in Sample","Exon",
    "gnomAD","ClinVar","dbSNP","SIGNAL",
    "Adenocarcinoma subtype WHO2015","Adjuvant","Adjuvant Chemotherapy",
    "Adjuvant Immunotherapy","Adjuvant Targeted","Adjuvant Therapy",
    "Adjuvant Treatment","Adjuvant XRT","Adrenal Met Status (Months)",
    "Adrenal Met Status","Age","Age at Death","Age at Surgery/Biopsy",
    "Age at First Mets Dx","Age at Last Contact","Age at Resection",
    "Age at Sequencing","Age at Which Sequencing was Reported (Years)",
    "Age At Surgery","Patient Current Age","Age Greater than Median","Age (yrs)",
    "Neoplasm Disease Stage American Joint Committee on Cancer Code",
    "American Joint Committee on Cancer Publication Version Type","Albumin",
    "ALK driver","ALK protein change","Aneuploidy Score","ARID1A driver",
    "ARID1A protein change","Best Overall Response","Biopsy contains MIP",
    "Biopsy contains Solid","BM occurrence","Day difference from BM to primary dx",
    "Bone Met Status (Months)","Bone Met Status","BOR","BRAF driver",
    "BRAF protein change","Buffa Hypoxia Score","TCGA PanCanAtlas Cancer Type Acronym",
    "Cell Cycle","cfDNA Input (ng)","Chemotherapy","Tobacco History",
    "Clinically reported PD-L1 score",
    "Neoplasm American Joint Committee on Cancer Clinical Group Stage",
    "Clinical Trial","Clinical N Stage Code","CNS Met Status (Months)",
    "CNS Met Status","Clinical Nodal Status","Cohort","Clinical Stage",
    "CT scan type","CT Size","CT Slice Spacing","Cytology fixation type",
    "Last Communication Contact from Initial Pathologic Diagnosis Date",
    "Birth from Initial Pathologic Diagnosis Date",
    "Death from Initial Pathologic Diagnosis Date",
    "Last Alive Less Initial Pathologic Diagnosis Date Calculated Day Value",
    "Days to Last Followup","Durable clinical response","Death","Date of death",
    "Disease Free (Months)","Disease Free Status",
    "Largest diameter of resected BM lesion, mm","Disease",
    "Diffusion lung capacity for carbon monoxide",
    "Distant Mets: Adrenal Gland","Distant Mets: Biliary tract",
    "Distant Mets: Bladder/UT","Distant Mets: Bone","Distant Mets: Bowel",
    "Distant Mets: Breast","Distant Mets: CNS/Brain","Distant Mets: Distant LN",
    "Distant Mets: Female Genital","Distant Mets: Head and Neck",
    "Distant Mets: Intra-Abdominal","Distant Mets: Kidney","Distant Mets: Liver",
    "Distant Mets: Lung","Distant Mets: Male Genital","Distant Mets: Mediastinum",
    "Distant Mets: Ovary","Distant Mets: Pleura","Distant Mets: PNS",
    "Distant Mets: Skin","Distant Mets: Unspecified","dNLR","Dosage",
    "Driver Mutations","Date of last drug administration","Drug start date",
    "Months of disease-specific survival","Disease-specific Survival status",
    "Durable Clinical Benefit","ECOG","ECOG performance status",
    "ECOG Performance Status","EGFR driver","EGFR protein change","ERBB2 driver",
    "ERBB2 protein change","Ethnicity Category","Ever Met to Site: Adrenal",
    "Ever Met to Site: Bone","Ever Met to Site: CNS","Ever Met to Site: Liver",
    "Ever Met to Site: Lymph Nodes","Ever Met to Site: Lung",
    "Ever Met to Site: Pleura","Exome Sequencing Status","Extrapulmonary","FEV1",
    "FGA","FGA Facets","Date of Last Contact","Form completion date",
    "Fraction Genome Altered","At Least 2 Years Follow Up","Genetic Ancestry Label",
    "Gene Panel","Genome Doubled","Neoplasm Histologic Grade","Group Number",
    "Had Surgery","Halo tumor quality score","HIPPO Pathway","Histological grade",
    "Histology","Neoadjuvant Therapy Type Administered Prior To Resection Text",
    "Prior Cancer Diagnosis Occurence","HLA_A1 allele","HLA_A2 alleles",
    "HLA_B1 alleles","HLA_B2 alleles","HLA_C1 alleles","HLA_C2 alleles",
    "ICD-10 Classification",
    "International Classification of Diseases for Oncology, Third Edition ICD-O-3 Histology Code",
    "International Classification of Diseases for Oncology, Third Edition ICD-O-3 Site Code",
    "Immunotherapy","Samples used in Metastatic Lesion (ML) group",
    "Was IMPACT done on the same tissue that PD-L1 IHC was done on?",
    "Group Assignment for Primary Samples","IMSIG B-cells","IMSIG interferon",
    "IMSIG macrophages","IMSIG monocytes","IMSIG neutrophils","IMSIG NK cells",
    "IMSIG plasma cells","IMSIG proliferation","IMSIG translation","IMSIG T-cells",
    "Informed consent verified","Institute Source","Intracranial disease progression",
    "Intracranial disease progression Type",
    "Samples used in Matched Analysis (Group 5)","In PanCan Pathway Analysis",
    "IO drug name","Line of therapy","IRB","Metastatic patient","Is WGD",
    "JS PD-L1 score","Ki67 Percentage","Karnofsky Performance Scale",
    "Lines of therapy prior to BM resection","Lines of treatment",
    "Number Treatment Lines Prior To Receiving Impact Results",
    "Liver Met Status (Months)","Liver Met Status","Lymph Node Involvement",
    "LN Met Status (Months)","LN Met Status","Lung Met Status (Months)",
    "Lung Met Status","Manual tumor annotation","Margin Status",
    "Mean Target Coverage Normal","Mean Target Coverage Tumor",
    "Metabolic Tumor Volume","Metastatic Burden","Metastatic Site","Met Count",
    "MET driver","MET protein change","Met Site Count","MGMT Status",
    "Molecular Smoking Signature","Monotherapy vs. Combination",
    "Months from Matched Primary","MSI Comment","MSI Score","MSI MANTIS Score",
    "MSIsensor Score","MSI Type","Mutation Rate","MYC Pathway","M Stage",
    "Neoadjuvant","Neoadjuvant Chemotherapy","Neoadjuvant Immunotherapy",
    "Neoadjuvant Targeted","Neoadjuvant XRT","Neoantigen Burden",
    "New Neoplasm Event Post Initial Therapy Indicator",
    "Nonsynonymous mutation burden","NOTCH Pathway","NRF2 Pathway",
    "NSCLC SubType","Number of BMS at Diagnosis","N Stage","Oncotree Code",
    "Organ System","Overall survival","Overall Survival (Months)",
    "Overall Survival Status","Other Patient ID","Overall Patient Histology",
    "Overall Response","Pack-year history","Pathologic Stage",
    "American Joint Committee on Cancer Metastasis Stage Code",
    "Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code",
    "American Joint Committee on Cancer Tumor Stage Code","Patient Display Name",
    "PD-L1 expression (Percentage)","PDL1 Expression","PD-L1 Score (%)","PD-L1 tissue site",
    "Percent Necrosis","Person Neoplasm Cancer Status",
    "Positron emission tomography, tumor background ratio","PFS date",
    "Progress Free Survival (Months)","CNS PFS Status","PI3K Pathway",
    "Pleural Invasion","Pleura Met Status (Months)","Pleura Met Status","Ploidy",
    "Pathologic Nodal Status","Post Sample Chemotherapy","Post Sample Immunotherapy",
    "Post Sample Targeted","Post Sample Tx","Post Sample XRT ",
    "Predicted neoantigen burden","Predominant Histologic Subtype",
    "Predominant Histologic Subtype","Pre Sample Chemotherapy",
    "Pre Sample Immunotherapy","Pre Sample Targeted","Pre Sample Tx","Pre Sample XRT ",
    "Primary Lymph Node Presentation Assessment","Primary Tumor Site",
    "Prior Diagnosis","Prior PCI","Prior TKI at any time","Prior Treatment",
    "Prior Treatment Of Cytotoxic ChemoTransversion High erapy","Prior Wbrt",
    "Clinical trial IRB","Pathologic Stage","Purity","Race Category",
    "Radiation Therapy","Ragnum Hypoxia Score","Recurrent/metastatic disease",
    "RET driver","RET protein change","Relapse Free Status (Months)",
    "Relapse Free Status","RNA Sequencing Status","ROS1 driver","ROS1 protein change",
    "RTK- RAS Pathway","Sample Class","Sample Collection Source",
    "Sample Collection Time","Number of Samples Per Patient","Sample coverage",
    "Sample pre any Lung Therapy","Sample Type","Sample type id","Sequencing Type",
    "Sex","Site","Smoker","Person Cigarette Smoking History Pack Year Value",
    "Smoking History","Person Cigarette Smoking History Pack Year Value",
    "Smoking Status","Somatic Status","Stage","Stage At Diagnosis","Stage at Draw",
    "Status","STK11 driver","STK11 protein change","Subtype","Subtype Abbreviation",
    "Subtype Group","Successful ctDx Lung","Standardized uptake values ",
    "Symptom at BM diagnosis?","Systemic therapy prior to BM resection?",
    "Target Therapy","Tumor Break Load","TGF-Beta Pathway",
    "Tissue Prospective Collection Indicator","Tissue Retrospective Collection Indicator",
    "Tissue Source Site","Tissue Source Site Code","Tissue Specimen Type",
    "TKI Treatment","Tumor Mutation Burden","TMB (nonsynonymous)",
    "Total Exonic Mutation Burden","TP53 Pathway","Chemotherapy",
    "Treatment prior to BM resection","Treatment Schedule","Treatment Type",
    "Tumor Morphologic Appearance on CT","Tumor Purity","Tumor Size (mm)",
    "Tumor Size on CT","Tumor Stage","Tumor SUV Max","Tumor Disease Anatomic Site",
    "Tumor Type","Tumor Volume (cm3)","Treatment Setting","Types of symptoms",
    "T Stage","Ubiquitous Assay Panel","Vascular invasion","Is the patient deceased?",
    "Patient Weight","WGD","Winter Hypoxia Score","WNT Pathway",
]
TARGET_COLS_DEDUP = list(dict.fromkeys(TARGET_COLS))

MUT_RENAME = {
    'gene_hugoGeneSymbol':   'Gene',
    'studyId':               'Study of Origin',
    'sampleId':              'Sample ID',
    'cancerType':            'Cancer Type',
    'cancerTypeDetailed':    'Cancer Type Detailed',
    'proteinChange':         'Protein Change',
    'functionalImpactScore': 'Functional Impact',
    'mutationType':          'Mutation Type',
    'variantType':           'Variant Type',
    'mutationStatus':        'MS',
    'validationStatus':      'VS',
    'center':                'Center',
    'chr':                   'Chromosome',
    'startPosition':         'Start Pos',
    'endPosition':           'End Pos',
    'referenceAllele':       'Ref',
    'variantAllele':         'Var',
    'tumorAltCount':         'Variant Reads',
    'tumorRefCount':         'Ref Reads',
    'normalAltCount':        'Variant Reads (Normal)',
    'normalRefCount':        'Ref Reads (Normal)',
    'dbSnpRsId':             'dbSNP',
}

CLIN_RENAME = {
    'Diagnosis Age':                 'Age',
    'Nonsynonymous Mutation Burden': 'Nonsynonymous mutation burden',
    'Progression Free Status':       'CNS PFS Status',
    'TMB':                           'Tumor Mutation Burden',
    'Mutation Count':                None,
}

def safe_freq(alt, ref):
    try:
        a, r = int(alt), int(ref)
        total = a + r
        return round(a / total, 6) if total > 0 else ''
    except (TypeError, ValueError, ZeroDivisionError):
        return ''

def deep_get(d, *keys, default=''):
    for k in keys:
        if not isinstance(d, dict):
            return default
        d = d.get(k)
        if d is None:
            return default
    return d if d is not None else default

def cbp_get(endpoint, params=None, timeout=120):
    r = requests.get(BASE_CBP + endpoint, params=params,
                     headers={"Accept": "application/json"}, timeout=timeout)
    if r.status_code != 200 or not r.text.strip():
        return None
    return r.json()

def cbp_post(endpoint, payload, params=None, timeout=300):
    r = requests.post(BASE_CBP + endpoint, json=payload, params=params,
                      headers={"Content-Type": "application/json"}, timeout=timeout)
    if r.status_code != 200 or not r.text.strip():
        print(f"  POST ERROR {r.status_code} ({endpoint}): {r.text[:200]}")
        return None
    return r.json()

import argparse

parser = argparse.ArgumentParser(
    description="Fetch gene mutations + clinical data from cBioPortal and annotate with Genome Nexus."
)
parser.add_argument("--session", required=True, metavar="SESSION_JSON",
                    help="Path to cBioPortal session JSON file")
parser.add_argument("--gene", required=True, metavar="GENE_SYMBOL",
                    help="Hugo gene symbol to fetch (e.g. ROS1, EGFR, TP53)")
parser.add_argument("--token", metavar="ONCOKB_TOKEN",
                    default=os.environ.get("ONCOKB_TOKEN", ""),
                    help="OncoKB API token. Can also be set via ONCOKB_TOKEN env var.")
args = parser.parse_args()

session_file = args.session
gene_symbol  = args.gene.upper()
oncokb_token = args.token

print(f"Gene: {gene_symbol}")

## Resolve gene symbol to Entrez ID via cBioPortal gene endpoint
gene_resp = requests.get(
    f"{BASE_CBP}/genes/{gene_symbol}",
    headers={"Accept": "application/json"},
    timeout=30
)
if gene_resp.status_code != 200 or not gene_resp.text.strip():
    print(f"ERROR: could not resolve gene symbol '{gene_symbol}'")
    print(f"  Response: {gene_resp.status_code} {gene_resp.text[:200]}")
    sys.exit(1)
gene_info = gene_resp.json()
entrez_id = gene_info['entrezGeneId']
print(f"  Entrez ID: {entrez_id}  ({gene_info.get('hugoGeneSymbol', gene_symbol)})")

if oncokb_token:
    print(f"OncoKB token provided: Annotation column will be populated.")
else:
    print("No OncoKB token: Annotation column will be empty.")
    print( "Provide as: --token YOUR_TOKEN or set env ONCOKB_TOKEN=YOUR_TOKEN")

## load session.json 
with open(session_file) as f:
    session = json.load(f)

case_ids = session['data']['case_ids']
entries  = [e for e in case_ids.split('+') if ':' in e]

study_samples = defaultdict(list)
for entry in entries:
    study, sample = entry.split(':', 1)
    study_samples[study].append(sample)

all_study_ids = sorted(study_samples.keys())
total_samples = sum(len(v) for v in study_samples.values())
print(f"\nStudies: {len(all_study_ids)}  |  Samples: {total_samples}")

## 1. mutations
print(f"\n[1/4] Fetching {gene_symbol} mutations...")
sample_mol_ids = [
    {"molecularProfileId": f"{s}_mutations", "sampleId": sid}
    for s, sids in study_samples.items() for sid in sids
]
t0 = time.time()
mutations = cbp_post("/mutations/fetch",
                     {"entrezGeneIds": [entrez_id],
                      "sampleMolecularIdentifiers": sample_mol_ids},
                     params={"projection": "DETAILED"})
if not mutations:
    sys.exit(1)
print(f"  -> {len(mutations)} records  ({time.time()-t0:.1f}s)")

for m in mutations:
    gene = m.pop('gene', {}) or {}
    for k, v in gene.items():
        m[f'gene_{k}'] = v

## 2. sample clinical data
print("[2/4] Fetching sample-level clinical data...")
sample_clin = defaultdict(dict)
t0 = time.time()
for study in all_study_ids:
    data = cbp_get(f"/studies/{study}/clinical-data",
                   params={"clinicalDataType": "SAMPLE", "projection": "SUMMARY"})
    if not data:
        continue
    our = set(study_samples[study])
    for row in data:
        if row.get('sampleId') in our:
            sample_clin[(study, row['sampleId'])][row['clinicalAttributeId']] = row['value']
    print(f"  {study}")
    time.sleep(0.1)
print(f"  done ({time.time()-t0:.1f}s)")

## 3. patient clinical data
print("[3/4] Fetching patient-level clinical data...")
patient_clin = defaultdict(dict)
study_patients = defaultdict(set)
for m in mutations:
    if m.get('patientId') and m.get('studyId'):
        study_patients[m['studyId']].add(m['patientId'])

t0 = time.time()
for study in all_study_ids:
    our = study_patients.get(study, set())
    if not our:
        continue
    data = cbp_get(f"/studies/{study}/clinical-data",
                   params={"clinicalDataType": "PATIENT", "projection": "SUMMARY"})
    if not data:
        continue
    for row in data:
        if row.get('patientId') in our:
            patient_clin[(study, row['patientId'])][row['clinicalAttributeId']] = row['value']
    print(f"  {study}")
    time.sleep(0.1)
print(f"  done ({time.time()-t0:.1f}s)")

## 4. clinical attribute display names
print("[4/4] Fetching clinical attribute display names...")
attr_display = {}
for study in all_study_ids:
    data = cbp_get(f"/studies/{study}/clinical-attributes",
                   params={"projection": "SUMMARY"})
    if data:
        for attr in data:
            name = attr['displayName']
            corrected = CLIN_RENAME.get(name, name)
            if corrected is not None:
                attr_display[attr['clinicalAttributeId']] = corrected
    time.sleep(0.1)
print(f"  -> {len(attr_display)} attributes resolved")

## 5. build base rows
print("\nBuilding rows...")
rows = []
for m in mutations:
    sid  = m.get('sampleId', '')
    stid = m.get('studyId',  '')
    pid  = m.get('patientId','')

    row = {col: '' for col in TARGET_COLS_DEDUP}

    for api_key, val in m.items():
        col = MUT_RENAME.get(api_key)
        if col and col in row:
            row[col] = '' if val is None else val

    row['Allele Freq (T)'] = safe_freq(m.get('tumorAltCount'), m.get('tumorRefCount'))
    row['Allele Freq (N)'] = safe_freq(m.get('normalAltCount'), m.get('normalRefCount'))

    sc = sample_clin.get((stid, sid), {})
    for attr_id, val in sc.items():
        col = attr_display.get(attr_id)
        if col and col in row:
            row[col] = val

    pc = patient_clin.get((stid, pid), {})
    for attr_id, val in pc.items():
        col = attr_display.get(attr_id)
        if col and col in row:
            row[col] = val

    rows.append(row)

## 6. write base TSV
base_outfile = f"{gene_symbol.upper()}_mutations_final.tsv"
with open(base_outfile, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=TARGET_COLS_DEDUP, delimiter='\t',
                            extrasaction='ignore', restval='')
    writer.writeheader()
    writer.writerows(rows)
print(f"\nWritten: {base_outfile}  ({len(rows)} rows)")

## 7. Genome Nexus annotation
## Build unique genomic location objects for POST /annotation/genomic
## Each item: {chromosome, start, end, referenceAllele, variantAllele}
## Result order matches input order exactly, no key matching needed.

print("\n[GN] Preparing Genome Nexus queries...")

## Deduplicate variants; track which unique index each row maps to
loc_list = []     # list of unique genomic location dicts
loc_index = {}    # (chrom, start, end, ref, alt) -> index in loc_list
row_loc_idx = []  # parallel to rows: index into loc_list (-1 if unusable)

for row in rows:
    chrom = str(row.get('Chromosome', '')).strip()
    start = str(row.get('Start Pos', '')).strip()
    end   = str(row.get('End Pos', '')).strip()
    ref   = str(row.get('Ref', '')).strip()
    alt   = str(row.get('Var', '')).strip()

    if not (chrom and start and end and ref and alt):
        row_loc_idx.append(-1)
        continue

    # normalise chromosome (strip 'chr' prefix if present)
    chrom_norm = chrom.replace('chr', '').replace('Chr', '')

    key = (chrom_norm, start, end, ref, alt)
    if key not in loc_index:
        loc_index[key] = len(loc_list)
        loc_list.append({
            "chromosome":    chrom_norm,
            "start":         int(start),
            "end":           int(end),
            "referenceAllele": ref,
            "variantAllele":   alt,
        })
    row_loc_idx.append(loc_index[key])

print(f"  {len(loc_list)} unique variants to annotate across {len(rows)} rows")

# Genome Nexus fields to request
GN_FIELDS = "annotation_summary,clinvar,my_variant_info,signal,hotspots,mutation_assessor"
if oncokb_token:
    GN_FIELDS += ",oncokb"

# Token param (JSON string as required by GN)
gn_token_param = ''
if oncokb_token:
    gn_token_param = json.dumps({"oncokb": oncokb_token})

GN_CHUNK = 100
gn_results = [None] * len(loc_list)   # pre-allocate by index

print(f"  POSTing to {BASE_GN}/annotation/genomic in chunks of {GN_CHUNK}...")
t0 = time.time()
errors = 0

for chunk_start in range(0, len(loc_list), GN_CHUNK):
    chunk_locs = loc_list[chunk_start:chunk_start + GN_CHUNK]
    chunk_end  = chunk_start + len(chunk_locs)

    sys.stdout.write(
        f"\r  chunk {chunk_start//GN_CHUNK + 1}/"
        f"{(len(loc_list)-1)//GN_CHUNK + 1}"
        f"  ({chunk_end}/{len(loc_list)})   "
    )
    sys.stdout.flush()

    params = {
        "isoformOverrideSource": "mskcc",
        "fields": GN_FIELDS,
    }
    if gn_token_param:
        params["token"] = gn_token_param

    try:
        r = requests.post(
            f"{BASE_GN}/annotation/genomic",
            json=chunk_locs,
            params=params,
            headers={"Content-Type": "application/json",
                     "Accept": "application/json"},
            timeout=120,
        )
        if r.status_code != 200 or not r.text.strip():
            print(f"\n  GN error {r.status_code}: {r.text[:200]}")
            errors += 1
            time.sleep(2)
            continue
        chunk_results = r.json()
    except Exception as ex:
        print(f"\n  GN request exception: {ex}")
        errors += 1
        time.sleep(2)
        continue

    # results are in same order as input
    for i, res in enumerate(chunk_results):
        gn_results[chunk_start + i] = res

    time.sleep(0.3)

print(f"\n  done ({time.time()-t0:.1f}s)  errors={errors}")
annotated_count = sum(1 for x in gn_results if x is not None)
print(f"  variants with GN result: {annotated_count}/{len(loc_list)}")

## 8. extract annotations from GN results
def extract_gn(res):
    """Extract annotation columns from a single GN /annotation/genomic result."""
    out = {}
    if res is None:
        return out

    ## HGVSg: top-level
    out['HGVSg'] = res.get('hgvsg', '') or ''

    ## HGVSc, Exon
    anno   = res.get('annotation_summary', {}) or {}
    tx_sum = anno.get('transcriptConsequenceSummary', {}) or {}
    out['HGVSc'] = tx_sum.get('hgvsc', '') or ''
    out['Exon']  = tx_sum.get('exon',  '') or ''

    ## reVUE
    out['_isVue'] = 'yes' if tx_sum.get('isVue') else 'no'

    ## Functional Impact (browser format: "MutationAssessor: X;SIFT: X;Polyphen-2: X;AlphaMissense: X")
    ma_pred  = deep_get(res,    'mutation_assessor', 'functionalImpactPrediction') or 'NA'
    sift     = tx_sum.get('siftPrediction',    '') or 'NA'
    polyphen = tx_sum.get('polyphenPrediction','') or 'NA'
    alpha    = deep_get(tx_sum, 'alphaMissense', 'pathogenicity') or 'NA'
    out['Functional Impact'] = (
        f"MutationAssessor: {ma_pred};"
        f"SIFT: {sift};"
        f"Polyphen-2: {polyphen};"
        f"AlphaMissense: {alpha}"
    )

    ## gnomAD
    mvi = (res.get('my_variant_info', {}) or {}).get('annotation', {}) or {}
    af_exome  = deep_get(mvi, 'gnomadExome',  'alleleFrequency', 'af')
    af_genome = deep_get(mvi, 'gnomadGenome', 'alleleFrequency', 'af')
    out['gnomAD'] = af_exome if af_exome not in ('', None) else af_genome

    ## ClinVar
    clinsig = ''
    cv_mvi = mvi.get('clinVar', {}) or {}
    rcv = cv_mvi.get('rcv', [])
    if isinstance(rcv, list) and rcv:
        sigs = list(dict.fromkeys(
            r.get('clinicalSignificance', '') for r in rcv
            if r.get('clinicalSignificance')
        ))
        clinsig = '; '.join(sigs)
    elif isinstance(rcv, dict) and rcv.get('clinicalSignificance'):
        clinsig = rcv['clinicalSignificance']
    if not clinsig:
        clinsig = deep_get(res, 'clinvar', 'annotation', 'clinicalSignificance')
    out['ClinVar'] = clinsig

    ## dbSNP
    out['dbSNP'] = deep_get(mvi, 'dbsnp', 'rsid')

    ## SIGNAL: generalPopulationStats.frequencies.impact is the germline carrier
    ## frequency in the MSK-IMPACT cohort (e.g. 0.000058), matching the browser.
    sig_list = (res.get('signalAnnotation', {}) or {}).get('annotation', []) or []
    if sig_list:
        freq = deep_get(sig_list[0], 'generalPopulationStats', 'frequencies', 'impact')
        out['SIGNAL'] = str(freq) if freq not in ('', None) else ''

    ## CancerHotspot and 3DHotspot
    hs = (res.get('hotspots', {}) or {}).get('annotation', []) or []
    has_hotspot = any(bool(h) for h in hs)
    out['_cancerHotspot'] = 'yes' if has_hotspot else 'no'
    has_3d = False
    for h_group in hs:
        items = h_group if isinstance(h_group, list) else [h_group]
        for h in items:
            if isinstance(h, dict) and '3d' in str(h.get('type', '')).lower():
                has_3d = True
                break
        if has_3d:
            break
    out['_3dHotspot'] = 'yes' if has_3d else 'no'

    ## OncoKB
    if oncokb_token:
        okb = (res.get('oncokb', {}) or {}).get('annotation', {}) or {}
        parts = []
        oncogenic = okb.get('oncogenic', '') or ''
        if oncogenic:
            parts.append(f"OncoKB: {oncogenic}")
        sens = okb.get('highestSensitiveLevel', '') or ''
        if sens:
            parts.append(sens.replace('LEVEL_', 'level_').lower())
        resist = okb.get('highestResistanceLevel', '') or ''
        resist_str = resist.replace('LEVEL_R', 'R').lower() if resist else 'NA'
        parts.append(f"resistance {resist_str}")
        effect = deep_get(okb, 'mutationEffect', 'knownEffect')
        if effect:
            parts.append(f"{effect}")
        out['_oncokbStr'] = ', '.join(parts)
    else:
        out['_oncokbStr'] = ''

    return out


## 9. apply GN annotations to rows
print("Applying annotations to rows...")
for i, row in enumerate(rows):
    idx = row_loc_idx[i]
    ann = extract_gn(gn_results[idx]) if idx >= 0 else {}

    ## direct column fills
    for col, val in ann.items():
        if col.startswith('_'):
            continue
        if col in row and val not in ('', None):
            if col == 'dbSNP' and not val and row.get('dbSNP'):
                continue
            row[col] = val

    ## build Annotation string
    parts = []
    oncokb_str = ann.get('_oncokbStr', '')
    if oncokb_str:
        parts.append(oncokb_str)
    parts.append(f"reVUE: {ann.get('_isVue', 'no')}")
    parts.append(f"CancerHotspot: {ann.get('_cancerHotspot', 'no')}")
    parts.append(f"3DHotspot: {ann.get('_3dHotspot', 'no')}")
    row['Annotation'] = ';'.join(parts)


## 10. write annotated TSV
ann_outfile = f"annotated.{gene_symbol.upper()}_mutations_final.tsv"
with open(ann_outfile, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=TARGET_COLS_DEDUP, delimiter='\t',
                            extrasaction='ignore', restval='')
    writer.writeheader()
    writer.writerows(rows)

print(f"\nDone:")
print(f"  {base_outfile}             --> cBioPortal only")
print(f"  {ann_outfile}              --> + Genome Nexus")
print(f"  Rows: {len(rows)}  |  Columns: {len(TARGET_COLS_DEDUP)}")
if not oncokb_token:
    print("\nAnnotation column is empty. Re-run with your OncoKB token to populate it.")
print("Still empty: Custom Driver, Custom Driver Tiers, Copy #, # Mut in Sample")
print("CIViC intentionally excluded for now - add separately if needed.")