import urllib.request
from urllib.parse import quote
from bs4 import BeautifulSoup
from xml.etree.ElementTree import Element, SubElement, ElementTree
import sys
import re
from itertools import zip_longest


# 引数で生物名（学名でも可）を受け取る
species_list = sys.argv[1:]

if not species_list:
    sys.stderr.write(
        "Usage: python API_RNA-seq.py \"pearl millet\" \"Cenchrus americanus\" ...\n"
    )
    sys.exit(1)

# 生物名ごとに (name[Organism] OR name[All Fields]) を作る
organism_terms = []
for sp in species_list:
    organism_terms.append(f"({sp}[Organism] OR {sp}[All Fields])")

# 全生物名を OR でつなぐ
organism_part = " OR ".join(organism_terms)

# RNA-seq/transcriptome 条件を AND で付与
search_term = f"({organism_part}) AND (RNA-seq[All Fields] OR transcriptome[All Fields])"

print(search_term, file=sys.stderr)

encoded_term = quote(search_term)

base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
query_url = base + "esearch.fcgi?db=sra&term={}&usehistory=y".format(encoded_term)

data = urllib.request.urlopen(query_url)
xml = data.read()
soup = BeautifulSoup(xml, "xml")

QueryKey = soup.find("QueryKey").text
WebEnv = soup.find("WebEnv").text

query_url = base + (
    "esummary.fcgi?db=sra&query_key={}&WebEnv={}&RetMax=1000000"
    "&rettype=abstract&retmode=text&idtype=acc"
).format(QueryKey, WebEnv)

data = urllib.request.urlopen(query_url)
xml = data.read()
soup = BeautifulSoup(xml, "xml")

out = soup.find_all("Item")
out = list(out)


with open("tmp", mode="w") as f:
    for i in out:
        f.write(str(i) + "\n")


def _extract_first(pattern: str, text: str) -> str:
    m = re.search(pattern, text)
    return m.group(1) if m else ""

def _extract_last(pattern: str, text: str) -> str:
    ms = re.findall(pattern, text)
    return ms[-1] if ms else ""


with open("tmp", "r", encoding="utf-8", errors="replace") as f:
    lines = [ln.rstrip("\n") for ln in f]

exp_lines = [ln for ln in lines if "ExpXml" in ln]
run_lines = [ln for ln in lines if "Runs" in ln]

with open("tmp_2", "w", encoding="utf-8") as f:
    for exp, run in zip_longest(exp_lines, run_lines, fillvalue=""):
        f.write(f"{exp} {run}\n")

SRR_list = []
biosample_list = []
bioproject_list = []
sample_name_list = []
library_type_list = []

with open("tmp_2", "r", encoding="utf-8", errors="replace") as f:
    for ln in f:
        ln = ln.rstrip("\n")

        SRR_list.append(_extract_last(r'Run acc="([^"]*)"', ln))

        biosample_list.append(
            _extract_first(r".*&lt;Biosample&gt;([^<]*)&lt;\/Biosample&gt;.*", ln)
        )
        bioproject_list.append(
            _extract_first(r".*&lt;Bioproject&gt;([^<]*)&lt;\/Bioproject&gt;.*", ln)
        )
        sample_name_list.append(
            _extract_first(r".*&lt;LIBRARY_NAME&gt;([^<]*)&lt;\/LIBRARY_NAME&gt;.*", ln)
        )
        library_type_list.append(
            _extract_first(r".*&lt;LIBRARY_LAYOUT&gt;&lt;(PAIRED|SINGLE).*", ln)
        )


def _write_list(path: str, arr):
    with open(path, "w", encoding="utf-8") as f:
        for x in arr:
            f.write(f"{x}\n")

_write_list("SRR_list_tmp", SRR_list)
_write_list("biosample_list_tmp", biosample_list)
_write_list("PRJNA_list_tmp", bioproject_list)
_write_list("sample_name_list_tmp", sample_name_list)
_write_list("Library_type_list_tmp", library_type_list)



rows = list(zip(SRR_list, biosample_list, bioproject_list, sample_name_list, library_type_list))


sorted_lines = sorted("\t".join(map(str, r)) for r in rows)

with open("tmp_list", "w", encoding="utf-8") as f:
    for line in sorted_lines:
        f.write(line + "\n")