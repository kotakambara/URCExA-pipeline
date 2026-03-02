#!/usr/bin/env python
# coding: utf-8
"""
Usage:
    python get_expt_from_srr.py SRR12345678

標準出力:
    SRX / ERX / DRX 系の Experiment アクセッション ID を 1 行で出力します。
    見つからない場合は標準エラーにメッセージを出して終了コード 1 で終了します。
"""

import sys
import re
import html
import urllib.request
from urllib.error import URLError, HTTPError
from bs4 import BeautifulSoup

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

def fetch_xml(url: str) -> str:
    """指定 URL から XML テキストを取得"""
    try:
        with urllib.request.urlopen(url) as res:
            return res.read().decode("utf-8")
    except (URLError, HTTPError) as e:
        raise RuntimeError(f"HTTP error while fetching {url}: {e}")

def get_uid_from_srr(srr: str) -> str:
    """SRR / ERR / DRR から SRA の UID を取得"""
    url = f"{BASE}esearch.fcgi?db=sra&term={srr}[Accession]&retmode=xml"
    xml = fetch_xml(url)
    soup = BeautifulSoup(xml, "xml")

    id_tag = soup.find("Id")
    if not id_tag or not id_tag.text.strip():
        raise RuntimeError(f"Could not find UID for accession {srr}")
    return id_tag.text.strip()

def _validate_expt_id(expt: str) -> str:
    """Experiment ID の形式をチェックして返す"""
    expt = expt.strip()
    if not re.match(r"^[SED]RX[0-9]+$", expt):
        raise RuntimeError(f"Invalid Experiment ID extracted: {expt!r}")
    return expt

def get_experiment_from_uid(uid: str) -> str:
    """SRA UID から Experiment アクセッション (SRX/ERX/DRX) を取得"""
    url = f"{BASE}esummary.fcgi?db=sra&id={uid}&retmode=xml"
    xml = fetch_xml(url)
    soup = BeautifulSoup(xml, "xml")

    # 1) まずは <Item Name="Experiment">ERX11882943</Item> があればそれを使う
    exp_item = soup.find("Item", {"Name": "Experiment"})
    if exp_item and exp_item.text and exp_item.text.strip():
        try:
            return _validate_expt_id(exp_item.text)
        except RuntimeError:
            # 形式がおかしければ ExpXml から取り直す
            pass

    # 2) 次に <Item Name="ExpXml"> 内の XML から <Experiment acc="..."> を探す
    expxml_item = soup.find("Item", {"Name": "ExpXml"})
    if not expxml_item:
        raise RuntimeError(f"ExpXml item not found for UID {uid}")

    # ExpXml は &lt;Summary&gt;... のようなエスケープ文字列であることが多い
    # text を一旦取り出して &lt; / &gt; を復元する
    expxml_raw = expxml_item.text or ""
    if not expxml_raw.strip():
        # まれに text ではなく子要素として XML が入っているケースもある
        expxml_raw = "".join(str(child) for child in expxml_item.children)

    if not expxml_raw.strip():
        raise RuntimeError(f"Empty ExpXml for UID {uid}")

    expxml_str = html.unescape(expxml_raw)

    # ExpXml は <Summary>...</Summary><Experiment .../>... のように
    # ルート要素を持たない断片なので、ダミーの root タグで囲んでから XML としてパースする
    wrapped = f"<root>{expxml_str}</root>"
    expxml_soup = BeautifulSoup(wrapped, "xml")

    exp_tag = expxml_soup.find("Experiment")
    if not exp_tag:
        raise RuntimeError(f"<Experiment> tag not found inside ExpXml for UID {uid}")

    expt_acc = (exp_tag.get("acc") or "").strip()
    if not expt_acc:
        raise RuntimeError(f"Experiment acc attribute not found in ExpXml for UID {uid}")

    return _validate_expt_id(expt_acc)

def main():
    if len(sys.argv) != 2:
        print("Usage: python get_expt_from_srr.py SRR_or_ERR", file=sys.stderr)
        sys.exit(1)

    srr = sys.argv[1].strip()
    try:
        uid = get_uid_from_srr(srr)
        expt = get_experiment_from_uid(uid)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    print(expt)

if __name__ == "__main__":
    main()
