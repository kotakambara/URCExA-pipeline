#!/usr/bin/env python
# coding: utf-8

import sys
import time
import urllib.request
from urllib.error import URLError, HTTPError
from bs4 import BeautifulSoup
from xml.etree.ElementTree import Element, SubElement, ElementTree  # 使っていないが互換性のため残す


BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def fetch_with_retry(url, retries=3, wait=5):
    """
    NCBI eutils からの取得をリトライ付きで行う。
    最終的に失敗した場合は SAMN_info_tmp を空で作成し、exit code 0 で終了する。
    （mapping_script.sh のパイプラインを止めないため）
    """
    for i in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                return resp.read()
        except (URLError, HTTPError) as e:
            sys.stderr.write(f"[WARN] Failed to fetch {url} (attempt {i+1}/{retries}): {e}\n")
            if i < retries - 1:
                time.sleep(wait)
            else:
                sys.stderr.write(f"[WARN] Giving up fetching {url}. Use NA metadata.\n")
                # 空ファイルを作って正常終了（bash側で NA として扱われる）
                with open("SAMN_info_tmp", "w") as f:
                    f.write("")
                sys.exit(0)


def main():
    if len(sys.argv) < 2:
        sys.stderr.write("Usage: findSAMN_info.py <SAMN_ID>\n")
        # SAMN が無い場合も、空ファイルを出して 0 終了にしておく
        with open("SAMN_info_tmp", "w") as f:
            f.write("")
        sys.exit(0)

    search_term = sys.argv[1]

    # --- esearch: QueryKey / WebEnv を取得 ---
    esearch_url = (
        BASE
        + f"esearch.fcgi?db=biosample&term={search_term}&usehistory=y"
    )

    xml = fetch_with_retry(esearch_url)
    soup = BeautifulSoup(xml, "xml")

    qkey_tag = soup.find("QueryKey")
    webenv_tag = soup.find("WebEnv")

    if qkey_tag is None or webenv_tag is None:
        sys.stderr.write(f"[WARN] QueryKey/WebEnv not found for {search_term}. Use NA metadata.\n")
        with open("SAMN_info_tmp", "w") as f:
            f.write("")
        sys.exit(0)

    QueryKey = qkey_tag.text
    WebEnv = webenv_tag.text

    # --- esummary: BioSample の中身を取得 ---
    esummary_url = (
        BASE
        + "esummary.fcgi?db=biosample"
        + f"&query_key={QueryKey}"
        + f"&WebEnv={WebEnv}"
        + "&RetMax=10000&rettype=abstract&retmode=text&idtype=acc"
    )

    xml = fetch_with_retry(esummary_url)
    soup = BeautifulSoup(xml, "xml")

    # もともとのスクリプトと互換性を保つため、soup を list にして 1 行ずつ書き出す
    out = list(soup)

    with open("SAMN_info_tmp", mode="w") as f:
        for i in out:
            f.write(str(i) + "\n")


if __name__ == "__main__":
    main()
