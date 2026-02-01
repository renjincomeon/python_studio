#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Parse HTML table, extract Symbol column,
map to human HGNC official symbols,
remove NA and duplicates.
"""

from bs4 import BeautifulSoup
import mygene
from pathlib import Path
import csv


HTML_FILE = "showFullTableHTML.html"


def extract_symbols(html_file):
    """Extract Symbol column from HTML table"""
    with open(html_file, encoding="utf-8") as f:
        soup = BeautifulSoup(f, "html.parser")

    symbols = []
    for row in soup.select("table tbody tr"):
        tds = row.find_all("td")
        if not tds:
            continue
        symbol = tds[0].get_text(strip=True)
        if symbol:
            symbols.append(symbol)

    # 去重，保留顺序
    symbols = list(dict.fromkeys(symbols))
    return symbols


def map_to_human_symbols(symbols):
    """Map gene symbols to human HGNC official symbols"""
    mg = mygene.MyGeneInfo()

    results = mg.querymany(
        symbols,
        scopes=["symbol", "alias"],
        fields="symbol",
        species="human"
    )

    mapping = {}
    for r in results:
        query = r.get("query")
        if r.get("notfound", False):
            mapping[query] = "NA"
        else:
            mapping[query] = r.get("symbol", "NA")

    return mapping


def write_outputs(symbols, mapping):
    """Write output files"""

    # 1. raw symbols
    with open("symbol_raw.txt", "w") as f:
        for s in symbols:
            f.write(s + "\n")

    # 2. mapping table
    with open("symbol_human.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Original_Symbol", "Human_Official_Symbol"])
        for s in symbols:
            writer.writerow([s, mapping.get(s, "NA")])

    # 3. final gene list
    final_genes = []
    for v in mapping.values():
        if v != "NA" and v not in final_genes:
            final_genes.append(v)

    with open("symbol_human_final.txt", "w") as f:
        for g in final_genes:
            f.write(g + "\n")

    print(f"Raw symbols: {len(symbols)}")
    print(f"Final human genes: {len(final_genes)}")


def main():
    html_path = Path(HTML_FILE)
    if not html_path.exists():
        raise FileNotFoundError(f"{HTML_FILE} not found")

    symbols = extract_symbols(html_path)
    mapping = map_to_human_symbols(symbols)
    write_outputs(symbols, mapping)


if __name__ == "__main__":
    main()
