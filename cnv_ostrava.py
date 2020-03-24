#!/usr/local/bin/python3
# -*- coding: utf-8 -*-

from pathlib import Path
import csv

source = Path('all_samples_max_cov_noheader.tsv')

with source.open() as src:
    csvreader = csv.reader(src, dialect='excel-tab')

    # get number of columns and rewind
    cols = len(next(csvreader)[1:])
    src.seek(0)

    csvwriters = []

    # create a csv.writer for each column
    for i in range(cols):
        # output_col_01.tsv, output_col_02.tsv ...
        csvwriters.append(
            csv.writer(
                Path(f'output_col_{i + 1:02d}.tsv').open('w'),
                dialect='excel-tab'
            )
        )

    nan = float('nan')

    for name, *cols in csvreader:
        for i, a in enumerate(cols):
            row = [name]
            for j, b in enumerate(cols):
                # skip the quotient of a col by itself
                if i != j:
                    a = float(a)
                    b = float(b)
                    # nan if division by zero
                    row.append(round(a / b, 4) if b else nan)

            csvwriters[i].writerow(row)