#!/usr/bin/env python3
import sys
import biometal

bam_path = sys.argv[1]
count = 0
for record in biometal.BamReader.from_path(bam_path):
    count += 1
