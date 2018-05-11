#!/usr/bin/env python
# TODO: Proper API, including option to set paths for the input and output
# locations to make dealing with the filenames easier

import os
import sys

try:
    wanted = set(int(_) for _ in sys.argv[1:])
except ValueError:
    wanted = None

if not wanted:
    sys.exit("ERROR: Supply one or more NCBI taxonomy identifiers, space separated.\n")


print("Filtering NCBI taxonomy files /tmp/nodes.dmp and names.dmp etc")
print("Will create ./nodes.dmp and ./names.dmp etc using just the given")
print("%i entries and their parent nodes." % len(wanted))

tree = dict()
with open("/tmp/nodes.dmp") as handle:
    for line in handle:
        part = line.split("\t|\t", 2)
        taxid = int(part[0].strip())
        parent = int(part[1].strip())
        tree[taxid] = parent
print("Loaded %i entries from /tmp/nodes.dmp" % len(tree))

include = set()
for taxid in wanted:
    include.add(taxid)
    while True:
        parent = tree[taxid]
        if parent == taxid:
            # Reached root node
            break
        if parent in include:
            # Short cut
            break
        include.add(parent)
        taxid = parent
print("Expanded %i given TaxID to a list of %i including ancestors"
      % (len(wanted), len(include)))

# Can ignore delnodes.dmp
with open("delnodes.dmp", "w"):
    pass

for name in ("merged.dmp", "names.dmp", "nodes.dmp"):
    print("Filtering %s" % name)
    # Using binary mode to handle encoding of citations.dmp
    # Python 2 doesn't support open(..., encoding="latin1")
    with open(os.path.join("/tmp", name), "rb") as handle:
        with open(name, "wb") as output:
            for line in handle:
                part = line.decode("latin1").split("\t|\t", 1)
                taxid = int(part[0].strip())
                if taxid in include:
                    output.write(line)

print("Checking which genetic codes and divisions are needed")
wanted_gencode = set()
wanted_division = set()
with open("/tmp/nodes.dmp") as handle:
    for line in handle:
        part = line.split("\t|\t")
        taxid = int(part[0].strip())
        if taxid in include:
            wanted_gencode.add(int(part[6].strip()))
            wanted_division.add(int(part[4].strip()))
print("Only want these genetic codes: "
      + ",".join(str(_) for _ in sorted(wanted_gencode)))
print("Only want these divisions: "
      + ",".join(str(_) for _ in sorted(wanted_division)))

print("Filtering gencode.dmp")
with open("/tmp/gencode.dmp", "rb") as handle:
    with open("gencode.dmp", "wb") as output:
        for line in handle:
            part = line.decode("latin1").split("\t|\t", 1)
            if int(part[0].strip()) in wanted_gencode:
                output.write(line)

print("Filtering division.dmp")
with open("/tmp/division.dmp", "rb") as handle:
    with open("division.dmp", "wb") as output:
        for line in handle:
            part = line.decode("latin1").split("\t|\t", 1)
            if int(part[0].strip()) in wanted_division:
                output.write(line)

print("Filtering citations.dmp")
with open("/tmp/citations.dmp", "rb") as handle:
    with open("citations.dmp", "wb") as output:
        for line in handle:
            # Columns 6 is space-separated taxids,
            # but expect trailing \t|\n as well.
            field = line.decode("latin1").split("\t|\t")[6].strip("\t\n| ")
            # Note at time of writing, there was a corrupt entry
            # 54234 | Dissanayake et al. (2017c) | ...
            # which had line breaks within the URL
            # http://www.mycosphere.org/pdf/Mycosphere_8_2_7.pdf
            # I opted to fix the input file manually.
            if include.intersection(int(_) for _ in field.split()):
                output.write(line)

print("Done")
