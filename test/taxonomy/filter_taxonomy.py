#!/usr/bin/env python
# TODO: Proper API, including option to set paths for the input and output
# locations to make dealing with the filenames easier

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

with open("/tmp/nodes.dmp") as handle:
    with open("nodes.dmp", "w") as output:
        for line in handle:
            part = line.split("\t|\t", 1)
            taxid = int(part[0].strip())
            if taxid in include:
                output.write(line)
print("Created nodes.dmp")


with open("/tmp/names.dmp") as handle:
    with open("names.dmp", "w") as output:
        for line in handle:
            part = line.split("\t|\t", 1)
            taxid = int(part[0].strip())
            if taxid in include:
                output.write(line)
print("Created names.dmp")
