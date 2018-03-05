import sys

try:
    wanted = set(int(_) for _ in sys.argv[1:])
except ValueError:
    wanted = None

if not wanted:
    sys.exit("ERROR: Supply one or more NCBI taxonomy identifiers, space separated.\n")


print("Filtering NCBI taxonomy files nodes.dmp and names.dmp")
print("Will create nodes_.dmp and names_.dmp using just the given")
print("%i entries and their parent nodes." % len(wanted))

tree = dict()
with open("nodes.dmp") as handle:
    for line in handle:
        part = line.split("\t|\t", 2)
        taxid = int(part[0].strip())
        parent = int(part[1].strip())
        tree[taxid] = parent
print("Loaded %i entries from nodes.dmp" % len(tree))

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

with open("nodes.dmp") as handle:
    with open("nodes_.dmp", "w") as output:
        for line in handle:
            part = line.split("\t|\t", 1)
            taxid = int(part[0].strip())
            if taxid in include:
                output.write(line)
print("Created nodes_.dmp")


with open("names.dmp") as handle:
    with open("names_.dmp", "w") as output:
        for line in handle:
            part = line.split("\t|\t", 1)
            taxid = int(part[0].strip())
            if taxid in include:
                output.write(line)
print("Created names_.dmp")
