names = "/home/ae42909/Scratch/kaijuDB_viral/names.dmp"
nodes = "/home/ae42909/Scratch/kaijuDB_viral/nodes.dmp"


with open(nodes, 'r') as in_file:
    for line in in_file:
        if line[:6] == "75079\t":
            print line


            
    head = [next(in_file) for x in xrange(1)]
 
