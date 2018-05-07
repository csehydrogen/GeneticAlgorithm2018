a = ["unweighted_50", "unweighted_100", "unweighted_500", "weighted_500", "weighted_chimera_297"]

for aa in a:
    d = []
    for b in range(30):
        with open("{}_result_{}.txt".format(aa, b)) as f:
            d.append(f.read().strip())
    print ",".join(d)
