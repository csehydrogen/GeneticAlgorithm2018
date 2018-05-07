a = ["unweighted_50", "unweighted_100", "unweighted_500", "weighted_500", "weighted_chimera_297"]

from subprocess import call
for aa in a:
    for b in range(30):
        call(["thorq", "--add", "./a", "{}.txt".format(aa), "{}_result_{}.txt".format(aa, b)])
