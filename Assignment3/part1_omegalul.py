from dataclasses import dataclass
from math import sqrt
from math import log
from collections import defaultdict


#dataclass for an atom
@dataclass
class Atom:
    x: float
    y: float
    z: float
    i: int


def handlefile(filename: str) -> list:
    dataraw = []

    with open(filename, "r") as reader:
        for line in reader:
            line = line.strip()
            dataraw.append(line)

    data = []
    for line in dataraw:
        if line[0:4] == "ATOM" or line[0:6] == "HETATM":
            data.append(line)
    atoms = []
    for line in data:
        tmp = line.split(" ")
        newline = []
        for s in tmp:
            if not s == '':
                newline.append(s)
        atoms.append(Atom(float(newline[6]), float(newline[7]), float(newline[8]), int(newline[1])))
    return atoms


atoms1 = handlefile("1cdh.pdb")
atoms2 = handlefile("2csn.pdb")

def distance(a1: Atom, a2: Atom) -> float:
    x = (a2.x - a1.x)**2
    y = (a2.y - a1.y)**2
    z = (a2.z - a1.z)**2
    return sqrt(x + y + z)

overlaps = {}

xcomps = defaultdict(list)
ycomps = defaultdict(list)
zcomps = defaultdict(list)

for atom in atoms1:
    tmp = round(atom.x)
    for i in range(tmp-6, tmp+6):
        xcomps[tmp+i].append(atom)
    tmp = round(atom.y)
    for i in range(tmp-6, tmp+6):
        ycomps[tmp+i].append(atom)
    tmp = round(atom.z)
    for i in range(tmp-6, tmp+6):
        zcomps[tmp+i].append(atom)

    
xmin = min(xcomps.keys())
xmax = max(xcomps.keys())


for atomne in atoms2:
    xi =  round(atomne.x)
    yi =  round(atomne.y)
    zi =  round(atomne.z)
    for atom in xcomps[xi] + ycomps[yi] + zcomps[zi]:
        if distance(atomne, atom) < 4.0:
            overlaps[atomne.i] = (atomne, atom)
            break

print(len(overlaps.keys()))
#print(comps)