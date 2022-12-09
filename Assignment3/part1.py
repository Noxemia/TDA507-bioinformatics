from dataclasses import dataclass
from math import sqrt
import sys

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


atoms1 = handlefile(sys.argv[2])
atoms2 = handlefile(sys.argv[1])

def distance(a1: Atom, a2: Atom) -> float:
    x = (a2.x - a1.x)**2
    y = (a2.y - a1.y)**2
    z = (a2.z - a1.z)**2
    return sqrt(x + y + z)

overlaps = {}


comps = 0
for atom1 in atoms1:
    for atom2 in atoms2:
        comps += 1
        if distance(atom1, atom2) < 4.0:
            overlaps[atom1.i] = (atom1, atom2)
            break

        
for atomi in sorted(overlaps.keys()):
    print(atomi)
print(f"Number of overlaps: {len(overlaps.keys())}")
print(f"Number of comparissons: {comps}")