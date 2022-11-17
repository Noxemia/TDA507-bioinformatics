from copy import copy
from dataclasses import dataclass
from math import sqrt

data = []

with open("data_q2.txt", "r") as reader:
    for line in reader:
        line = line.strip()
        data.append(line)


    #dataclass for an atom
@dataclass
class Atom:
    x: float
    y: float
    z: float
    i: int
    nc: list 


#convert raw str data into atom class list
atoms = []
for line in data:
    line = list(line)
    
    i = ''
    for c in line:
        if c == '\t': break
        i += c
    line = line[len(i)+1:]

    x = ''
    for c in line:
        if c == ' ': break
        x += c
    line = line[len(x)+1:]

    y = ''
    for c in line:
        if c == ' ': break
        y += c
    line = line[len(y)+1:]

    z = ''
    for c in line:
        if c == ' ': break
        z += c
    atoms.append(Atom(float(x), float(y), float(z), int(i), []))

def distance(a1: Atom, a2: Atom) -> float:
    x = (a2.x - a1.x)**2
    y = (a2.y - a1.y)**2
    z = (a2.z - a1.z)**2
    return sqrt(x + y + z)

CDIST_MIN = 3.75
CDIST_MAX = 3.86

candidates = []

# add all candidate a-carbon neighbours to an atoms list
for i, atom in enumerate(atoms):

    for atom2 in atoms:
        if atom == atom2: continue
        dist = distance(atom, atom2)
        if dist > CDIST_MIN and dist < CDIST_MAX:
            atom.nc.append(atom2)

for atom in atoms:
    if len(atom.nc) != 0:
        candidates.append(atom)

def walk(atom: Atom):
    chains = []
    

    for _nc in atom.nc:
        chains.append([atom, _nc])

    
for i in range(5):
    print(candidates[i], "\n")
         