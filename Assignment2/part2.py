from copy import copy, deepcopy
from dataclasses import dataclass
from math import sqrt

data = []

with open("data_q2.txt", "r") as reader:
    for line in reader:
        line = line.strip()
        data.append(line)

    # dataclass for an atom


@dataclass
class Atom:
    x: float
    y: float
    z: float
    i: int
    nc: list
    bc: list

    def __repr__(self) -> str:
        return (
            f"Atom(x:{self.x}, y: {self.y}, z:{self.z}, i:{self.i}, nc:{len(self.nc)})"
        )

    def getId(self) -> int:
        return self.i


# convert raw str data into atom class list
atoms = []
for line in data:
    line = list(line)

    i = ""
    for c in line:
        if c == "\t":
            break
        i += c
    line = line[len(i) + 1 :]

    x = ""
    for c in line:
        if c == " ":
            break
        x += c
    line = line[len(x) + 1 :]

    y = ""
    for c in line:
        if c == " ":
            break
        y += c
    line = line[len(y) + 1 :]

    z = ""
    for c in line:
        if c == " ":
            break
        z += c
    atoms.append(Atom(float(x), float(y), float(z), int(i), [], []))


def distance(a1: Atom, a2: Atom) -> float:
    x = (a2.x - a1.x) ** 2
    y = (a2.y - a1.y) ** 2
    z = (a2.z - a1.z) ** 2
    return sqrt(x + y + z)


CDIST_MIN = 3.77
CDIST_MAX = 3.89

candidates = []

# add all candidate a-carbon neighbours to an atoms list
for i, atom in enumerate(atoms):

    for atom2 in atoms:
        if not atom.getId() == atom2.getId():
            dist = distance(atom, atom2)
            if  dist < CDIST_MIN:
                atom.bc.append(atom2)

# add all candidate a-carbon neighbours to an atoms list
for i, atom in enumerate(atoms):

    for atom2 in atoms:
        if atom == atom2:
            continue
        dist = distance(atom, atom2)
        if dist > CDIST_MIN and dist < CDIST_MAX:
            cnt = 0
            for _bc in atom.bc:
                if _bc in atom2.bc:
                    cnt += 1
            if cnt > 1:
                atom.nc.append(atom2)


for atom in atoms:
    if len(atom.nc) != 0:
        candidates.append(atom)

chains = []

candFilt = []
for cand in candidates:
    if len(cand.nc) == 1:
        candFilt.append(cand)


def walk(atom: Atom, chain: list = []):

    
    # While we only have 2 neighbours we continue
    while len(atom.nc) == 2:
        chain.append(atom)
        if not atom.nc[0] in chain:
            atom = atom.nc[0]
        elif not atom.nc[1] in chain:
            atom = atom.nc[1]
        else:
            chains.append(chain)
            return
    
    # if the 2s ended in a 1
    if len(atom.nc) == 1 and len(chain) != 0:
        chain.append(atom)
        chains.append(chain)
        return
    
    
    chain.append(atom)       
    # When we have more than 2 neighbours we recurse on each of them 
    for _nc in atom.nc:
        if _nc in chain: continue
        walk(_nc, copy(chain))
        
        
for cand in candFilt:
    walk(cand, [])
    
maxlen = 0
for chain in chains:
    if (tmp := len(chain)) > maxlen:
        maxlen = tmp

print(maxlen)