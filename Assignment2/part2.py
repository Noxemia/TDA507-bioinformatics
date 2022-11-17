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

CDIST_MIN = 3.79
CDIST_MAX = 3.81

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

def walk(atom: Atom, dir: bool, visited: list, chain: list):
    # either add the current atom number to the start or end of the chain
    if dir:
        chain.append(atom.i)
    else:
        chain.insert(0, atom.i)
    # add the current atom to the visited
    visited.append(atom.i)

    # see if there still are unvisited elements in nc
    cont = False
    for nc in atom.nc:
        if not nc in visited:
            cont = True
            break

    if not cont or len(chain) > 50:
        #if dir:
          #walk(atoms[chain[0] - 1], False, [cand.i], chain)
        #else:
            chains.append(chain)
            return


    # if the first closest atom is in visited we continue the walk on the other, and vice versa
    for _nc in atom.nc:
        # if we have seen this atom specifically, skip
        if _nc in visited: 
            continue
        # if not recurse down
        else:
            walk(_nc, dir, copy(visited), copy(chain))


chains = []

for cand in candidates:
    # walk in both direction for a candiate
    chain = walk(cand, True, [cand.i], [cand.i] )
    #pop last chain and use that with walk
    #chain = walk(cand, False, [cand.i], chains.pop())



maxlen = 0
for chain in chains:
    if len(chain) > maxlen:
        maxlen = len(chain)

print(maxlen)
maxlen = 0
for cand in candidates:
    if len(cand.nc) > maxlen:
        maxlen = len(cand.nc)

print(maxlen)