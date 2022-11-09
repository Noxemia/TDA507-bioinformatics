from cmath import nan
from dataclasses import dataclass
from math import sqrt
from xmlrpc.client import MAXINT
# Read in data in forms of lines
data = []

with open("test_q1.txt", "r") as reader:
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
    closest: list = (None, None)

    def updateClosest(self, c1, c2):
        self.closest = (c1, c2)




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
    atoms.append(Atom(float(x), float(y), float(z), int(i)))

def distance(a1: Atom, a2: Atom) -> float:
    x = (a2.x - a1.x)**2
    y = (a2.y - a1.y)**2
    z = (a2.z - a1.z)**2
    return sqrt(x + y + z)

for atom in atoms:
    distances = []
    
    #Get the distance between the current atom and all others
    for atom2 in atoms:
        if atom2.i == atom.i: continue
        dist = distance(atom, atom2)
        distances.append((atom2.i, dist))

    close1 = [MAXINT, MAXINT]
    close2 = [MAXINT, MAXINT]

    # Find the closest two, unordered
    for distPair in distances:
        if distPair[1] < close1[1] or distPair[1] < close2[1]:
            if close1[1] < close2[1]:
                close2 = distPair
            else:
                close1 = distPair

    #update atom with its closest ones
    atom.updateClosest(close1[0], close2[0])
    

# we track through the order by starting at atom 1 and then walking up the chain by keeping track of which atoms we have visited
# dir is used to see if we are walking "up" or "down" the chain
order = [1]
distances = []
last = 1
def walk(atom: Atom, dir: bool, visited: list):
    # either add the current atom number to the start or end of the chain
    if dir:
        order.append(atom.i)
    else:
        order.insert(0, atom.i)
    # add the current atom to the visited
    visited.append(atom.i)


    # if the both the closest atoms are visited we stop walking the chain
    if atom.closest[0] in visited and atom.closest[1] in visited: return

    # if the first closest atom is in visited we continue the walk on the other, and vice versa
    if atom.closest[0] in visited:
        walk(atoms[atom.closest[1]-1 ], dir, visited)
    elif atom.closest[1] in visited:
        walk(atoms[atom.closest[0]-1 ], dir, visited)


# walk both directions with atom 1 as the start
walk(atoms[atoms[0].closest[0] - 1], True, [1] )
walk(atoms[atoms[0].closest[1] - 1], False, [1] )
print(order)

#find the average lengths
for i in range(len(order) -2):
    print(distance(atoms[order[i] -1], atoms[order[i+1] -1]))
    