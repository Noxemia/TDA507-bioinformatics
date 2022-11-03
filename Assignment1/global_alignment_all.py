from copy import copy
# All set variables
X = "PAWHEAE"
Y = "HDAGAWGHEQ"
STOP = 0
UP = 1
LEFT = 2
DIAG = 3
UPLEFT = 4
UPDIAG = 5
LEFTDIAG = 6
ALL = 7
matchScore = 2
mismatchScore = -1
gapPenatly = 2
m = len(X)
n = len(Y)

#The score matrix as well as the trace matrix to be able to backtrace
matrix = []
trace = []

#init both matricies to zeroes
for _ in range(m+1):
    matrix.append([0] * (n+1))
for _ in range(m+1):
    trace.append([0] * (n+1))

#function to print a matrix in the same format as the example programs
def matrixPrint(mat=matrix):
    print("\t", end="")
    print("\t", end="")
    for c in list(Y):
        print("{}\t".format(c), end="")
    print("\n")
    for i, line in enumerate(mat):
        if(i != 0):
            print("{}\t".format(X[i-1]), end="")
        else:
            print("\t", end="")
        for c in line:
            print("{}\t".format(c), end="")
        print("\n")
        

#intialize the score matrix with the gap penatlies
for x in range(1, m+1):
    matrix[x][0] = matrix[x-1][0] - gapPenatly
for y in range(1, n+1):
    matrix[0][y] = matrix[0][y-1] - gapPenatly


# perform Needleman-Wunsch algorithm
# as well as update the traceback matrix
for i in range(1, m+1):
    for j in range(1,n+1):
        # calculate match/mismatch score
        diagScore = 0
        if X[i-1] == Y[j-1]:
            diagScore = matrix[i-1][j-1] + matchScore
        else:
            diagScore = matrix[i-1][j-1] + mismatchScore

        upScore = matrix[i-1][j] - gapPenatly
        leftScore = matrix[i][j-1] - gapPenatly

        # case all
        if diagScore == upScore and upScore == leftScore:
            matrix[i][j] = diagScore
            trace[i][j] = ALL
        # case left and diag
        elif diagScore == leftScore and diagScore > upScore:
            matrix[i][j] = diagScore
            trace[i][j] = LEFTDIAG
        # case up and diag
        elif diagScore == upScore and diagScore > leftScore:
            matrix[i][j] = diagScore
            trace[i][j] = UPDIAG
        # case up and left
        elif upScore == leftScore and upScore > diagScore:
            matrix[i][j] = upScore
            trace[i][j] = UPLEFT
        # case diag
        elif diagScore > leftScore and diagScore > upScore:
            matrix[i][j] = diagScore
            trace[i][j] = DIAG
        #case left
        elif leftScore > diagScore and leftScore > upScore:
            matrix[i][j] = leftScore
            trace[i][j] = LEFT
        #case uo
        elif upScore > diagScore and upScore > leftScore:
            matrix[i][j] = upScore
            trace[i][j] = UP

matrixPrint(trace)

allTraces = []
def fnTrace(ti=m, tj=n, alignedX=[], alignedY=[], command=None):
    
# step back through trackback matrix and determine aligned strings
    if trace[ti][tj] == STOP:
        for i in range(ti):
            alignedX.insert(0, X[ti-i])
            alignedY.insert(0, "-")
        for i in range(tj):
            alignedY.insert(0, X[tj-i])
            alignedX.insert(0, "-")
        allTraces.append((alignedX, alignedY))
        return
    elem = trace[ti][tj]
    # forward propagated command from branching
    if command: elem = command

    # standard functionality
    if elem == DIAG:
        alignedX.insert(0, X[ti-1])
        alignedY.insert(0, Y[tj-1])
        fnTrace(ti-1, tj-1, copy(alignedX), copy(alignedY) )
    elif elem == LEFT:
        alignedX.insert(0, "-")
        alignedY.insert(0, Y[tj-1])
        fnTrace(ti, tj-1, copy(alignedX), copy(alignedY) )
    elif elem == UP:
        alignedX.insert(0, X[ti-1])
        alignedY.insert(0, "-")
        fnTrace(ti-1, tj, copy(alignedX), copy(alignedY) )
    elif elem == UPLEFT:
        fnTrace(ti, tj, copy(alignedX), copy(alignedY), UP)
        fnTrace(ti, tj, copy(alignedX), copy(alignedY), LEFT)
    elif elem == UPDIAG:
        fnTrace(ti, tj, copy(alignedX), copy(alignedY), UP)
        fnTrace(ti, tj, copy(alignedX), copy(alignedY), DIAG)
    elif elem == LEFTDIAG:
        fnTrace(ti, tj, copy(alignedX), copy(alignedY), DIAG)
        fnTrace(ti, tj, copy(alignedX), copy(alignedY), LEFT)
    elif elem == ALL:
        fnTrace(ti, tj, copy(alignedX), copy(alignedY), UP)
        fnTrace(ti, tj, copy(alignedX), copy(alignedY), LEFT)
        fnTrace(ti, tj, copy(alignedX), copy(alignedY), DIAG)
        
fnTrace()

for i, pair in enumerate(allTraces):
    print(" ---------- ALIGNMENT {} ----------".format(i+1))
    print("".join(pair[0]))
    print("".join(pair[1]))
    
print("Total number of alignments", len(allTraces))