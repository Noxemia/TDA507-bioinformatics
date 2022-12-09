# All set variables
X = "AGTCC"
Y = "CGCTCA"
STOP = 0
UP = 1
LEFT = 2
DIAG = 3
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

score = 0
# perform Smith-Waterman algorithm
# as well as update the traceback matrix
for i in range(1, m+1):
    for j in range(1,n+1):
        if X[i-1] == Y[j-1]:
            score = matrix[i-1][j-1] + matchScore
        else:
            score = matrix[i-1][j-1] + mismatchScore
        trace[i][j] = DIAG

        #compare match/mismatch score to gapped score and store the largest
        tmp = matrix[i-1][j] - gapPenatly
        if tmp > score:
            score = tmp
            trace[i][j] = UP
        tmp = matrix[i][j-1] - gapPenatly
        if tmp > score:
            score = tmp
            trace[i][j] = LEFT
        
        # update matrix with the new calculated max score as well as see if the score is less than 0 in accordance to
        # the Smith-Waterman extention on Needleman-Wunsch
        if score < 0: score = 0
        matrix[i][j] = score

matrixPrint()

#Coordinates for traceback
ti = 0
tj = 0

#Find the highest score in the score matrix and set the begning trace coords to that cell.
highest = 0
for ci in range(len(matrix)-1):
    for cj in range(len(matrix[0])-1):
        if matrix[ci][cj] > highest:
            highest = matrix[ci][cj]
            ti = ci
            tj = cj

alignedX = []
alignedY = []


# step back through trackback matrix and determine locally aligned strings
while trace[ti][tj] != STOP:
    if matrix[ti][tj] == 0: break
    elem = trace[ti][tj]
    if elem == DIAG:
        alignedX.insert(0, X[ti-1])
        alignedY.insert(0, Y[tj-1])
        ti -= 1
        tj -= 1
    elif elem == LEFT:
        alignedX.insert(0, "-")
        alignedY.insert(0, Y[tj-1])
        tj = tj -  1
    elif elem == UP:
        alignedX.insert(0, X[ti-1])
        alignedY.insert(0, "-")
        ti = ti - 1
    
print("".join(alignedX))
print("".join(alignedY))