# All set variables
X = "ATCGAT"
Y = "ATACGT"
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
# perform Needleman-Wunsch algorithm
# as well as update the traceback matrix
for i in range(1, m+1):
    for j in range(1,n+1):
        # calculate match/mismatch score
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
        
        # update matrix with the new calculated max score
        matrix[i][j] = score

matrixPrint()

ti = m
tj = n

alignedX = []
alignedY = []

# step back through trackback matrix and determine aligned strings
while trace[ti][tj] != STOP:
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

# Add the rest of the string with gaps in the direction of the longer sequence
for i in range(ti):
    alignedX.insert(0, X[ti-i])
    alignedY.insert(0, "-")
for i in range(tj):
    alignedY.insert(0, X[tj-i])
    alignedX.insert(0, "-")

pipes = []
count = 0
#Construct new string with | between each matching character
for i in range(len(alignedX)):
    if alignedX[i]==alignedY[i]:
        pipes.append("|")
        count +=1
    else:
        pipes.append(" ")

print("".join(alignedX))
print("".join(pipes))
print("".join(alignedY))

# for me it makes the most intuitive sense to count the number of matching characters
# and divide it by the full string length. An aligned sequence is expected to have gaps and with the point being aligning the strings 
# I think the characters that now match should be counted. 
# This would result in two strings that are the exact same would match but nothing else
print("{}{}{}".format((count/len(alignedX))*100, '%',' identity between the two strings'))
print("The hamming distance is {}".format(len(alignedX)-count))