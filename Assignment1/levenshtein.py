# All set variables
X = "AGTCC"
Y = "CGCTCA"
xlist = list(X)
ylist = list(Y)

#Follows the implementation given in the wikipedia article on levenshtein distance
def lev(a:str, b:str):
    # if one string is shorter than the other, we will need the length difference in edits so we return that
    if len(a) == 0: return len(b)
    if len(b) == 0: return len(a)
    # if the letters are the same we just progress since no operation is needed
    if a[0] == b[0]: return lev(a[1:], b[1:])
    # We then check the recursive chain with the least amount of operations but we add one since we would need an edit where we are
    return min(
        #remove a letter in either of the sequences
        lev(a[1:], b) +1, 
        lev(a, b[1:]) +1 , 
        # edit the letter and therefore continue
        lev(a[1:], b[1:]) + 1
    )
 
print(X, " --- ", Y)     
print(lev(xlist,ylist))
