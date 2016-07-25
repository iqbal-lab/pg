def rotations(t):
    ''' Return list of rotations of input string t '''
    tt = t * 2
    return [ tt[i:i+len(t)] for i in range(0, len(t)) ]

def bwm(t):
    ''' Return lexicographically sorted list of t’s rotations '''
    return sorted(rotations(t))

def bwtViaBwm(t):
    ''' Given T, returns BWT(T) by way of the BWM '''
    return ''.join(map(lambda x: x[-1], bwm(t)))

#t = 'abaaba$'
#b = bwtViaBwm(t)
#b
#'abba$aa'


def rankBwt(bw):
    ''' Given BWT string bw, return parallel list of B-ranks.  Also
        returns tots: map from character to # times it appears. '''
    tots = dict()
    ranks = []
    for c in bw:
        if c not in tots:
            tots[c] = 0
        ranks.append(tots[c])
        tots[c] += 1
    return ranks, tots

#ranks, tots = rankBwt(b)
#print zip(b, ranks) # print characters of BWT(T) in order, along with rank

# [('a', 0), ('b', 0), ('b', 1), ('a', 1), ('$', 0), ('a', 2), ('a', 3)]

def firstCol(tots):
    ''' Return map from character to the range of rows prefixed by
        the character. '''
    first = {}
    totc = 0
    for c, count in sorted(tots.iteritems()):
        first[c] = (totc, totc + count)
        totc += count
    return first

#firstCol(tots)
#{'$': (0, 1), 'a': (1, 5), 'b': (5, 7)}


# confirm that the representation of the first column above is sensible
#print('\n'.join(bwm(t)))
# $abaaba
# a$abaab
# aaba$ab
# aba$aba
# abaaba$
# ba$abaa
# baaba$a


def reverseBwt(bw):
    ''' Make T from BWT(T) '''
    ranks, tots = rankBwt(bw)
    first = firstCol(tots)
    rowi = 0 # start in first row
    t = '$' # start with rightmost character
    while bw[rowi] != '$':
        c = bw[rowi]
        t = c + t # prepend to answer
        # jump to row that starts with c of same rank
        rowi = first[c][0] + ranks[rowi]
    return t

#reverseBwt(b)
# 'abaaba$'

#reverseBwt(bwtViaBwm('In_the_jingle_jangle_morning$'))

'In_the_jingle_jangle_morning$'

