
def rotations(t):
    # Return list of rotations of input string t
    tt = t * 2
    return [ tt[i:i+len(t)] for i in range(0, len(t)) 

print(rotations('cat'))
#['cat', 'atc', 'tca']


def bwm(t):
    # Return lexicographically sorted list of t’s rotations
    return sorted(rotations(t))

# bwm('abaaba$')


##get matrix:   print('\n'.join(bwm('abaaba$')))

def bwtViaBwm(t):
    # Given T, returns BWT(T) by way of the BWM
    return ''.join(map(lambda x: x[-1], bwm(t)))


def suffixArray(s):
    satups = sorted([(s[i:], i) for i in xrange(0, len(s))] )
    return map(lambda x: x[1], satups)

def bwtViaSa(t):
    # Given T, returns BWT(T) by way of the suffix array
    bw = []
    for si in suffixArray(t):
        if si == 0:
            bw.append('$')
        else:
            bw.append(t[si-1])
    return ''.join(bw) # return string-ized version of list bw


#  bwtViaBwm('abaaba$'), bwtViaSa('abaaba$') # same result
#    ('abba$aa', 'abba$aa')
