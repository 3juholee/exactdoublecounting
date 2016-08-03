which  = [31,33,35,37,38,39]
which = range(28,40)[::2]
which += range(0,28)[::4]
which += range(40,52)[::4]
which += range(52,70)[::6]
which += range(70,210)[::10]
which += range(90,110)[::2]
which += range(110,130)[::2]
which += range(130,150)[::2]
which += range(150,170)[::2]
which += range(170,190)[::2]


#which += range(190,300)[::10]
which.sort()
print len(which)
#print which
