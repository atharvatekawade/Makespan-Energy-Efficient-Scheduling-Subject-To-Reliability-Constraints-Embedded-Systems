import math
import config, utils

def mecrg(R, l = []):
    if(l == []):
        l, _ = utils.urv()
    
    alloc = []
    avail = []
    rels = []
    rmax = []

    for i in range(len(config.graph)):
        alloc.append(-1)
        rels.append(1)
        r = 0
        for pr in range(len(config.times)):
            r = max(r, math.exp(-config.L[pr]*config.times[pr][i]))
        rmax.append(r)
    
    for i in range(len(config.times)):
        avail.append(0)
    
    max_rel = math.prod(rmax)

    for i in range(len(l)):
        max_rel /= rmax[l[i]]
        r_goal = R/(math.prod(rels)*max_rel)
        emin = float('inf')
        pmin = -1
        for pr in range(len(config.times)):
            if(math.exp(-config.L[pr]*config.times[pr][l[i]]) >= r_goal):
                e = (config.PS[pr] + config.C[pr]*(1**config.A[pr]))*config.times[pr][l[i]]/config.F[pr]
                if(e < emin):
                    pmin = pr
                    emin = e
        
        if(pmin == -1):
            print("Wrong processor rel constraint violated")
            quit()
        
        else:
            rels[l[i]] = math.exp(-config.L[pmin]*config.times[pmin][l[i]])
            alloc[l[i]] = pmin
    
    return alloc
