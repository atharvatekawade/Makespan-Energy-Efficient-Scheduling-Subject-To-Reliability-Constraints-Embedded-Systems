import math
import config, utils

def max_re(R, l = []):
    if(l == []):
        l, _ = utils.urv()
    
    rels = []
    target = []
    for i in range(len(l)):
        r = 1
        target.append(0)
        rels.append(1)
    
    alloc = []
    for i in range(len(l)):
        alloc.append([])
    
    energy = 0
    exp = 1
    for i in range(len(l)):
        r_goal = R**(1/len(config.graph))
        metrics = []
        for pr in range(len(config.times)):
            e = (config.PS[pr] + config.C[pr]*(config.F[pr]**config.A[pr]))*config.times[pr][l[i]]/config.F[pr]
            lf = config.L[pr]*10**(config.D[pr]*(1-config.F[pr])/(1-config.F_MIN[pr]))
            r = 1-math.exp(-lf*config.times[pr][l[i]]/config.F[pr])
            metrics.append([r, e, pr])
        
        metrics.sort(key=lambda x: x[0])
        r = 1
        for k in range(len(metrics)):
            if(1-r >= r_goal):
                target[l[i]] = 1
                break

            r = r*metrics[k][0]
            energy += metrics[k][1]
            alloc[l[i]].append(metrics[k][2])

        
        rels[l[i]] = 1-r

    return alloc, math.prod(rels), energy