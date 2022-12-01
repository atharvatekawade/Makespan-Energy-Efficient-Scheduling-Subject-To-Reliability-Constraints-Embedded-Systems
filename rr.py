import math
import config, utils

def rr(R, l = []):
    if(l == []):
        l, _ = utils.urv()
    
    rels = []
    rexp = []
    target = []
    r1 = 1
    for i in range(len(l)):
        r = 1
        target.append(0)
        for j in range(len(config.times)):
            lf = config.L[j]*10**(config.D[j]*(1-config.F[j])/(1-config.F_MIN[j]))
            r = r * (1-math.exp(-lf*config.times[j][l[i]]/config.F[j]))

        r1 = r1 * (1-r)
        rexp.append(math.log(1-r))
    
    s = sum(rexp)
    for i in range(len(config.graph)):
        rexp[i] = 1/len(config.graph)
        rels.append(1)
    
    alloc = []
    for i in range(len(l)):
        alloc.append([])
    
    energy = 0
    exp = 1
    for i in range(len(l)):
        exp = exp - rexp[l[i]]
        r_goal = (R**(1-exp))/math.prod(rels)
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