import numpy as np
import math
import config, utils

def esrg(R, l = []):
    if(l == []):
        l, _ = utils.urv()
    
    rels = []
    alloc = []
    freq = []

    for i in range(len(l)):
        rels.append(1)
        alloc.append(-1)
        freq.append(1)

    for i in range(len(l)):
        r_goal = (R**(1 + (i+1)/len(l)))/(math.prod(rels) * R)
        # print("Rel goal:",r_goal)
        emin = float('inf')
        f1 = -1
        pr1 = -1
        dmin = float('inf')
        f2 = -1
        pr2 = -1
        for pr in range(len(config.times)):
            fr = np.arange(config.F_MIN[pr], config.F[pr], 10**(-4))
            for f in fr:
                lf = config.L[pr]*10**(config.D[pr]*(1-f)/(1-config.F_MIN[pr]))
                r = math.exp(-lf*config.times[pr][l[i]]/f)
                e = (config.PS[pr] + config.C[pr]*(f**config.A[pr]))*config.times[pr][l[i]]/f
                if(r >= r_goal and e < emin):
                    emin = e
                    f1 = f
                    pr1 = pr
                
                elif(r < r_goal and (r_goal-r) < dmin):
                    dmin = r_goal - r
                    f2 = f
                    pr2 = pr


        if(pr1 != -1): 
            freq[l[i]] = f1   
            alloc[l[i]] = pr1           
        
        else:
            freq[l[i]] = f2  
            alloc[l[i]] = pr2          
        
        lf = config.L[alloc[l[i]]]*10**(config.D[alloc[l[i]]]*(1-freq[l[i]])/(1-config.F_MIN[alloc[l[i]]]))
        r = math.exp(-lf*config.times[alloc[l[i]]][l[i]]/freq[l[i]])
        rels[l[i]] = r
    
    return alloc, freq

def efsrg(R):
    l, _ = utils.urv()
    alloc = []
    freq = []
    energy = 0
    rels = []

    for i in range(len(l)):
        alloc.append([])
        freq.append([])
        rels.append(1)
        for _ in range(len(config.times)):
            freq[i].append(1)

    for i in range(len(l)):
        r_goal = (R**(1 + (i+1)/len(l)))/(math.prod(rels) * R)
        metrics = []
        for pr in range(len(config.times)):
            f = config.F_MIN[pr]
            while(f <= 1):
                e = (config.PS[pr] + config.C[pr]*(f**config.A[pr]))*config.times[pr][l[i]]/f
                lf = config.L[pr]*10**(config.D[pr]*(1-f)/(1-config.F_MIN[pr]))
                r = math.exp(-lf*config.times[pr][l[i]]/f)
                metrics.append([r, e, pr, f])
                f +=  0.0001
        
        metrics.sort(key=lambda x: x[1])
        ll = []
        ff = []
        rr = []
        ee = []
        for j in range(len(metrics)):
            idx = -1
            if(metrics[j][2] in ll):
                idx = ll.index(metrics[j][2])
            
            if(idx >= 0):
                ll[idx] = metrics[j][2]
                ff[idx] = metrics[j][3]
                rr[idx] = metrics[j][0]
                ee[idx] = metrics[j][1]
            else:
                ll.append(metrics[j][2])
                ff.append(metrics[j][3])
                rr.append(metrics[j][0])
                ee.append(metrics[j][1])

            r = 1
            for k in rr:
                r = r * (1-k)
            
            if(1 - r >= r_goal):
                alloc[l[i]] = [ll[i] for i in range(len(ll)) ]
                rels[l[i]] = 1 - r
                energy += sum(ee)
                for pr in range(len(alloc[l[i]])):
                    freq[l[i]][alloc[l[i]][pr]] = ff[pr] 

                break

        if(alloc[l[i]] == []):
            r = 1
            for k in rr:
                r = r * (1-k)
            alloc[l[i]] = [ll[i] for i in range(len(ll)) ]
            rels[l[i]] = 1 - r
            energy += sum(ee)
            for pr in range(len(alloc[l[i]])):
                freq[l[i]][alloc[l[i]][pr]] = ff[pr] 

    return alloc, freq, math.prod(rels), energy