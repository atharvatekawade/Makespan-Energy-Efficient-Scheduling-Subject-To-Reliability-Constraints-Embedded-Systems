import math
import numpy as np
import config, utils

def ods(sp, fp, step, l = []):
    th_list = np.arange(sp, fp + step, step).tolist()
    if(l == []):
        l, _ = utils.urv()

    od = []
    idx = []
    for i in range(len(config.graph)):
        idx.append(0)
        d = 0
        for j in range(len(config.graph)):
            if(config.graph[i][j] > 0):
                d += 1
        od.append(d)

    u = [[i, od[i]] for i in range(len(config.graph))]
    
    for i in range(len(config.graph)):
        for j in range(len(config.graph)-i-1):
            if(u[j][1] < u[j + 1][1]):
                tempo = u[j]
                u[j]= u[j + 1]
                u[j + 1]= tempo
            
            elif(u[j][1] == u[j + 1][1] and l[j] < l[j+1]):
                tempo = u[j]
                u[j]= u[j + 1]
                u[j + 1]= tempo
  
    v = [u[i][0] for i in range(len(config.graph))]



    for i in range(len(v)):
        if(i not in v):
            print("Vertex missing")
            quit()

        idx[v[i]] = i

    for a in range(len(config.graph)):
        for th in th_list:
            fins = []
            avail = []
            alloc =[]

            e = 0
            r = 1

            for i in range(len(config.times)):
                avail.append(0)
            for i in range(len(l)):
                fins.append(0)
                alloc.append(-1)

            for i in range(len(l)):
                min_metric = float('inf')

                for j in range(len(config.times)):
                    s = avail[j]
                    for k in range(len(config.graph)):
                        if(config.graph[k][l[i]] > 0):
                            comm = config.graph[k][l[i]]
                            if(j == alloc[k]):
                                comm = 0

                            s = max(s, fins[k]+comm)

                    f = s + config.times[j][l[i]]
                    e1 = (config.PS[j] + config.C[j]*config.F[j]**config.A[j])*config.times[j][l[i]]/config.F[j]
                    lf = config.L[j]*10**(config.D[j]*(1-config.F[j])/(1-config.F_MIN[j]))
                    r1 = math.exp(-lf*config.times[j][l[i]]/config.F[j])

                    if(idx[l[i]] <= a):
                        metric = f + th*(1-r1)*config.times[j][l[i]]
                    else:
                        metric = e1

                    if(metric < min_metric):
                        min_metric = metric
                        alloc[l[i]] = j
                        fins[l[i]] = f

                avail[alloc[l[i]]] = fins[l[i]]
                e = e + (config.PS[alloc[l[i]]] + config.C[alloc[l[i]]]*(config.F[alloc[l[i]]]**config.A[alloc[l[i]]]))*config.times[alloc[l[i]]][l[i]]/config.F[alloc[l[i]]]
                lf = config.L[alloc[l[i]]]*10**(config.D[alloc[l[i]]]*(1-config.F[alloc[l[i]]])/(1-config.F_MIN[alloc[l[i]]]))
                r = r * math.exp(-lf*config.times[alloc[l[i]]][l[i]]/config.F[alloc[l[i]]])

            f5, r5, e5 = utils.simulate(alloc)
            config.ods_alloc.append(alloc)
            config.ods_results.append([f5, r5, e5])