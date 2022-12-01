import math
import config, utils


def mert(a, th, R, l = []):
    assert(th <= 1)
    if(l == []):
        l, _ = utils.urv()

    od = []
    idx = []

    for i in range(len(config.graph)):
        idx.append(0)
        t = 0
        for j in range(len(config.times)):
            t += config.times[j][l[i]]
        
        m = 0
        for j in range(len(config.graph)):
            if(config.graph[j][i] > 0):
                m = max(m, config.graph[j][i])

        od.append(t/len(config.times)+m)     
    
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
    assert(len(v) == len(config.graph))

    for i in range(len(config.graph)):
        if(i not in v):
            print("Vertex missing")
            quit()

        idx[v[i]] = i

    alloc = []
    avail = []
    fins = []
    rels = [1]
    energy = 0
    rexp = []

    for i in range(len(config.graph)):
        rels.append(1)
        rmax = 0
        for pr in range(len(config.times)):
            rmax = max(rmax, math.exp(-config.L[pr]*config.times[pr][i]))
        rexp.append(math.log(rmax))
    
    s = sum(rexp)
    for i in range(len(config.graph)):
        rexp[i] = rexp[i]/s
    
    for i in range(len(config.times)):
        avail.append(0)

    for i in range(len(l)):
        alloc.append(-1)
        fins.append(0)

    exp = 1
    for i in range(len(l)):
        exp = exp - rexp[l[i]]
        r_goal = (R**(1-exp))/math.prod(rels)
        metrics = []
        tmax = 0
        smax = 0
        emax = 0
        tmin = float('inf')
        smin = float('inf')
        emin = float('inf')

        for pr in range(len(config.times)):
            s = avail[pr]
            for j in range(len(config.graph)):
                if(config.graph[j][l[i]] > 0):
                    comm_time = config.graph[j][l[i]]
                    if(pr == alloc[j]):
                        comm_time = 0
                        
                    s = max(s, fins[j]+comm_time)

            r = math.exp(-config.L[pr]*config.times[pr][l[i]])
            e = (config.PS[pr] + config.C[pr]*(1**config.A[pr]))*config.times[pr][l[i]]/1

            smax = max(smax, s+config.times[pr][l[i]])
            emax = max(emax, e)
            tmax = max(tmax, config.times[pr][l[i]])

            smin = min(smin, s+config.times[pr][l[i]])
            emin = min(emin, e)
            tmin = min(tmin, config.times[pr][l[i]])
            metrics.append([s+config.times[pr][l[i]], r, e])
            
        best_metric = float('inf')
        proc = -1
        rmax = 0
        proc1 = -1

        for pr in range(len(config.times)):
            if(metrics[pr][1] >= r_goal):
                metric = th*(metrics[pr][0]-smin)/(smax-smin) + (1-th)*(config.times[pr][l[i]]-tmin)/(tmax-tmin)
                if(idx[l[i]] > a):
                    metric = metrics[pr][2]

                if(metric < best_metric):
                    proc = pr
                    best_metric = metric
            
            elif(metrics[pr][1] > rmax):
                proc1 = pr
                rmax = metrics[pr][1]

        if(proc == -1):
            proc = proc1
            print("Wrong vertex",l[i])
            quit()

        alloc[l[i]] = proc
        fins[l[i]] = metrics[proc][0]
        avail[proc] = fins[l[i]]
        energy += metrics[proc][2]
        rels[l[i]] = metrics[proc][1]
    
    config.my_alloc.append(alloc)
    config.my_results.append([max(fins), math.prod(rels), energy])

