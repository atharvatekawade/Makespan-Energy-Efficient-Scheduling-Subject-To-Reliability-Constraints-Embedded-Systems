import random
import math
import numpy as np

graph = []
times = []

PS = []
C = []
F = []
F_MIN = []
A = []
L = []
D = []

t1 = 10
t2 = 100
best_energy = 0
epsilon = 10**(-7)
ods_results = []
ods_alloc = []
my_results = []
my_alloc = []


def urv():
    x = [[i, 0] for i in range(len(graph))]

    for i in range(len(graph)-1, -1, -1):
        mt = 0
        for j in range(len(graph)):
            if(graph[i][j] > 0):
                mt = max(mt, x[j][1]+graph[i][j])
        
        rnk = 0

        for j in range(len(times)):
            rnk += times[j][i]
        
        x[i][1] = rnk/len(times) + mt

    x.sort(key=lambda x: x[1], reverse= True) 
    return [x[i][0] for i in range(len(x))], [x[i][1] for i in range(len(x))]

def simulate(alloc, freq = [], l = []):
    r = 1
    e = 0

    if(l == []):
        l, _ = urv()

    
    if(freq == []):
        for _ in range(len(graph)):
            freq.append(1)
    
    for i in range(len(l)):
        for j in range(len(l)):
            if(graph[l[i]][l[j]] > 0 and i >= j):
                print("Wrong graph levels")
                quit()

    starts = []
    fins = []
    avail = []

    for i in range(len(times)):
        avail.append(0)

    for i in range(len(l)):
        starts.append(0)
        fins.append(0)

    for i in range(len(l)):
        pr = alloc[l[i]]
        s = avail[pr]
        for j in range(len(graph)):
            if(graph[j][l[i]] > 0):
                comm = graph[j][l[i]]
                if(pr == alloc[j]):
                    comm = 0

                s = max(s, fins[j]+comm)

        f = s + times[pr][l[i]]/freq[l[i]]
        

        starts[l[i]] = s
        fins[l[i]] = f  
        avail[pr] = fins[l[i]]

        e = e + (PS[pr] + C[pr]*(freq[l[i]]**A[pr]))*times[pr][l[i]]/freq[l[i]]

        lf = L[pr]*10**(D[pr]*(1-freq[l[i]])/(1-F_MIN[pr]))
        r = r * math.exp(-lf*times[pr][l[i]]/freq[l[i]])

    
    return max(fins), r, e


def mr(l = []):
    if(l == []):
        l, _ = urv()
    
    for i in range(len(l)):
        for j in range(len(l)):
            if(graph[l[i]][l[j]] > 0 and i >= j):
                print("Wrong graph levels")
                quit()

    starts = []
    fins = []
    avail = []
    alloc = []

    for i in range(len(times)):
        avail.append(0)

    for i in range(len(l)):
        starts.append(0)
        fins.append(0)
        alloc.append(-1)

    for i in range(len(l)):
        rmax = 0
        for j in range(len(times)):
            s = avail[j]
            for k in range(len(graph)):
                if(graph[k][l[i]] > 0):
                    comm = graph[k][l[i]]
                    if(j == alloc[k]):
                        comm = 0

                    s = max(s, fins[k]+comm)

            f = s + times[j][l[i]]
            lf = L[j]*10**(D[j]*(1-F[j])/(1-F_MIN[j]))
            r1 = math.exp(-lf*times[j][l[i]]/F[j])

            if(r1 > rmax):
                rmax = r1
                alloc[l[i]] = j
                fins[l[i]] = f
 
        avail[alloc[l[i]]] = fins[l[i]]

    return alloc, max(fins)

def lec(l = []):
    if(l == []):
        l, _ = urv()
    
    for i in range(len(l)):
        for j in range(len(l)):
            if(graph[l[i]][l[j]] > 0 and i >= j):
                print("Wrong graph levels")
                quit()

    starts = []
    fins = []
    avail = []
    alloc = []

    for i in range(len(times)):
        avail.append(0)

    for i in range(len(l)):
        starts.append(0)
        fins.append(0)
        alloc.append(-1)

    for i in range(len(l)):
        emin = float('inf')
        for j in range(len(times)):
            s = avail[j]
            for k in range(len(graph)):
                if(graph[k][l[i]] > 0):
                    comm = graph[k][l[i]]
                    if(j == alloc[k]):
                        comm = 0

                    s = max(s, fins[k]+comm)

            f = s + times[j][l[i]]
            e1 = (PS[j] + C[j]*(F[j]**A[j]))*times[j][l[i]]/F[j]

            if(e1 < emin):
                emin = e1
                alloc[l[i]] = j
                fins[l[i]] = f
 
        avail[alloc[l[i]]] = fins[l[i]]

    return alloc, max(fins)


def mr_fault():
    rels = []
    alloc = []
    l, _ = urv()
    energy = 0
    for i in range(len(l)):
        procs = []
        r = 1
        for j in range(len(times)):
            procs.append(j)
            lf = L[j]*10**(D[j]*(1-F[j])/(1-F_MIN[j]))
            r1 = math.exp(-lf*times[j][l[i]]/F[j])
            r = r * (1-r1)
            energy += (PS[j] + C[j]*(F[j]**A[j]))*times[j][l[i]]/F[j]

        rels.append(1-r)
        alloc.append(procs)
    
    return alloc, math.prod(rels), energy

def soea1(alloc, R):
    y1 = []
    y2 = []

    for pr in range(len(times)):
        num = C[pr]*(A[pr]-1)*(F_MIN[pr]**(A[pr]-1)) - PS[pr]/F_MIN[pr]
        lf = L[pr]*10**(D[pr]*(1-F_MIN[pr])/(1-F_MIN[pr]))
        den = D[pr]*math.log(10)/(1-F_MIN[pr]) + 1/F_MIN[pr]
        y1.append(num/(R*lf*den))

        num = C[pr]*(A[pr]-1)*(F[pr]**(A[pr]-1)) - PS[pr]/F[pr]
        lf = L[pr]*(10**(D[pr]*(1-F[pr])/(1-F_MIN[pr])))
        den = D[pr]*math.log(10)/(1-F_MIN[pr]) + 1/F[pr]
        y2.append(num/(lf*den))

    lb = min(y1)
    ub = max(y2)
    freq = []

    freq = []
    for i in range(len(graph)):
        freq.append(1)

    while(ub - lb > epsilon):
        mid = (ub+lb)/2

        for pr in range(len(times)):
            flb = F_MIN[pr]
            fub = F[pr]
            f = (flb+fub)/2
            while(fub - flb > epsilon):
                num = C[pr]*(A[pr]-1)*(f**(A[pr]-1)) - PS[pr]/f
                lf = L[pr]*(10**(D[pr]*(1-f)/(1-F_MIN[pr])))
                den = D[pr]*math.log(10)/(1-F_MIN[pr]) + 1/f
                y = num/(R*lf*den)

                if(y < mid):
                    flb = f
                else:
                    fub = f
                f = (flb+fub)/2

            for i in range(len(graph)):
                if(alloc[i] == pr):
                    freq[i] = f   

        _, r, _ = simulate(alloc, freq=freq)

        if(r < R):
            lb = mid
        else:
            ub = mid
    
    return freq


def simulate_fault(alloc, freq = []):
    if(freq == []):
        # print("Internally supplied freq",freq)
        for i in range(len(alloc)):
            freq.append([])
            for j in range(len(times)):
                freq[i].append(1)

    l, _ = urv()
    e = 0
    rel = 1

    for i in range(len(l)):
        r = 1
        for pr in alloc[l[i]]:
            # print(pr, l[i], len(freq[l[i]]), len(times))
            # print(freq[l[i]][pr])
            e = e + (PS[pr] + C[pr]*(freq[l[i]][pr]**A[pr]))*times[pr][l[i]]/freq[l[i]][pr]
            lf = L[pr]*10**(D[pr]*(1-freq[l[i]][pr])/(1-F_MIN[pr]))
            r = r * (1-math.exp(-lf*times[pr][l[i]]/freq[l[i]][pr]))

        rel = rel * (1-r)
    
    return rel, e

def mert_fault(R, l = []):
    if(l == []):
        l, _ = urv()
    
    rels = []
    rexp = []
    target = []
    r1 = 1
    for i in range(len(l)):
        r = 1
        target.append(0)
        for j in range(len(times)):
            lf = L[j]*10**(D[j]*(1-F[j])/(1-F_MIN[j]))
            r = r * (1-math.exp(-lf*times[j][l[i]]/F[j]))

        r1 = r1 * (1-r)
        rexp.append(math.log(1-r))
    
    s = sum(rexp)
    for i in range(len(graph)):
        rexp[i] = 1/len(graph)
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
        for pr in range(len(times)):
            e = (PS[pr] + C[pr]*(F[pr]**A[pr]))*times[pr][l[i]]/F[pr]
            lf = L[pr]*10**(D[pr]*(1-F[pr])/(1-F_MIN[pr]))
            r = 1-math.exp(-lf*times[pr][l[i]]/F[pr])
            metrics.append([r, e, pr])
        
        metrics.sort(key=lambda x: x[1])
        r = 1
        for k in range(len(metrics)):
            if(1-r >= r_goal):
                target[l[i]] = 1
                break

            r = r*metrics[k][0]
            energy += metrics[k][1]
            alloc[l[i]].append(metrics[k][2])

        rels[l[i]] = 1-r
    
    rel, energy = simulate_fault(alloc)

    while(rel > R):
        items = []
        for i in range(len(alloc)):
            if(len(alloc[i]) > 1):
                for pr in alloc[i]:
                    new_alloc = [row[:] for row in alloc]
                    new_alloc[i].remove(pr)
                    rel, _ = simulate_fault(new_alloc)
                    e = (PS[pr] + C[pr]*(F[pr]**A[pr]))*times[pr][i]/F[pr]
                    if(rel > R):
                        items.append([i, pr, e])

        if(len(items) == 0):
            break

        emax = 0
        idx = -1
        for i in range(len(items)):
            if(items[i][2] > emax):
                idx = i
                emax = items[i][2]

        alloc[items[idx][0]].remove(items[idx][1])
        rel, _ = simulate_fault(alloc)
    
    rel, energy = simulate_fault(alloc)

    return alloc, rel, energy


def mert_fault1(R, l = []):
    if(l == []):
        l, _ = urv()
    
    rels = []
    rexp = []
    target = []
    r1 = 1
    for i in range(len(l)):
        r = 1
        target.append(0)
        for j in range(len(times)):
            lf = L[j]*10**(D[j]*(1-F[j])/(1-F_MIN[j]))
            r = r * (1-math.exp(-lf*times[j][l[i]]/F[j]))

        r1 = r1 * (1-r)
        rexp.append(math.log(1-r))
    
    s = sum(rexp)
    for i in range(len(graph)):
        rexp[i] = 1/len(graph)
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
        for pr in range(len(times)):
            e = (PS[pr] + C[pr]*(F[pr]**A[pr]))*times[pr][l[i]]/F[pr]
            lf = L[pr]*10**(D[pr]*(1-F[pr])/(1-F_MIN[pr]))
            r = 1-math.exp(-lf*times[pr][l[i]]/F[pr])
            metrics.append([r, e, pr])
        
        metrics.sort(key=lambda x: x[1])
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

def mert_fault2(R, l = []):
    if(l == []):
        l, _ = urv()
    
    rels = []
    rexp = []
    target = []
    r1 = 1
    for i in range(len(l)):
        r = 1
        target.append(0)
        for j in range(len(times)):
            lf = L[j]*10**(D[j]*(1-F[j])/(1-F_MIN[j]))
            r = r * (1-math.exp(-lf*times[j][l[i]]/F[j]))

        r1 = r1 * (1-r)
        rexp.append(math.log(1-r))
    
    s = sum(rexp)
    for i in range(len(graph)):
        rexp[i] = 1/len(graph)
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
        rmin = float('inf')
        rmax = 0
        emin = float('inf')
        emax = 0

        for pr in range(len(times)):
            e = (PS[pr] + C[pr]*(F[pr]**A[pr]))*times[pr][l[i]]/F[pr]
            lf = L[pr]*10**(D[pr]*(1-F[pr])/(1-F_MIN[pr]))
            r = 1-math.exp(-lf*times[pr][l[i]]/F[pr])
            metrics.append([r, e, pr])
            rmin = min(r, rmin)
            rmax = max(r, rmax)
            emin = min(e, emin)
            emax = max(e, emax)
        
        for pr in range(len(metrics)):
            qty = 0.8*(metrics[pr][1] - emin)/(emax - emin) + 0.2*(rmax - metrics[pr][0])/(rmax - rmin)
            metrics[pr].append(qty)
        
        metrics.sort(key=lambda x: x[3])
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

def rr(R, l = []):
    if(l == []):
        l, _ = urv()
    
    rels = []
    rexp = []
    target = []
    r1 = 1
    for i in range(len(l)):
        r = 1
        target.append(0)
        for j in range(len(times)):
            lf = L[j]*10**(D[j]*(1-F[j])/(1-F_MIN[j]))
            r = r * (1-math.exp(-lf*times[j][l[i]]/F[j]))

        r1 = r1 * (1-r)
        rexp.append(math.log(1-r))
    
    s = sum(rexp)
    for i in range(len(graph)):
        rexp[i] = 1/len(graph)
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
        for pr in range(len(times)):
            e = (PS[pr] + C[pr]*(F[pr]**A[pr]))*times[pr][l[i]]/F[pr]
            lf = L[pr]*10**(D[pr]*(1-F[pr])/(1-F_MIN[pr]))
            r = 1-math.exp(-lf*times[pr][l[i]]/F[pr])
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

def max_re(R, l = []):
    if(l == []):
        l, _ = urv()
    
    rels = []
    # rexp = []
    target = []
    # r1 = 1
    for i in range(len(l)):
        r = 1
        target.append(0)
        rels.append(1)

    #     for j in range(len(times)):
    #         lf = L[j]*10**(D[j]*(1-F[j])/(1-F_MIN[j]))
    #         r = r * (1-math.exp(-lf*times[j][l[i]]/F[j]))

    #     r1 = r1 * (1-r)
    #     rexp.append(math.log(1-r))
    
    # s = sum(rexp)
    # for i in range(len(graph)):
    #     rexp[i] = rexp[i]/s
    #     rels.append(1)
    
    alloc = []
    for i in range(len(l)):
        alloc.append([])
    
    energy = 0
    exp = 1
    for i in range(len(l)):
        # exp = exp - rexp[l[i]]
        r_goal = R**(1/len(graph))
        metrics = []
        for pr in range(len(times)):
            e = (PS[pr] + C[pr]*(F[pr]**A[pr]))*times[pr][l[i]]/F[pr]
            lf = L[pr]*10**(D[pr]*(1-F[pr])/(1-F_MIN[pr]))
            r = 1-math.exp(-lf*times[pr][l[i]]/F[pr])
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
        if(rels[l[i]] < r_goal):
            print("Failed to achieve target",l[i], rels[l[i]], r_goal)

    return alloc, math.prod(rels), energy


def find(freq, task, pr, rt):

    f = 0.01
    lf = L[pr]*10**(D[pr]*(1-f)/(1-F_MIN[pr]))
    r2 = math.exp(-lf*times[pr][task]/f)

    f = 0.02
    lf = L[pr]*10**(D[pr]*(1-f)/(1-F_MIN[pr]))
    r3 = math.exp(-lf*times[pr][task]/f)

    f = freq[task][pr]
    lf = L[pr]*10**(D[pr]*(1-f)/(1-F_MIN[pr]))
    r = math.exp(-lf*times[pr][task]/f)
    r1 = r

    epsilon = 10**(-5)
    while(r-rt > epsilon and f >= F_MIN[pr]):
        lf = L[pr]*10**(D[pr]*(1-f)/(1-F_MIN[pr]))
        r = math.exp(-lf*times[pr][task]/f)
        f = f -  0.0001

    # print(f"Original rel: {r1} Obtained rel: {r} Target rel: {rt} Middle Rel: {r3} Min rel: {r2} Original freq: {freq[task][pr]} Freq: {f}")
    return f

def adjust(alloc, R):
    freq = []
    for i in range(len(alloc)):
        freq.append([])
        for j in range(len(times)):
            freq[i].append(1)
    
    r, _ = simulate_fault(alloc, freq)
    epsilon = 10**(-9)

    while(abs(R-r) > epsilon):
        task = -1
        proc = -1
        fr = -1
        min_qty = float('inf')

        for i in range(len(alloc)): 
            for j in alloc[i]:
                r1 = 1
                r2 = 1
                for k in alloc[i]:
                    lf = L[k]*10**(D[k]*(1-freq[i][k])/(1-F_MIN[k]))
                    r1 = r1 * (1-math.exp(-lf*times[k][i]/freq[i][k]))

                    if(k != j):
                        r2 = r2 * (1-math.exp(-lf*times[k][i]/freq[i][k]))
                
                rt = 1-(1-(1-r1)*R/r)/r2
                e1 = (PS[j] + C[j]*(freq[i][j]**A[j]))*times[j][i]/freq[i][j]
                fq = find(freq, i, j, rt)
                e2 = (PS[j] + C[j]*(fq**A[j]))*times[j][i]/fq

                if(e2 - e1 < min_qty):
                    task = i
                    proc = j
                    fr = fq
                    min_qty = e2 - e1  

        # print("Energy gain:",min_qty)
        fnew = (freq[task][proc]+fr)/2
        if(abs(fnew - freq[task][proc]) >=  0.0001):
            freq[task][proc] = fnew
        else:
            freq[task][proc] = fr
        r, _ = simulate_fault(alloc, freq)

        if(min_qty == 0):
            break
    
    return freq

def esrg(R, l = []):
    if(l == []):
        l, _ = urv()
    
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
        for pr in range(len(times)):
            fr = np.arange(F_MIN[pr], F[pr], 10**(-4))
            for f in fr:
                lf = L[pr]*10**(D[pr]*(1-f)/(1-F_MIN[pr]))
                r = math.exp(-lf*times[pr][l[i]]/f)
                e = (PS[pr] + C[pr]*(f**A[pr]))*times[pr][l[i]]/f
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
        
        lf = L[alloc[l[i]]]*10**(D[alloc[l[i]]]*(1-freq[l[i]])/(1-F_MIN[alloc[l[i]]]))
        r = math.exp(-lf*times[alloc[l[i]]][l[i]]/freq[l[i]])
        rels[l[i]] = r
    
    return alloc, freq


def ods(l = []):
    if(l == []):
        l, _ = urv()

    od = []
    idx = []
    for i in range(len(graph)):
        idx.append(0)
        d = 0
        for j in range(len(graph)):
            if(graph[i][j] > 0):
                d += 1
        od.append(d)

    u = [[i, od[i]] for i in range(len(graph))]
    
    for i in range(len(graph)):
        for j in range(len(graph)-i-1):
            if(u[j][1] < u[j + 1][1]):
                tempo = u[j]
                u[j]= u[j + 1]
                u[j + 1]= tempo
            
            elif(u[j][1] == u[j + 1][1] and l[j] < l[j+1]):
                tempo = u[j]
                u[j]= u[j + 1]
                u[j + 1]= tempo
  
    v = [u[i][0] for i in range(len(graph))]

    final = []
    fin = float('inf')
    rel = 0
    ene = float('inf')
    min_qty = float('inf')


    for i in range(len(v)):
        if(i not in v):
            print("Vertex missing")
            quit()

        idx[v[i]] = i

    for a in range(len(graph)):
        for th in range(11):
            fins = []
            avail = []
            alloc =[]

            e = 0
            r = 1

            for i in range(len(times)):
                avail.append(0)
            for i in range(len(l)):
                fins.append(0)
                alloc.append(-1)

            for i in range(len(l)):
                min_metric = float('inf')

                for j in range(len(times)):
                    s = avail[j]
                    for k in range(len(graph)):
                        if(graph[k][l[i]] > 0):
                            comm = graph[k][l[i]]
                            if(j == alloc[k]):
                                comm = 0

                            s = max(s, fins[k]+comm)

                    f = s + times[j][l[i]]
                    e1 = (PS[j] + C[j]*F[j]**A[j])*times[j][l[i]]/F[j]
                    lf = L[j]*10**(D[j]*(1-F[j])/(1-F_MIN[j]))
                    r1 = math.exp(-lf*times[j][l[i]]/F[j])

                    if(idx[l[i]] <= a):
                        metric = f + th*100*(1-r1)*times[j][l[i]]
                    else:
                        metric = e1

                    if(metric < min_metric):
                        min_metric = metric
                        alloc[l[i]] = j
                        fins[l[i]] = f

                avail[alloc[l[i]]] = fins[l[i]]
                e = e + (PS[alloc[l[i]]] + C[alloc[l[i]]]*(F[alloc[l[i]]]**A[alloc[l[i]]]))*times[alloc[l[i]]][l[i]]/F[alloc[l[i]]]
                lf = L[alloc[l[i]]]*10**(D[alloc[l[i]]]*(1-F[alloc[l[i]]])/(1-F_MIN[alloc[l[i]]]))
                r = r * math.exp(-lf*times[alloc[l[i]]][l[i]]/F[alloc[l[i]]])

            f5, r5, e5 = simulate(alloc)
            ods_alloc.append(alloc)
            ods_results.append([f5, r5, e5])


def my(a, th, R, l = []):
    assert(th <= 1)
    if(l == []):
        l, _ = urv()

    od = []
    idx = []

    for i in range(len(graph)):
        idx.append(0)
        t = 0
        for j in range(len(times)):
            t += times[j][l[i]]
        
        m = 0
        for j in range(len(graph)):
            if(graph[j][i] > 0):
                m = max(m, graph[j][i])

        od.append(t/len(times)+m)     
    
    u = [[i, od[i]] for i in range(len(graph))]
    
    for i in range(len(graph)):
        for j in range(len(graph)-i-1):
            if(u[j][1] < u[j + 1][1]):
                tempo = u[j]
                u[j]= u[j + 1]
                u[j + 1]= tempo
            
            elif(u[j][1] == u[j + 1][1] and l[j] < l[j+1]):
                tempo = u[j]
                u[j]= u[j + 1]
                u[j + 1]= tempo
  
    v = [u[i][0] for i in range(len(graph))]
    assert(len(v) == len(graph))

    for i in range(len(graph)):
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

    for i in range(len(graph)):
        rels.append(1)
        rmax = 0
        for pr in range(len(times)):
            rmax = max(rmax, math.exp(-L[pr]*times[pr][i]))
        rexp.append(math.log(rmax))
    
    s = sum(rexp)
    for i in range(len(graph)):
        rexp[i] = rexp[i]/s
    
    for i in range(len(times)):
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

        for pr in range(len(times)):
            s = avail[pr]
            for j in range(len(graph)):
                if(graph[j][l[i]] > 0):
                    comm_time = graph[j][l[i]]
                    if(pr == alloc[j]):
                        comm_time = 0
                        
                    s = max(s, fins[j]+comm_time)

            r = math.exp(-L[pr]*times[pr][l[i]])
            e = (PS[pr] + C[pr]*(1**A[pr]))*times[pr][l[i]]/1

            smax = max(smax, s+times[pr][l[i]])
            emax = max(emax, e)
            tmax = max(tmax, times[pr][l[i]])

            smin = min(smin, s+times[pr][l[i]])
            emin = min(emin, e)
            tmin = min(tmin, times[pr][l[i]])
            metrics.append([s+times[pr][l[i]], r, e])
            
        best_metric = float('inf')
        proc = -1
        rmax = 0
        proc1 = -1

        for pr in range(len(times)):
            if(metrics[pr][1] >= r_goal):
                metric = th*(metrics[pr][0]-smin)/(smax-smin) + (1-th)*(times[pr][l[i]]-tmin)/(tmax-tmin)
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
    
    my_alloc.append(alloc)
    my_results.append([max(fins), math.prod(rels), energy])
    # return alloc, max(fins), math.prod(rels), energy

def rel_constraint1(R, l = []):
    if(l == []):
        l, _ = urv()

    alloc = []
    avail = []
    fins = []
    rels = [1]
    energy = 0
    rexp = []

    for i in range(len(graph)):
        rels.append(1)
        rmax = 0
        for pr in range(len(times)):
            rmax = max(rmax, math.exp(-L[pr]*times[pr][i]))
        rexp.append(math.log(rmax))
    
    s = sum(rexp)
    for i in range(len(graph)):
        rexp[i] = rexp[i]/s
    
    for i in range(len(times)):
        avail.append(0)

    for i in range(len(graph)):
        alloc.append(-1)
        fins.append(0)

    exp = 1
    for i in range(len(l)):
        exp -= rexp[l[i]]
        r_goal = R**(1-exp)/math.prod(rels)
        metrics = []
        smax = 0
        emax = 0
        for pr in range(len(times)):
            s = avail[pr]
            for j in range(len(graph)):
                if(graph[j][l[i]] > 0):
                    if(pr == alloc[j]):
                        comm_time = 0
                    else:
                        comm_time = graph[j][l[i]]

                    s = max(s, fins[j]+comm_time)

            r = math.exp(-L[pr]*times[pr][l[i]])
            e = (PS[pr] + C[pr]*(1**A[pr]))*times[pr][l[i]]/1

            smax = max(smax, s+times[pr][l[i]])
            emax = max(emax, e)
            metrics.append([s+times[pr][l[i]], r, e])
            
        
        best_metric = float('inf')
        proc = -1
        rmax = 0
        proc1 = -1

        for pr in range(len(times)):
            if(metrics[pr][1] >= r_goal):
                metric = 0.7*metrics[pr][0]/smax + 0.3*metrics[pr][2]/emax
                if(metric < best_metric):
                    proc = pr
                    best_metric = metric
            
            elif(metrics[pr][1] > rmax):
                proc1 = pr
                rmax = metrics[pr][1]

        if(proc == -1):
            proc = proc1

        alloc[l[i]] = proc
        fins[l[i]] = metrics[proc][0]
        avail[proc] = fins[l[i]]
        energy += metrics[proc][2]
        rels[l[i]] = metrics[proc][1]
    
    return alloc, max(fins), math.prod(rels), energy


def mecrg(R, l = []):
    if(l == []):
        l, _ = urv()
    
    alloc = []
    avail = []
    rels = []
    rmax = []

    for i in range(len(graph)):
        alloc.append(-1)
        rels.append(1)
        r = 0
        for pr in range(len(times)):
            r = max(r, math.exp(-L[pr]*times[pr][i]))
        rmax.append(r)
    
    for i in range(len(times)):
        avail.append(0)
    
    max_rel = math.prod(rmax)

    for i in range(len(l)):
        max_rel /= rmax[l[i]]
        r_goal = R/(math.prod(rels)*max_rel)
        emin = float('inf')
        pmin = -1
        for pr in range(len(times)):
            if(math.exp(-L[pr]*times[pr][l[i]]) >= r_goal):
                e = (PS[pr] + C[pr]*(1**A[pr]))*times[pr][l[i]]/F[pr]
                if(e < emin):
                    pmin = pr
                    emin = e
        
        if(pmin == -1):
            print("Wrong processor rel constraint violated")
            quit()
        
        else:
            rels[l[i]] = math.exp(-L[pmin]*times[pmin][l[i]])
            alloc[l[i]] = pmin
    
    return alloc


def efsrg(R):
    l, _ = urv()
    alloc = []
    freq = []
    energy = 0
    rels = []

    for i in range(len(l)):
        alloc.append([])
        freq.append([])
        rels.append(1)
        for _ in range(len(times)):
            freq[i].append(1)

    for i in range(len(l)):
        r_goal = (R**(1 + (i+1)/len(l)))/(math.prod(rels) * R)
        metrics = []
        for pr in range(len(times)):
            f = F_MIN[pr]
            while(f <= 1):
                e = (PS[pr] + C[pr]*(f**A[pr]))*times[pr][l[i]]/f
                lf = L[pr]*10**(D[pr]*(1-f)/(1-F_MIN[pr]))
                r = math.exp(-lf*times[pr][l[i]]/f)
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


rels = [(1000+i)/1000 for i in range(1, 21)]
results = [[0, 0, 0, 0, 0] for i in range(len(rels))]
# results2 = [[0, 0] for i in range(len(rels))]
itr = int(input("Iterations: "))
max_val = 0

for g in range(itr):
    print("Itr:",g+1)
    graph = []
    n = 20
    v = int(n*(n+1)/2) - 1

    for i in range(v):
        graph.append([])
        for _ in range(v):
            graph[i].append(0)

    curr_level = []
    for i in range(1, n):
        curr_level.append(i)
        graph[0][i] = random.randint(t1, t2)

    for i in curr_level:
        graph[i][i+n] = random.randint(t1, t2)
        
    for i in range(n+1, 2*n):
        graph[n][i] = random.randint(t1, t2)

    curr_level = [curr_level[i]+n-1 for i in range(len(curr_level))]

    l = len(curr_level)
    while(l > 2):
        for i in range(1, len(curr_level)):
            graph[curr_level[i]][curr_level[i]+l-1] = random.randint(t1, t2)

        curr_level = [curr_level[i]+l-1 for i in range(len(curr_level))][1:]

        l -= 1
        if(l > 2):
            for i in range(1, len(curr_level)):
                graph[curr_level[0]][curr_level[i]] = random.randint(t1, t2)

        assert(len(curr_level) == l)

    # n = 5
    # v = (2+n)*2**n - 1

    # for i in range(v):
    #     graph.append([])
    #     for _ in range(v):
    #         graph[i].append(0)

    # curr_level = [0]
    # next_level = []

    # for _ in range(n):
    #     next_level = []
    #     for i in range(len(curr_level)):
    #         graph[curr_level[i]][2*curr_level[i]+1] = random.randint(t1, t2)
    #         next_level.append(2*curr_level[i]+1)

    #         graph[curr_level[i]][2*curr_level[i]+2] = random.randint(t1, t2)
    #         next_level.append(2*curr_level[i]+2)
        
    #     curr_level = next_level[:]
    
    # for i in range(n):
    #     next_level = []
    #     turn = []
    #     curr = 0
    #     ctr = 0
    #     for j in range(len(curr_level)):
    #         next_level.append(curr_level[j]+2**n)
    #         turn.append(curr)
    #         ctr += 1
    #         if(ctr == 2**i):
    #             ctr = 0
    #             curr = 1 - curr
        
    #     for j in range(len(curr_level)):
    #         graph[curr_level[j]][next_level[j]] = random.randint(t1, t2)
    #         if(turn[j] == 0):
    #             graph[curr_level[j]][next_level[j]+2**i] = random.randint(t1, t2)
    #         else:
    #             graph[curr_level[j]][next_level[j]-2**i] = random.randint(t1, t2)

    #     curr_level = next_level[:]

    PS = []
    C = []
    F = []
    F_MIN = []
    A = []
    L = []
    D = []

    pr = 32
    times = []
    for i in range(pr):
        times.append([])
        for _ in range(len(graph)):
            times[i].append(random.randint(t1, t2))
                

    for i in range(len(times)):
        PS.append(random.randint(400, 800)/1000)
        C.append(random.randint(800, 1300)/1000)
        F.append(1)
        F_MIN.append(random.randint(300, 700)/1000)
        A.append(random.randint(2700,3000)/1000)
        L.append(random.randint(100, 1000)/10**8)
        D.append(random.randint(1000, 3000)/1000)


    # ods_results = []
    # ods_alloc = []
    # print("Running ODS")
    # ods()
    # assert(len(ods_results) == len(ods_alloc) == 11*len(graph))

    alloc, _ = mr()
    _, max_rel1, _ = simulate(alloc)
    _, max_rel2, max_energy = mr_fault()
    ratio = max_rel2/max_rel1
    print("Ratio:", ratio)

    if(ratio >= 1.0118 and ratio < 1.012):
        for i in range(len(rels)):
            R = rels[i]*max_rel1
            print("Rel constraint:",rels[i], R)
            if(R > max_rel2):
                break

            # alloc, _ = mr()
            # freq = soea1(alloc, R)
            # f1, r1, e1 = simulate(alloc, freq = freq)

            # results1[i][0] += f1
            # results2[i][0] += e1

            # print(f"Algo: MR No. of vertices: {len(graph)} Span: {f1} Rel: {r1} Energy: {e1}")

            # alloc, _ = lec()
            # f1, r1, e5 = simulate(alloc)
            # if(R - r1 > 10**(-5)):
            #     alloc, _ = mr()
            #     freq = soea1(alloc, R)
            #     f1, r1, e5 = simulate(alloc, freq = freq)

            # print(f"Algo: LEC No. of vertices: {len(graph)} Span: {f1} Rel: {r1} Energy: {e5}")

            # alloc, freq = esrg(R)
            # f1, r1, e2 = simulate(alloc, freq = freq)
            # if(R - r1 > 10**(-5)):
            #     alloc, _ = mr()
            #     freq = soea1(alloc, R)
            #     f1, r1, e2 = simulate(alloc, freq = freq)

            # results1[i][1] += f1
            # results2[i][1] += e2

            # print(f"Algo: ESRG No. of vertices: {len(graph)} Span: {f1} Rel: {r1} Energy: {e2}")

            # alloc = mecrg(R)
            # freq = soea1(alloc, R)
            # f1, r1, e3 = simulate(alloc, freq = freq)
            # results1[i][2] += f1
            # results2[i][2] += e3
            # print(f"Algo: MECRG No. of vertices: {len(graph)} Span: {f1} Rel: {r1} Energy: {e3}")

            # E = min(e1, e2, e3)

            # f2 = float('inf')
            # r2 = 0
            # e2 = float('inf')

            # f3 = float('inf')
            # r3 = 0
            # e3 = float('inf')

            # for j in range(len(ods_results)):
            #     if(ods_results[j][1] >= R):
            #         freq = soea1(ods_alloc[j], R)
            #         f, r, e = simulate(ods_alloc[j], freq = freq)
            #         if(e <= E):
            #             if(f < f2):
            #                 f2 = f
            #                 r2 = r
            #                 e2 = e
                        
            #             elif(f == f2 and e < e2):
            #                 f2 = f
            #                 r2 = r
            #                 e2 = e
                    
            #         else:
            #             if(e < e3):
            #                 f3 = f
            #                 r3 = r
            #                 e3 = e
                        
            #             elif(e == e3 and f < f3):
            #                 f3 = f
            #                 r3 = r
            #                 e3 = e

            # if(f2 == float('inf')):
            #     if(f3 ==float('inf')):
            #         alloc, _ = mr()
            #         freq = soea1(alloc, R)
            #         f2, r2, e2 = simulate(alloc, freq = freq)
                
            #     else:
            #         f2 = f3
            #         e2 = e3
            #         r2 = r3

            # results1[i][3] += f2
            # results2[i][3] += e2

            # print(f"Algo: ODS No. of vertices: {len(graph)} Span: {f2} Rel: {r2} Energy: {e2}")

            # my_results = []
            # my_alloc = []
            # for a in range(-1, len(graph)):
            #     for th in range(11):
            #         my(a = a, th = th/10, R = R)

            # assert(len(my_results) == len(my_alloc))

            # f2 = float('inf')
            # r2 = 0
            # e2 = float('inf')

            # f3 = float('inf')
            # r3 = 0
            # e3 = float('inf')

            # for j in range(len(my_results)):
            #     if(my_results[j][1] >= R):
            #         freq = soea1(my_alloc[j], R)
            #         f, r, e = simulate(my_alloc[j], freq = freq)
                    
            #         if(e <= E):
            #             if(f < f2):
            #                 f2 = f
            #                 r2 = r
            #                 e2 = e
                        
            #             elif(f == f2 and e < e2):
            #                 f2 = f
            #                 r2 = r
            #                 e2 = e
                    
            #         else:
            #             if(e < e3):
            #                 f3 = f
            #                 r3 = r
            #                 e3 = e
                        
            #             elif(e == e3 and f < f3):
            #                 f3 = f
            #                 r3 = r
            #                 e3 = e
            #     else:
            #         print("Something went wrong in results of mert")
            #         quit()

            # if(f2 == float('inf')):
            #     f2 = f3
            #     e2 = e3
            #     r2 = r3

            # results1[i][4] += f2
            # results2[i][4] += e2

            # print(f"Algo: MERT No. of vertices: {len(graph)} Span: {f2} Rel: {r2} Energy: {e2}")
            # print("\n")

            # alloc, _, _ = max_re(R)
            # freq = adjust(alloc, R)
            # rel, energy = simulate_fault(alloc, freq)
            # results[i][0] += energy
            # print(f"Graph: {len(graph)} Algo: Max-Re fault tolerance: Rel:{rel} Energy: {energy}")

            print(f"Graph: {len(graph)} Algo: Max-Re fault tolerance: Rel:{max_rel2} Energy: {max_energy}")

            alloc, _, _ = rr(R)
            freq = adjust(alloc, R)
            rel, energy = simulate_fault(alloc, freq)
            results[i][1] += energy
            print(f"Graph: {len(graph)} Algo: RR fault tolerance: Rel:{rel} Energy: {energy}")


            alloc, freq, rel, energy = efsrg(R)
            results[i][2] += energy
            print(f"Graph: {len(graph)} Algo: ESFRG fault tolerance: Rel:{rel} Energy: {energy}")


            alloc, _, _ = mert_fault1(R)
            freq = adjust(alloc, R)
            rel, energy = simulate_fault(alloc, freq)
            results[i][3] += energy
            print(f"Graph: {len(graph)} Algo: EAFTS fault tolerance: Rel:{rel} Energy: {energy}")

            alloc, _, _ = mert_fault2(R)
            freq = adjust(alloc, R)
            rel, energy = simulate_fault(alloc, freq)
            results[i][4] += energy
            print(f"Graph: {len(graph)} Algo: EAFTS1 fault tolerance: Rel:{rel} Energy: {energy}")

            print("\n")
        
        break

# for i in range(len(results)):
#     for j in range(len(results[i])):
#         results[i][j] = results[i][j]/itr

print("Energy:",results)



