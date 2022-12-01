import numpy as np
import math
import config


def urv():
    x = [[i, 0] for i in range(len(config.graph))]

    for i in range(len(config.graph)-1, -1, -1):
        mt = 0
        for j in range(len(config.graph)):
            if(config.graph[i][j] > 0):
                mt = max(mt, x[j][1]+config.graph[i][j])
        
        rnk = 0

        for j in range(len(config.times)):
            rnk += config.times[j][i]
        
        x[i][1] = rnk/len(config.times) + mt

    x.sort(key=lambda x: x[1], reverse= True) 
    return [x[i][0] for i in range(len(x))], [x[i][1] for i in range(len(x))]

def simulate(alloc, freq = [], l = []):
    r = 1
    e = 0

    if(l == []):
        l, _ = urv()

    
    if(freq == []):
        for _ in range(len(config.graph)):
            freq.append(1)
    
    for i in range(len(l)):
        for j in range(len(l)):
            if(config.graph[l[i]][l[j]] > 0 and i >= j):
                print("Wrong config.graph levels")
                quit()

    starts = []
    fins = []
    avail = []

    for i in range(len(config.times)):
        avail.append(0)

    for i in range(len(l)):
        starts.append(0)
        fins.append(0)

    for i in range(len(l)):
        pr = alloc[l[i]]
        s = avail[pr]
        for j in range(len(config.graph)):
            if(config.graph[j][l[i]] > 0):
                comm = config.graph[j][l[i]]
                if(pr == alloc[j]):
                    comm = 0

                s = max(s, fins[j]+comm)

        f = s + config.times[pr][l[i]]/freq[l[i]]
        

        starts[l[i]] = s
        fins[l[i]] = f  
        avail[pr] = fins[l[i]]

        e = e + (config.PS[pr] + config.C[pr]*(freq[l[i]]**config.A[pr]))*config.times[pr][l[i]]/freq[l[i]]

        lf = config.L[pr]*10**(config.D[pr]*(1-freq[l[i]])/(1-config.F_MIN[pr]))
        r = r * math.exp(-lf*config.times[pr][l[i]]/freq[l[i]])

    
    return max(fins), r, e


def mr(l = []):
    if(l == []):
        l, _ = urv()
    
    for i in range(len(l)):
        for j in range(len(l)):
            if(config.graph[l[i]][l[j]] > 0 and i >= j):
                print("Wrong config.graph levels")
                quit()

    starts = []
    fins = []
    avail = []
    alloc = []

    for i in range(len(config.times)):
        avail.append(0)

    for i in range(len(l)):
        starts.append(0)
        fins.append(0)
        alloc.append(-1)

    for i in range(len(l)):
        rmax = 0
        for j in range(len(config.times)):
            s = avail[j]
            for k in range(len(config.graph)):
                if(config.graph[k][l[i]] > 0):
                    comm = config.graph[k][l[i]]
                    if(j == alloc[k]):
                        comm = 0

                    s = max(s, fins[k]+comm)

            f = s + config.times[j][l[i]]
            lf = config.L[j]*10**(config.D[j]*(1-config.F[j])/(1-config.F_MIN[j]))
            r1 = math.exp(-lf*config.times[j][l[i]]/config.F[j])

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
            if(config.graph[l[i]][l[j]] > 0 and i >= j):
                print("Wrong config.graph levels")
                quit()

    starts = []
    fins = []
    avail = []
    alloc = []

    for i in range(len(config.times)):
        avail.append(0)

    for i in range(len(l)):
        starts.append(0)
        fins.append(0)
        alloc.append(-1)

    for i in range(len(l)):
        emin = float('inf')
        for j in range(len(config.times)):
            s = avail[j]
            for k in range(len(config.graph)):
                if(config.graph[k][l[i]] > 0):
                    comm = config.graph[k][l[i]]
                    if(j == alloc[k]):
                        comm = 0

                    s = max(s, fins[k]+comm)

            f = s + config.times[j][l[i]]
            e1 = (config.PS[j] + config.C[j]*(config.F[j]**config.A[j]))*config.times[j][l[i]]/config.F[j]

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
        for j in range(len(config.times)):
            procs.append(j)
            lf = config.L[j]*10**(config.D[j]*(1-config.F[j])/(1-config.F_MIN[j]))
            r1 = math.exp(-lf*config.times[j][l[i]]/config.F[j])
            r = r * (1-r1)
            energy += (config.PS[j] + config.C[j]*(config.F[j]**config.A[j]))*config.times[j][l[i]]/config.F[j]

        rels.append(1-r)
        alloc.append(procs)
    
    return alloc, math.prod(rels), energy

def soea1(alloc, R):
    y1 = []
    y2 = []

    for pr in range(len(config.times)):
        num = config.C[pr]*(config.A[pr]-1)*(config.F_MIN[pr]**(config.A[pr]-1)) - config.PS[pr]/config.F_MIN[pr]
        lf = config.L[pr]*10**(config.D[pr]*(1-config.F_MIN[pr])/(1-config.F_MIN[pr]))
        den = config.D[pr]*math.log(10)/(1-config.F_MIN[pr]) + 1/config.F_MIN[pr]
        y1.append(num/(R*lf*den))

        num = config.C[pr]*(config.A[pr]-1)*(config.F[pr]**(config.A[pr]-1)) - config.PS[pr]/config.F[pr]
        lf = config.L[pr]*(10**(config.D[pr]*(1-config.F[pr])/(1-config.F_MIN[pr])))
        den = config.D[pr]*math.log(10)/(1-config.F_MIN[pr]) + 1/config.F[pr]
        y2.append(num/(lf*den))

    lb = min(y1)
    ub = max(y2)
    freq = []

    freq = []
    for i in range(len(config.graph)):
        freq.append(1)

    while(ub - lb > config.epsilon):
        mid = (ub+lb)/2

        for pr in range(len(config.times)):
            flb = config.F_MIN[pr]
            fub = config.F[pr]
            f = (flb+fub)/2
            while(fub - flb > config.epsilon):
                num = config.C[pr]*(config.A[pr]-1)*(f**(config.A[pr]-1)) - config.PS[pr]/f
                lf = config.L[pr]*(10**(config.D[pr]*(1-f)/(1-config.F_MIN[pr])))
                den = config.D[pr]*math.log(10)/(1-config.F_MIN[pr]) + 1/f
                y = num/(R*lf*den)

                if(y < mid):
                    flb = f
                else:
                    fub = f
                f = (flb+fub)/2

            for i in range(len(config.graph)):
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
            for j in range(len(config.times)):
                freq[i].append(1)

    l, _ = urv()
    e = 0
    rel = 1

    for i in range(len(l)):
        r = 1
        for pr in alloc[l[i]]:
            # print(pr, l[i], len(freq[l[i]]), len(config.times))
            # print(freq[l[i]][pr])
            e = e + (config.PS[pr] + config.C[pr]*(freq[l[i]][pr]**config.A[pr]))*config.times[pr][l[i]]/freq[l[i]][pr]
            lf = config.L[pr]*10**(config.D[pr]*(1-freq[l[i]][pr])/(1-config.F_MIN[pr]))
            r = r * (1-math.exp(-lf*config.times[pr][l[i]]/freq[l[i]][pr]))

        rel = rel * (1-r)
    
    return rel, e

# def find(freq, task, pr, rt):

#     f = 0.01
#     lf = config.L[pr]*10**(config.D[pr]*(1-f)/(1-config.F_MIN[pr]))
#     r2 = math.exp(-lf*config.times[pr][task]/f)

#     f = 0.02
#     lf = config.L[pr]*10**(config.D[pr]*(1-f)/(1-config.F_MIN[pr]))
#     r3 = math.exp(-lf*config.times[pr][task]/f)

#     f = freq[task][pr]
#     lf = config.L[pr]*10**(config.D[pr]*(1-f)/(1-config.F_MIN[pr]))
#     r = math.exp(-lf*config.times[pr][task]/f)
#     r1 = r

#     epsilon = 10**(-5)
#     while(r-rt > epsilon and f >= config.F_MIN[pr]):
#         print(r, rt)
#         lf = config.L[pr]*10**(config.D[pr]*(1-f)/(1-config.F_MIN[pr]))
#         r = math.exp(-lf*config.times[pr][task]/f)
#         f = f -  0.0001

#     return f

def find(freq, task, pr, rt):
    freq_list = np.arange(config.F_MIN[pr], freq[task][pr] + config.freq_step, config.freq_step).tolist()
    f = freq[task][pr] 
    lf = config.L[pr]*10**(config.D[pr]*(1-f)/(1-config.F_MIN[pr]))
    r = math.exp(-lf*config.times[pr][task]/f)

    s = 0
    fin = len(freq_list) - 1
    
    while(abs(r-rt > config.epsilon) and s < fin):
        # print(r, rt, s, fin)
        mid = int((s+fin)/2)
        f = freq_list[mid]
        lf = config.L[pr]*10**(config.D[pr]*(1-f)/(1-config.F_MIN[pr]))
        r = math.exp(-lf*config.times[pr][task]/f)

        if(r > rt):
            fin = mid
        
        elif(r < rt):
            s = mid
        
        else:
            break

    
    if(r < rt):
        f += config.freq_step
        
    return f

def fa(alloc, R):
    freq = []
    for i in range(len(alloc)):
        freq.append([])
        for j in range(len(config.times)):
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
                    lf = config.L[k]*10**(config.D[k]*(1-freq[i][k])/(1-config.F_MIN[k]))
                    r1 = r1 * (1-math.exp(-lf*config.times[k][i]/freq[i][k]))

                    if(k != j):
                        r2 = r2 * (1-math.exp(-lf*config.times[k][i]/freq[i][k]))
                
                rt = 1-(1-(1-r1)*R/r)/r2
                e1 = (config.PS[j] + config.C[j]*(freq[i][j]**config.A[j]))*config.times[j][i]/freq[i][j]
                fq = find(freq, i, j, rt)
                e2 = (config.PS[j] + config.C[j]*(fq**config.A[j]))*config.times[j][i]/fq

                if(e2 - e1 < min_qty):
                    task = i
                    proc = j
                    fr = fq
                    min_qty = e2 - e1  

        fnew = (freq[task][proc]+fr)/2
        if(abs(fnew - freq[task][proc]) >=  0.0001):
            freq[task][proc] = fnew
        else:
            freq[task][proc] = fr
        r, _ = simulate_fault(alloc, freq)

        if(min_qty == 0):
            break
    
    return freq