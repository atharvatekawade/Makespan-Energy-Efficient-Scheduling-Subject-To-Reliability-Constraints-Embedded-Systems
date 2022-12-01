import config, utils, max_re, rr, ods, mecrg, esrg, mert, eafts
import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

#-db DATABASE -u USERNAME -p PASSWORD -size 20000
parser.add_argument("-rho", dest = "rho", default = 3, type=int)
parser.add_argument("-pr", dest = "pr", default = 16, type=int)
parser.add_argument("-R", dest = "R", default = 0.9, type=float)
parser.add_argument("-smin", dest = "smin", default = 10, type=int)
parser.add_argument("-smax", dest = "smax", default = 10**2, type=int)
parser.add_argument("-ods_start", dest = "ods_start", default = 10, type=float)
parser.add_argument("-ods_end", dest = "ods_end", default = 100, type=float)
parser.add_argument("-ods_step", dest = "ods_step", default = 10, type=float)
parser.add_argument("-mert_start", dest = "mert_start", default = 0, type=float)
parser.add_argument("-mert_end", dest = "mert_end", default = 1, type=float)
parser.add_argument("-mert_step", dest = "mert_step", default = 0.1, type=float)


args = parser.parse_args()

print("Reliability constraint:",args.R)

config.t1 = args.smin
config.t2 = args.smax

config.init(args.pr, args.rho)

alloc, _ = utils.mr()
_, max_rel1, _ = utils.simulate(alloc)
_, max_rel2, max_energy = utils.mr_fault()

spans = [0, 0, 0, 0, 0]
energies = [0, 0, 0, 0, 0]
rels = [0, 0, 0, 0, 0]

energies_ft = [0, 0, 0, 0]
rels_ft = [0, 0, 0, 0]
R = args.R

if(R <= max_rel1):
    print("Running ODS")
    ods.ods(args.ods_start, args.ods_end, args.ods_step)

    alloc, _ = utils.mr()
    freq = utils.soea1(alloc, R)
    f1, r1, e1 = utils.simulate(alloc, freq = freq)

    spans[0] += f1
    energies[0] += e1
    rels[0] += r1

    print(f"Algo: MR No. of vertices: {len(config.graph)} Span: {f1} Rel: {r1} Energy: {e1}")

    alloc, _ = utils.lec()
    f1, r1, e5 = utils.simulate(alloc)
    if(R - r1 > config.epsilon):
        alloc, _ = utils.mr()
        freq = utils.soea1(alloc, R)
        f1, r1, e5 = utils.simulate(alloc, freq = freq)

    print(f"Algo: LEC No. of vertices: {len(config.graph)} Span: {f1} Rel: {r1} Energy: {e5}")

    alloc, freq = esrg.esrg(R)
    f1, r1, e2 = utils.simulate(alloc, freq = freq)
    if(R - r1 > config.epsilon):
        alloc, _ = utils.mr()
        freq = utils.soea1(alloc, R)
        f1, r1, e2 = utils.simulate(alloc, freq = freq)

    spans[1] += f1
    energies[1] += e2
    rels[1] += r1

    print(f"Algo: ESRG No. of vertices: {len(config.graph)} Span: {f1} Rel: {r1} Energy: {e2}")

    alloc = mecrg.mecrg(R)
    freq = utils.soea1(alloc, R)
    f1, r1, e3 = utils.simulate(alloc, freq = freq)
    spans[2] += f1
    energies[2] += e3
    rels[2] += r1
    print(f"Algo: MECRG No. of vertices: {len(config.graph)} Span: {f1} Rel: {r1} Energy: {e3}")

    E = min(e1, e2, e3)

    f2 = float('inf')
    r2 = 0
    e2 = float('inf')

    f3 = float('inf')
    r3 = 0
    e3 = float('inf')

    for j in range(len(config.ods_results)):
        if(config.ods_results[j][1] >= R):
            freq = utils.soea1(config.ods_alloc[j], R)
            f, r, e = utils.simulate(config.ods_alloc[j], freq = freq)
            if(e <= E):
                if(f < f2):
                    f2 = f
                    r2 = r
                    e2 = e
                        
                elif(f == f2 and e < e2):
                    f2 = f
                    r2 = r
                    e2 = e
                    
            else:
                if(e < e3):
                    f3 = f
                    r3 = r
                    e3 = e
                        
                elif(e == e3 and f < f3):
                    f3 = f
                    r3 = r
                    e3 = e

    if(f2 == float('inf')):
        if(f3 ==float('inf')):
            alloc, _ = utils.mr()
            freq = utils.soea1(alloc, R)
            f2, r2, e2 = utils.simulate(alloc, freq = freq)
                
        else:
            f2 = f3
            e2 = e3
            r2 = r3

    spans[3] += f2
    energies[3] += e2
    rels[3] += r1

    print(f"Algo: ODS No. of vertices: {len(config.graph)} Span: {f2} Rel: {r2} Energy: {e2}")

    th_list = np.arange(args.mert_start, args.mert_end + args.mert_step, args.mert_step).tolist()
    for a in range(-1, len(config.graph)):
        for th in th_list:
            mert.mert(a = a, th = th, R = R)

    f2 = float('inf')
    r2 = 0
    e2 = float('inf')

    f3 = float('inf')
    r3 = 0
    e3 = float('inf')

    for j in range(len(config.my_results)):
        if(config.my_results[j][1] >= R):
            freq = utils.soea1(config.my_alloc[j], R)
            f, r, e = utils.simulate(config.my_alloc[j], freq = freq)
                    
            if(e <= E):
                if(f < f2):
                    f2 = f
                    r2 = r
                    e2 = e
                        
                elif(f == f2 and e < e2):
                    f2 = f
                    r2 = r
                    e2 = e
                    
            else:
                if(e < e3):
                    f3 = f
                    r3 = r
                    e3 = e
                        
                elif(e == e3 and f < f3):
                    f3 = f
                    r3 = r
                    e3 = e
        else:
            print("Something went wrong in results of mert")
            quit()

    if(f2 == float('inf')):
        f2 = f3
        e2 = e3
        r2 = r3

    spans[4] += f2
    energies[4] += e2
    rels[4] += r1

    print(f"Algo: MERT No. of vertices: {len(config.graph)} Span: {f2} Rel: {r2} Energy: {e2}")

    plt.subplot(1, 3, 1)
    algos = ["MR" , "LEC", "ESRG", "MECRG", "ODS", "MERT"]
    plt.bar(algos, spans, color ='blue')
    # plt.xlabel("Algorithms")
    plt.ylabel("Makespan")
    plt.title("Makespan comparison")

    plt.subplot(1, 3, 2)
    plt.bar(algos, energies, color ='red')
    # plt.xlabel("Algorithms")
    plt.ylabel("Energy")
    plt.title("Energy comparison")

    plt.subplot(1, 3, 3)
    plt.bar(algos + ["Constraint"], rels + [R], color ='maroon')
    # plt.xlabel("Algorithms")
    plt.ylabel("Rel")
    plt.title("Reliability comparison")

    plt.show()


elif(R <= max_rel2):
    alloc, _, _ = max_re.max_re(R)
    freq = utils.fa(alloc, R)
    rel, energy = utils.simulate_fault(alloc, freq)
    energies_ft[0] += energy
    rels_ft[0] += rel
    print(f"Graph: {len(config.graph)} Algo: Max-Re Rel:{rel} Energy: {energy}")

    alloc, _, _ = rr.rr(R)
    freq = utils.fa(alloc, R)
    rel, energy = utils.simulate_fault(alloc, freq)
    energies_ft[1] += energy
    rels_ft[1] += rel
    print(f"Graph: {len(config.graph)} Algo: RR Rel:{rel} Energy: {energy}")


    alloc, freq, rel, energy = esrg.efsrg(R)
    energies_ft[2] += energy
    rels_ft[2] += rel
    print(f"Graph: {len(config.graph)} Algo: ESFRG Rel:{rel} Energy: {energy}")


    alloc, _, _ = eafts.eafts(R)
    freq = utils.fa(alloc, R)
    rel, energy = utils.simulate_fault(alloc, freq)
    energies_ft[3] += energy
    rels_ft[3] += rel
    print(f"Graph: {len(config.graph)} Algo: EAFTS Rel:{rel} Energy: {energy}")

    plt.subplot(1, 3, 1)
    algos = ["Max-Re" , "RR", "EFSRG", "EAFTS"]
    plt.bar(algos, energies_ft, color ='blue')
    # plt.xlabel("Algorithms")
    plt.ylabel("Energy")
    plt.title("Energy comparison")

    plt.subplot(1, 3, 2)
    plt.bar(algos + ["Constraint"], rels_ft + [R], color ='maroon')
    # plt.xlabel("Algorithms")
    plt.ylabel("Rel")
    plt.title("Reliability comparison")

    plt.show()


else:
    print("Maximum reliability achievable:", max_rel2)


