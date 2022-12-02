# Makespan-Energy-Efficient-Scheduling-Subject-To-Reliability-Constraints-Embedded-Systems
This problem is of scheduling applications modeled as a directed acyclic graph to minimize total execution time and energy consumption subject to reliability constraints in a multiprocessor embedded system. Depending upon the reliability constraint, we follow either of the two approaches:
a) Non-Fault Tolearnt: In this approach, we assign a each task to a single processor. The objectives we optimize for in this case are makespan and energy. Our proposed method named MERT is compared with the following state-of-art algorithms:
1) (ODS, MR, LEC) Dynamic DAG Scheduling on Multiprocessor Systems: Reliability, Energy, and Makespan - Jing Huang College of Computer Science and Electronic Engineering, Hunan University, Changsha, China,Key Laboratory for Embedded and Network Computing of Hunan Province, Hunan University, Changsha, China ; Renfa Li; Xun Jiao; Yu Jiang; Wanli Chang: https://ieeexplore.ieee.org/document/9211460
2) (ESRG) Energy-Efficient Fault-Tolerant Scheduling of Reliable Parallel Applications on Heterogeneous Distributed Embedded Systems - Guoqi Xie College of Computer Science and Electronic Engineering, Hunan University, Hunan, China ; Yuekun Chen; Xiongren Xiao; Cheng Xu; Renfa Li; Keqin Li: https://ieeexplore.ieee.org/document/7938375
3) (MECRG) Resource Consumption Cost Minimization of Reliable Parallel Applications on Heterogeneous Embedded Systems - Guoqi Xie
Key Laboratory for Embedded and Network Computing of Hunan Province, Hunan University, Changsha, China
; Yuekun Chen; Yan Liu; Yehua Wei; Renfa Li; Keqin Li: https://ieeexplore.ieee.org/document/7792673

![Figure_1](https://user-images.githubusercontent.com/64606981/205252506-810ecaa7-6a05-469f-b157-bdfd9eb02ed2.png)

b) Fault-Tolearnt: In this approach, we assign a each task to a multiple processors if the constraint is not satisfied by Non-Fault Tolerant approach. We only try to optimize energy in this case. Our proposed method named EAFTS is compared with the following state-of-art algorithms:

1) (Max-Re) Reliable workflow scheduling with less resource redundancy: LaipingZhao; YizhiRen; KouichiSakurai: https://www.sciencedirect.com/science/article/abs/pii/S0167819113000732
2) (RR) Fault-tolerant scheduling with dynamic number of replicas in heterogeneous systems - Laiping Zhao
Department of Informatics, Kyushu University, Fukuoka, Japan ; Yizhi Ren; Yang Xiang; Kouichi Sakurai: https://ieeexplore.ieee.org/document/5738914
2) (EFSRG) Energy-Efficient Fault-Tolerant Scheduling of Reliable Parallel Applications on Heterogeneous Distributed Embedded Systems - Guoqi Xie College of Computer Science and Electronic Engineering, Hunan University, Hunan, China ; Yuekun Chen; Xiongren Xiao; Cheng Xu; Renfa Li; Keqin Li: https://ieeexplore.ieee.org/document/7938375


![Figure_2](https://user-images.githubusercontent.com/64606981/205253318-b36c93d4-faf4-43d2-8313-aadb34756f6a.png)

## Usage
Clone the repositary and run the command: python main.py -rho -pr -R -smin -smax -ods_start -ods_end -ods_step -wanms_start -wanms_end -wanms_step, the arguments are explained below:

rho: Represents the number of parameter for number of nodes of FFT Task graph: n = (2+rho)*2**rho - 1.
pr: Represents the number of processors in our embedded system.
R: Represents the reliability constraint < 1.
smin: Represents the lower bound for task computation requirement and edge data.
smax: Represents the upper bound for task computation requirement and edge data.
ods_start: Represents the starting value of Θ in ODS.
ods_end: Represents the ending value of Θ in ODS.
ods_step: Represents the value of step from ods_start to ods_end.
mert_start: Represents the starting value of α in MERT.
mert_end: Represents the ending value of α in MERT.
mert_step: Represents the value of step from mert_start to mert_end.
Upon running the command and successful execution, we get plots for the cost and makespan of different algorithms. Depending on the constraint, we use either of the two appraches mentioned above. The reliability plot also includes the reliability constraint for reference. Leave all flags blank for running parameters with default values.
