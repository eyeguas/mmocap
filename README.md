# mmocap
Parallelization Strategies for Markerless Human Motion Capture

Markerless Motion Capture (MMOCAP) is the
problem of determining the pose of a person from images
captured by one or several cameras simultaneously without
using markers on the subject. Evaluation of the solutions
is frequently the most time-consuming task, making most
of the proposed methods inapplicable in real-time scenarios. This paper presents an efficient approach to parallelize
the evaluation of the solutions in CPUs and GPUs. Our proposal is experimentally compared on six sequences of the
HumanEva-I dataset using the CMAES algorithm. Multiple algorithm’s configurations were tested to analyze the
best trade-off in regard to the accuracy and computing time.
The proposed methods obtain speedups of 8× in multi-core
CPUs, 30× in a single GPU and up to 110× using 4 GPUs.

This paper presented an efficient and parallelizable approach
to evaluate the solutions in the MMOCAP problem. Our approach consists in approximating the triangle body meshes
by rectangular patches that are easily drawn and computed.
In addition, strategies to parallelize the computation both in
CPUs and GPUs were proposed. First, we proposed a strategy based on the CPU’s Streaming SIMD Extensions (SSE)
instruction set, which demonstrated to double performance.
Second, we proposed an strategy using multi-threading on
multi-core CPUs, which showed to speed up model evaluation up to 4 times. Third, we presented a GPU strategy
scalable to multiple devices. A total of 4 GPUs were used
to collaborate, providing a speedup of up to 110× faster
than the naive CPU code. The parallelization approaches
proposed also demonstrated better performance than an efficient OpenGL implementation.
Moreover, we experimented multiple algorithm’s configuration and body mesh resolution sizes, which allowed
for analyzing the performance impact of varying the number of evaluations and the number of body model vertices.
Accuracy and runtime were two conflicting objectives for
the MMOCAP problem, and we established a Pareto front
of solutions offering different levels of trade-off. Eventually, the user is allowed for selecting the configuration setup
which best matches his needs in terms of obtaining accurate but slow solutions, or fast and less accurate solutions. A
trade-off solution producing both accurate and fast results is
proposed as recommended configuration.
As for the future work, it would be interesting to compare the performance of the proposal with FPGA-based implementations that allow for efficient single-bit operations.
