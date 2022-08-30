# Renyi-Dimension-of-Chaotic-Systems
This code Calculates Renyi Dimension for chaotic systems. As chaotic system are often characterized by their multi fractal structure, one scaling law can not describe their complexity. Thus, their dynamic should be described by a more general definition. One the most prevalent methods for measering the complexity of chatic systems is Renyi Dimension. 
Renyi Dimension calculates the range of scaling laws for multifractal systems and is defined as:
$$D_q = {1 \over {q-1}} lim{ ln{ \sum_{i}P_i^q} \over ln{1 \over l} } $$
where P_i is the probability of a point in state space, l is the radius of the sphere and q is the order of Renyi entropy. Parameter q captures different fractal dimensions D_q in the multifractal structure. As a result, the general dimension is not just one value but a function of the parameter q.
# How to Run
D_alpha function calculates Renyi Dimension. The input arguments of this function are:
- signal:                     Input signal (There must be enough data inputs otherwise this function would not work properly)                       
- alpha:                      the order of the entropy     
- L0:                         Length of the initial box
- scaling_factor:             the scaling factor which shrinks the size of boxes at each iteration
- iteration:                  number of times scaling factor applies
- DimensionofSignal:          Defines the dimension of system. should be "1D" , "2D" or "3D"
- plotEntropyRadiusOption:    if true plots the logEntropy_logRadius graph
