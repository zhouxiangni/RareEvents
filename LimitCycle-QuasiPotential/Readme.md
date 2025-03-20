

The source code for the published work 

**Quasi-Potential Calculation and Minimum Action Method for Limit Cycle**
Ling Lin · Haijun Yu · Xiang Zhou
*Journal of Nonlinear Science (2019) 29:961–991*
https://doi.org/10.1007/s00332-018-9509-3
Abstract
We study the noise-induced escape from a stable limit cycle of a non-gradient dynam- ical system driven by a small additive noise. The fact that the optimal transition path in this case is infinitely long imposes a severe numerical challenge to resolve it in the minimum action method. We first consider the landscape of the quasi-potential near the limit cycle, which characterizes the minimal cost of the noise to drive the system far away form the limit cycle. We derive and compute the quadratic approximation of this quasi-potential near the limit cycle in the form of a positive definite solution to a matrix-valued periodic Riccati differential equation on the limit cycle. We then com- bine this local approximation in the neighborhood of the limit cycle with the minimum action method applied outside of the neighborhood. The neighborhood size is selected to be compatible with the path discretization error. By several numerical examples, we show that this strategy effectively improves the minimum action method to compute the spiral optimal escape path from limit cycles in various systems.

------------
The incomplete list of implemented codes
+ Algorithm to find the limit cycle with an initial guess of the point and the period, based on the Newton Raphson method
+ Geometric minimum action method
+ Solving the Riccati equation related to Floque theorem as a quadratic expansion near the limit cycle


