# FEM
1.Requirements analysis: 
Develop an FEM program to calculate the nonlinear elastic response of a bar subject to dynamic loading
using an implicit solution based on the Newmark-beta type solver. Include
a) An Ogden elasticity model for the material (plot the stress-strain response for the selected
parameters).
b) Verification (convergence) based on comparison to analytic solutions.
c) Provide the solution to sudden sinusoidal loading of amplitude in the nonlinear range. Assume the
bar is fixed at one end and the loading is suddenly applied at the other end.
2.Program design:
2.1 Basic idea:
The implicit solution based on the Newmark-beta type solver is similar to Newton-Raphson method, but it is used on solving differential equation of matrix.
2.2 Data structure:
Array.
2.4 Program structure: 
On the basis of last question, I add another class to achieve these requirements. 
2.4.1 
class node
{
public:
	double location;//location 
	double a;//cross-sectional area
	double E,rho,u,strain;//moudul,density,displacement
	double Constitution(double strain);//return stress
	double Etan;
};
Save the basic parameter of a material node.
