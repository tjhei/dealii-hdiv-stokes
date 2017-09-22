from sympy import *
from sympy.printing import print_ccode
from sympy.physics.vector import ReferenceFrame, gradient, divergence
from sympy.vector import CoordSysCartesian

R = ReferenceFrame('R');
x = R[0]; y = R[1];

visc = Symbol("nu")
# visc=1e-1;
# print(" visc=%f" % visc)

u=[0,0]
u[0]=pi*sin(pi*x)*sin(pi*x)*sin(2*pi*y);
u[1]=-pi*sin(2*pi*x)*sin(pi*y)*sin(pi*y);
p=cos(pi*x)*sin(pi*y);
# p=p - integrate(p, (x,a,b));

grad_p = gradient(p, R).to_matrix(R)
f0 = simplify(-divergence(visc*gradient(u[0], R), R)) + grad_p[0];
f1 = simplify(-divergence(visc*gradient(u[1], R), R)) + grad_p[1];
f2 = simplify(divergence(u[0]*R.x + u[1]*R.y, R));

print("f:")
print(f0)
print(f1)
print(f2)

print("\n * RHS:")
print(ccode(f0, assign_to = "values[0]"));
print(ccode(f1, assign_to = "values[1]"));
print(ccode(f2, assign_to = "values[2]"));


print("\n * ExactSolution:")
print(ccode(u[0], assign_to = "values[0]"));
print(ccode(u[1], assign_to = "values[1]"));
print(ccode(p, assign_to = "values[2]"));

print("")
print("pressure mean:", N(integrate(integrate(p,(x,0,1)),(y,0,1))))
