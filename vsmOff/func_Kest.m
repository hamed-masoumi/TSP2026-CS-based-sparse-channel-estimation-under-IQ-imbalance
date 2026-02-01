function Kest = func_Kest(z1,z2)

alpha = norm(z1, 2)^2;
beta = norm(z2, 2)^2;
gamma = sum(z1 .* z2);
a = alpha - beta + conj(gamma) - gamma;
b = 2*gamma + beta - alpha;
c = - gamma;
Delta = b^2 - 4*a*c;
Kest = (-b + sqrt(Delta)) / (2*a);

end