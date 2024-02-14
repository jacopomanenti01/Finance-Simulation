function y = plrnd(u, xmin, alpha)

y = xmin.*(1.-u).^(-1/(alpha-1));
end