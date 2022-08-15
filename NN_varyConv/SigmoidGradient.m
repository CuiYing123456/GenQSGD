function g = SigmoidGradient(z)

g = zeros(size(z));
g=(1./(1+exp(-z))).*(exp(-z)./(1+exp(-z)));

end
