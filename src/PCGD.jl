include("./cmd.jl")
include("./Potentials.jl")

using Random
using Plots 

# Bilinear problem:
f(x, y) = 2 * (x .- 0.00)' * (y .- 0.00) - sum((y .- 1.0).^2)
g(x, y) = -f(x, y)

m = 1
n = 1

α = 4.0

ψ_x = squared_distance(α)
ψ_y = squared_distance(α)
x = ones(m) * α
y = ones(n) * α

# ψ_x = shannon_entropy(α)
# ψ_y = shannon_entropy(α)
# x = zeros(m) * α
# y = zeros(n) * α


cmd = CMD(x, ψ_x, f, y, ψ_y, g)

max_iter = 100
xHist = zeros(m, max_iter + 1)
yHist = zeros(n, max_iter + 1)
for k = 1 : 100
  xHist[:, k], yHist[:, k] = iterate_primal(cmd)
  step!(cmd)

  cmd.x̄ .= max.(cmd.x̄, 0.0)
  cmd.ȳ .= max.(cmd.ȳ, 0.0)
end
xHist[:, end], yHist[:, end]= iterate_primal(cmd)

