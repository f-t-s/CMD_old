using Plots
using Random
using LaTeXStrings
using JLD 
include("./cmd.jl")
include("./simmd.jl")
include("./extramd.jl")
include("./pextramd.jl")
include("./AutoDiffUtils.jl")
include("./Potentials.jl")

@show α = 10
@show β = 1000

#reset past font scaling
Plots.scalefontsizes()
#scale fonts
Plots.scalefontsizes(1.75)

Random.seed!(123)

lightblue = colorant"rgb(63%,74%,78%)"
orange = colorant"rgb(85%,55%,13%)"
silver = colorant"rgb(69%,67%,66%)"
rust = colorant"rgb(72%,26%,6%)"
brown = colorant"rgb(57%,34%,16%)"

max_iter = 3000

m = 5000
n_terms = 50
n = 1

A = randn(n_terms, m)
b = vec(A[:, 1] + A[:, 2]) / 2 + 0.1 * randn(n_terms)

f(x, y) = sum((A * x - b).^2) + y[1] * (sum(x) - 1)
f_eval(x, y) = sum(abs.(A * x / sum(x) - b))
g(x, y) = - y[1] * (sum(x) - 1)

x0 = ones(m) / m
y0 = zeros(1)

let 
    global scatter_cmd = plot()
    global xHist_cmd = zeros(m, max_iter + 1)
    global yHist_cmd = zeros(n, max_iter + 1)
    global iters_cmd = zeros(Int, max_iter + 1)
    global m
    global n
    global A
    global f
    global g
    global x0
    global y0
   

    ψ_x = shannon_entropy(α)
    ψ_y = squared_distance(β)

    x = ψ_x[2](x0)
    y = ψ_y[2](y0)

    opt = CMD(x, ψ_x, f, y, ψ_y, g, cg!)

    for k = 1 : max_iter
      xHist_cmd[:, k], yHist_cmd[:, k] = iterate_primal(opt)
      iters_cmd[k+1] = iters_cmd[k] + step!(opt)

      opt.x̄ .= ψ_x[2](max.(iterate_primal(opt)[1], 1e-30))
      opt.x̄ .= ψ_x[2](min.(iterate_primal(opt)[1], 1e30))
    end
    xHist_cmd[:, end], yHist_cmd[:, end]= iterate_primal(opt)
end

let
    global xHist_extragd = zeros(m, max_iter + 1)
    global yHist_extragd = zeros(n, max_iter + 1)
    global iters_extragd = zeros(Int, max_iter + 1)
    global m
    global n
    global A
    global f
    global g
    global x0
    global y0
   

    ψ_x = positive_projection(α)
    ψ_y = squared_distance(β)

    x = ψ_x[2](x0)
    y = ψ_y[2](y0)


    opt = ExtraMD(x, ψ_x, f, y, ψ_y, g)

    for k = 1 : max_iter
      xHist_extragd[:, k], yHist_extragd[:, k] = iterate_primal(opt)
      iters_extragd[k+1] = iters_extragd[k] + step!(opt)

      opt.x̄ .= ψ_x[2](max.(iterate_primal(opt)[1], 1e-30))
      opt.x̄ .= ψ_x[2](min.(iterate_primal(opt)[1], 1e30))
    end
    xHist_extragd[:, end], yHist_extragd[:, end]= iterate_primal(opt)
end

let
    global xHist_extramd = zeros(m, max_iter + 1)
    global yHist_extramd = zeros(n, max_iter + 1)
    global iters_extramd = zeros(Int, max_iter + 1)
    global m
    global n
    global A
    global f
    global g
    global x0
    global y0
   


    ψ_x = shannon_entropy(α)
    ψ_y = squared_distance(β)

    x = ψ_x[2](x0)
    y = ψ_y[2](y0)


    opt = ExtraMD(x, ψ_x, f, y, ψ_y, g)

    for k = 1 : max_iter
      xHist_extramd[:, k], yHist_extramd[:, k] = iterate_primal(opt)
      iters_extramd[k+1] = iters_extramd[k] + step!(opt)

      opt.x̄ .= ψ_x[2](max.(iterate_primal(opt)[1], 1e-30))
      opt.x̄ .= ψ_x[2](min.(iterate_primal(opt)[1], 1e30))
    end
    xHist_extramd[:, end], yHist_extramd[:, end]= iterate_primal(opt)
end

evals_cmd = [f_eval(xHist_cmd[:, k], yHist_cmd[:, k]) for k in 1 : size(xHist_cmd, 2)]
evals_extramd = [f_eval(xHist_extramd[:, k], yHist_extramd[:, k]) for k in 1 : size(xHist_extramd, 2)]
evals_extragd = [f_eval(xHist_extragd[:, k], yHist_extragd[:, k]) for k in 1 : size(xHist_extragd, 2)]

plot(evals_cmd, label="cmd", color=lightblue)
plot!(evals_extramd, label="extramd", color=orange)
plot!(evals_extragd, label="extragd", color=rust)

save("./out/jld/regression_alpha_$(α)_beta_$(β).jld",
#      "f", f,
#      "f_eval", f_eval,
#      "g", g,
     "A", A,
     "b", b,
     "xHist_cmd", xHist_cmd,
     "yHist_cmd", yHist_cmd,
     "iters_cmd", iters_cmd,
     "evals_cmd", evals_cmd,
     "xHist_extramd", xHist_extramd,
     "yHist_extramd", yHist_extramd,
     "iters_extramd", iters_extramd,
     "evals_extramd", evals_extramd,
     "xHist_extragd", xHist_extragd,
     "yHist_extragd", yHist_extragd,
     "iters_extragd", iters_extragd,
     "evals_extragd", evals_extragd )