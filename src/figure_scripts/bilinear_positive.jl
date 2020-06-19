using Plots
using Random
using LaTeXStrings
include("../cmd.jl")
include("../simmd.jl")
include("../extramd.jl")
include("../pextramd.jl")
include("../AutoDiffUtils.jl")
include("../Potentials.jl")

lightblue = colorant"rgb(63%,74%,78%)"
orange = colorant"rgb(85%,55%,13%)"
silver = colorant"rgb(69%,67%,66%)"
rust = colorant"rgb(72%,26%,6%)"

m = 1
n = 1

true_x = [0.1]
# true_x[randperm(m)[1:20]] .= 0.0
true_y = [0.1]
# true_y[randperm(m)[1:20]] .= 0.0

A = ones(m, n) 

l3 = 3.0

f(x, y) = 0.1 * 4^l3 * sum((x .- true_x)' * A * (y .- true_y))
g(x, y) = -f(x, y)

x0 = [0.4]
y0 = [0.3]


let
    max_iter = 5000
    global scatter_cmd = plot()
    global xHist_cmd = zeros(m, max_iter + 1)
    global yHist_cmd = zeros(n, max_iter + 1)
    global m
    global n
    global A
    global f
    global g
    global x0
    global true_x
    global y0
    global true_y
   

    α = 1.0

    ψ_x = shannon_entropy(α)
    ψ_y = shannon_entropy(α)


    x = ψ_x[2](x0)
    y = ψ_y[2](y0)


    opt = CMD(x, ψ_x, f, y, ψ_y, g, cg!)

    grad_computations = zeros(Int, max_iter + 1)
    for k = 1 : max_iter
      xHist_cmd[:, k], yHist_cmd[:, k] = iterate_primal(opt)
      # grad_computations[k+1] = grad_computations[k] + step!(opt)
      step!(opt)

      opt.x̄ .= ψ_x[2](max.(iterate_primal(opt)[1], 0.0))
      opt.ȳ .= ψ_y[2](max.(iterate_primal(opt)[2], 0.0))
      # @show grad_computations[k+1] - grad_computations[k]
    end
    xHist_cmd[:, end], yHist_cmd[:, end]= iterate_primal(opt)

    @show grad_computations[end]/max_iter


    scatter!(scatter_cmd, xHist_cmd[1, :], yHist_cmd[1, :], xlims=(0.0, 0.6), ylims=(0.0, 0.6))
end

let
    max_iter = 5000
    global scatter_extramd = plot()
    global xHist_extramd = zeros(m, max_iter + 1)
    global yHist_extramd = zeros(n, max_iter + 1)
    global m
    global n
    global A
    global f
    global g
    global x0
    global true_x
    global y0
    global true_y
   

    α = 1.0

    ψ_x = positive_projection(α)
    ψ_y = positive_projection(α)

    x = ψ_x[2](x0)
    y = ψ_y[2](y0)


    opt = ExtraMD(x, ψ_x, f, y, ψ_y, g)

    iters = zeros(Int, max_iter + 1)
    for k = 1 : max_iter
      xHist_extramd[:, k], yHist_extramd[:, k] = iterate_primal(opt)
      # iters[k+1] = iters[k] + step!(opt)
      step!(opt)

      opt.x̄ .= ψ_x[2](max.(iterate_primal(opt)[1], 0.0))
      opt.ȳ .= ψ_y[2](max.(iterate_primal(opt)[2], 0.0))
    end
    @show xHist_extramd[:, end], yHist_extramd[:, end]= iterate_primal(opt)


     scatter!(scatter_extramd, xHist_extramd[1, :], yHist_extramd[1, :], xlims=(0.0, 0.6), ylims=(0.0, 0.6))
end

out_plot = scatter(xHist_cmd[1, :], yHist_cmd[1, :], xlims=(0.0, 0.6), ylims=(0.0, 0.6), markercolor=silver, aspect_ratio=:equal, size=(500,500), label=L"\mathrm{CMW}", thickness_scaling=2.00, markerstrokewidth=0.1)
scatter!(out_plot, xHist_extramd[1, :], yHist_extramd[1, :], xlims=(0.0, 0.6), ylims=(0.0, 0.6), markercolor=rust, aspect_ratio=:equal, label=L"\mathrm{Extra gradient}", markerstrokewidth=0.1)

savefig(out_plot, "./out/figures/bilinear_positive_3to$(Int(l3)).pdf")
run(`pdfcrop --margins '5 5 5 5' ./out/figures/bilinear_positive_3to$(Int(l3)).pdf ./out/figures/bilinear_positive_3to$(Int(l3)).pdf`)

