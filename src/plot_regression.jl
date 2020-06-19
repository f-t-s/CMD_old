using Plots
# using Random
using LaTeXStrings
using JLD 
# include("./cmd.jl")
# include("./simmd.jl")
# include("./extramd.jl")
# include("./pextramd.jl")
# include("./AutoDiffUtils.jl")
# include("./Potentials.jl")

lightblue = colorant"rgb(63%,74%,78%)"
orange = colorant"rgb(85%,55%,13%)"
silver = colorant"rgb(69%,67%,66%)"
rust = colorant"rgb(72%,26%,6%)"
brown = colorant"rgb(57%,34%,16%)"

α = 100
β = 1

#reset past font scaling
Plots.scalefontsizes()
#scale fonts
Plots.scalefontsizes(2.00)

#Extracting stored varibles
ld = load("./out/jld/regression_alpha_$(α)_beta_$(β).jld")
for key in keys(ld)
  symb = Symbol(key)
  data = ld[key]
  @eval $(symb) = $(data)
end

iterplot = plot(evals_cmd, label=L"\mathrm{CMW}", color=silver, title=latexstring("\\alpha = $(α), \\beta = $(β)"), linewidth=5, ylims=(0,40))
plot!(iterplot, evals_extragd, label=L"\mathrm{PX}", color=orange, linestyle=:dot, linewidth=5)
plot!(iterplot, evals_extramd, label=L"\mathrm{PXM}", color=rust, linestyle=:dash, linewidth=5)

savefig(iterplot, "./out/figures/regression_iterplot_alpha_$(α)_beta_$(β).pdf")

gradplot = plot(Int.(iters_cmd),  evals_cmd, label=L"\mathrm{CMW}", color=silver, title=latexstring("\\alpha = $(α), \\beta = $(β)"), linewidth=5, ylims=(0,40))
plot!(gradplot, Int.(iters_extragd), evals_extragd, label=L"\mathrm{PX}", color=orange, linestyle=:dot, linewidth=5)

plot!(gradplot, Int.(iters_extramd), evals_extramd, label=L"\mathrm{PXM}", color=rust, linestyle=:dash, linewidth=5)

savefig(gradplot, "./out/figures/regression_gradplot_alpha_$(α)_beta_$(β).pdf")


