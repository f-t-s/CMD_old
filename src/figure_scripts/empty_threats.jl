using Plots
using LaTeXStrings

# backend(:pgfplotsx)
backend(:pyplot)
# backend(:gr)

#reset past font scaling
Plots.scalefontsizes()
#scale fonts
Plots.scalefontsizes(5.00)
#Extracting stored varibles

#Add colors
lightblue = colorant"rgb(63%,74%,78%)"
orange = colorant"rgb(85%,55%,13%)"
silver = colorant"rgb(69%,67%,66%)"
rust = colorant"rgb(72%,26%,6%)"

f(x,y) = 2 * x * y - 4 * sum((1 - y).^2) + 100 * Float64(x < 0)
flocal(x,y) = f(0.0, 1.0) + 2 * (x - 0.0) + 2 * (x - 0.0) * (y - 1.0)

x1 = -0.5 : 0.05: 1.5
y1 = - 0.5 : 0.05 : 1.5 


x2 = -0.5 : 0.10: 1.5
y2 = -0.5 : 0.10 : 1.5 

srf = wireframe(x1, y1, f, zlims=(-3,5), camera=(15,35), linecolor=silver, linealpha=0.2, xlabel=L"x", ylabel=L"y", label=L"\mathrm{true} \ f)", thickness_scaling=1.5, linewidth=0.75)

wireframe!(srf, x2, y2, flocal, linecolor=rust, label=L"\mathrm{approximate} \ f)", legend=:topright)

savefig(srf, "./out/figures/empty_threats.pdf")
run(`pdfcrop --margins '5 5 5 5' ./out/figures/empty_threats.pdf ./out/figures/empty_threats.pdf`)
