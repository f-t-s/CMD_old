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

x0 = 1.0
y0 = 1.0

KL(x,y) = x * log(x / x0) + y * log(y / y0) - x - y + x0 + y0
squared_distance(x,y) = 0.5 * ((x - x0)^2 + (y - y0)^2) 

ep = 0.0001

x1 = (0.0 + ep) : 0.05 : (2.0 + ep) 
y1 = (0.0 + ep) : 0.05 : (2.0 + ep) 

x2 = (0.0 + ep) : 0.10 : (2.0 + ep)  
y2 = (0.0 + ep) : 0.10 : (2.0 + ep)

srf = wireframe(x1, y1, squared_distance, camera=(60,30), linecolor=silver, linealpha=0.2, xlabel=L"x", ylabel=L"y", zticks=0.0:0.5:2.0, label=L"\mathrm{true} \ f)", thickness_scaling=1.5, linewidth=0.75)

wireframe!(srf, x2, y2,  KL, linecolor=rust, label=L"\mathbb{D}\left(x\middle|y\right)", legend=:topright)

savefig(srf, "./out/figures/divergences.pdf")
run(`pdfcrop --margins '5 5 5 5' ./out/figures/divergences.pdf ./out/figures/divergences.pdf`)
