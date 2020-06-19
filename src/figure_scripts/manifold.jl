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

N = 32
u = range(0, stop=2π, length=N)
v = range(0, stop=π, length=N)
x = cos.(u) * sin.(v)'
y = sin.(u) * sin.(v)'
z = repeat(cos.(v)',outer=[N, 1])

x_t = -0.75:0.15:0.75
y_t = -0.75:0.15:0.75
z_t = one.(x_t)


# 
ν = 0.92 * 2 * π 
v_e = range(0, stop=1.25 * π/2, length=100)
x_e = cos(ν) * sin.(v_e)
y_e = sin(ν) * sin.(v_e)
z_e = cos.(v_e)




srf = wireframe(x, y, z .-1.0, camera=(25,40), linecolor=silver, linealpha=0.2,linewidth=0.50, xticks=[], yticks=[], zticks=[], xlims=[-1.5,1.5], ylims=[-1.5,1.5], zlims=[-2.0,0.0], legend=false)

wireframe!(srf, x_t, y_t, (x,y) -> 0.0, linecolor=lightblue, label=L"\mathrm{approximate} \ f)", legend=:topright)

scatter!(srf, [0.0], [0.0], [0.0], color=orange, label="")
quiver!(srf, [0.0], [0.0], [0.0], color=orange, label="", quiver=(1.50 * [cos(ν)], 1.50 * [sin(ν)]), arrow = (0.5, 0.2), linewidth=2.0)
plot!(srf, x_e, y_e, z_e .- 1.00, color=rust, label="", linewidth=2.0)


savefig(srf, "./out/figures/manifold.pdf")
run(`pdfcrop --margins '-65 -35 -65 -55' ./out/figures/manifold.pdf ./out/figures/manifold.pdf`)
