using Plots
using LaTeXStrings

# backend(:pgfplotsx)
# backend(:pyplot)
backend(:gr)

#reset past font scaling
Plots.scalefontsizes()
#scale fonts
Plots.scalefontsizes(1.25)
#Extracting stored varibles

#Add colors
lightblue = colorant"rgb(63%,74%,78%)"
orange = colorant"rgb(85%,55%,13%)"
silver = colorant"rgb(69%,67%,66%)"
rust = colorant"rgb(72%,26%,6%)"

x1 = 0.65
y1 = 0.65

u1 = -0.5 
v1 = -0.3 

u1_log = u1 / x1  
v1_log = v1 / y1 

x2 = 0.35
y2 = 0.75

u2 = -0.8 
v2 = -0.8 

u2_log = u2 / x2  
v2_log = v2 / y2 


x3 = 0.2
y3 = 0.2

u3 = 0.6 
v3 = -0.40 

u3_log = u3 / x3  
v3_log = v3 / y3 


times = 0.0: 0.1:1.0
times = times[2:end]


trajectory_1 = hcat([[x1 + t * u1, y1 + t * v1] for t in times]...)'
derivatives_1 = hcat([[u1, v1] for t in times]...)'

start_1 = Matrix([x1, y1]')
derivatives_start_1 = Matrix([u1, v1]')

trajectory_2 = hcat([[x2 + t * u2, y2 + t * v2] for t in times]...)'
derivatives_2 = hcat([[u2, v2] for t in times]...)'


start_2 = Matrix([x2, y2]')
derivatives_start_2 = Matrix([u2, v2]')

trajectory_3 = hcat([[x3 + t * u3, y3 + t * v3] for t in times]...)'
derivatives_3 = hcat([[u3, v3] for t in times]...)'

start_3 = Matrix([x3, y3]')
derivatives_start_3 = Matrix([u3, v3]')



trajectory_1_log = hcat([exp.([log(x1) + t * u1_log, log(y1) + t * v1_log]) for t in times]...)'

derivatives_1_log = hcat([[u1_log * trajectory_1_log[k, 1], v1_log * trajectory_1_log[k, 2]] for k = 1 : length(times)]...)'

trajectory_2_log = hcat([exp.([log(x2) + t * u2_log, log(y2) + t * v2_log]) for t in times]...)'

derivatives_2_log = hcat([[u2_log * trajectory_2_log[k, 1], v2_log * trajectory_2_log[k, 2]] for k = 1 : length(times)]...)'

trajectory_3_log = hcat([exp.([log(x3) + t * u3_log, log(y3) + t * v3_log]) for t in times]...)'

derivatives_3_log = hcat([[u3_log * trajectory_3_log[k, 1], v3_log * trajectory_3_log[k, 2]] for k = 1 : length(times)]...)'

α = 0.10

plt = scatter(xlims=(0,1), ylims=(0,1), aspect_ratio=:equal, grid=false, thickness_scaling=1.25)


scatter!(plt, start_1[:, 1], start_1[:, 2], color=lightblue, label="")
quiver!(plt, start_1[:, 1], start_1[:, 2], quiver=(α * derivatives_start_1[: , 1], α * derivatives_start_1[:, 2]), color=lightblue, label="")

scatter!(plt, trajectory_1[:, 1], trajectory_1[:, 2], color=silver, label="")
quiver!(plt, trajectory_1[:, 1], trajectory_1[:, 2], quiver=(α * derivatives_1[: , 1], α * derivatives_1[:, 2]), color=silver, label="")

scatter!(plt, trajectory_1_log[:, 1], trajectory_1_log[:, 2], color=rust, label="")
quiver!(plt, trajectory_1_log[:, 1], trajectory_1_log[:, 2], quiver=(α * derivatives_1_log[: , 1], α * derivatives_1_log[:, 2]), color=rust, label="")



scatter!(plt, start_2[:, 1], start_2[:, 2], color=lightblue, label="")
quiver!(plt, start_2[:, 1], start_2[:, 2], quiver=(α * derivatives_start_2[: , 1], α * derivatives_start_2[:, 2]), color=lightblue, label="")

scatter!(plt, trajectory_2[:, 1], trajectory_2[:, 2], color=silver, label="")
quiver!(plt, trajectory_2[:, 1], trajectory_2[:, 2], quiver=(α * derivatives_2[: , 1], α * derivatives_2[:, 2]), color=silver, label="")
scatter!(plt, trajectory_2_log[:, 1], trajectory_2_log[:, 2], color=rust, label="")
quiver!(plt, trajectory_2_log[:, 1], trajectory_2_log[:, 2], quiver=(α * derivatives_2_log[: , 1], α * derivatives_2_log[:, 2]), color=rust, label="")



scatter!(plt, start_3[:, 1], start_3[:, 2], color=lightblue, label="")
quiver!(plt, start_3[:, 1], start_3[:, 2], quiver=(α * derivatives_start_3[: , 1], α * derivatives_start_3[:, 2]), color=lightblue, label="")

scatter!(plt, trajectory_3[:, 1], trajectory_3[:, 2], color=silver, label="")
quiver!(plt, trajectory_3[:, 1], trajectory_3[:, 2], quiver=(α * derivatives_3[: , 1], α * derivatives_3[:, 2]), color=silver, label="")
scatter!(plt, trajectory_3_log[:, 1], trajectory_3_log[:, 2], color=rust, label="")
quiver!(plt, trajectory_3_log[:, 1], trajectory_3_log[:, 2], quiver=(α * derivatives_3_log[: , 1], α * derivatives_3_log[:, 2]), color=rust, label="")

scatter!(plt, [], [], color=lightblue, label=L"\mathrm{initialization}")
scatter!(plt, [], [], color=silver, label=L"\mathrm{primal} \ \mathrm{geometry}")
scatter!(plt, [], [], color=rust, label=L"\mathrm{dual} \ \mathrm{geometry}")

savefig(plt, "./out/figures/dual_geometry.pdf")
run(`pdfcrop --margins '5 5 5 5' ./out/figures/dual_geometry.pdf ./out/figures/dual_geometry.pdf`)
