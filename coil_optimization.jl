using LinearAlgebra

using JuMP, Juniper, Ipopt
ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
optimizer = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt)

# Voltage
V = 5.0


# Resistivity 
ρ = 1.724 * 10^-8


PCBside = 0.1 # 10cm
Max_side_coil = PCBside
 

model = Model(optimizer)

min_feature_width = 0.000089 # 1oz thickness - 0.089 mm
min_feature_width = 0.0002 # 2oz thickness - 0.2 mm

@variable(model, N, Int)
@variable(model, 0 <= I <= 1.0)
@variable(model, min_feature_width <= t_w <= 0.01)
@variable(model, min_feature_width <= g_w <= 0.01)

w = t_w + g_w
A = (PCBside - N * w)^2
l = 4 * (PCBside - N * w)
R = ρ * l / A

@objective(model, Max,  N * I * A - I^2 * R) # - I^2 * R +


@constraint(model, N * w  <= 0.5 * PCBside)
#@constraint(model, V >= I * R)


optimize!(model)


println("====================================")
println(":::Optimization results:::")
println(termination_status(model))
println("Objective_value: ", objective_value(model))
N = value.(N)
I = value.(I)
t_w = value.(t_w)
g_w = value.(g_w)
w = t_w + g_w
A = (PCBside - N * w)^2

println("N: ", N)
println("I: ", I)
println("trace width: ", t_w)
println("gap width: ", g_w)

R = ρ * 4 * (PCBside - N * w) / ((PCBside - N * w)^2)

println("Cross section area (m^2): ", A)
println("Magnetic moment (A m^2): ", N * I * A)
println("Resistance (Ohm): ", R)
println("Resistive Loss (W): ", I^2 * R)