using LinearAlgebra

using JuMP, Juniper, Ipopt
ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
optimizer = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt)

# Resistive loss limit 
resistive_loss_limit = 1.0

# Voltage
V = 5.0

# Copper Resistivity 
ρ = 1.724 * 10^-8

# PCB edge length (square)
pcb_side = 0.1 # 10cm

# Outermost coil edge length
coil_side_max = pcb_side

# PCB trace manufacturing constraint
min_feature_width = 0.00009 # 1oz thickness - 0.09 mm
#min_feature_width = 0.0002 # 2oz thickness - 0.2 mm

# PCB trace thickness
trace_thickness =  3.556 * 10^-5 # 1oz copper - 35um = 1.4 mils
#trace_thickness = 7.112 * 10^-5 # 2oz copper - 70um = 2.8 mils
 
# PCB layers
pcb_layers = 2


model = Model(optimizer)


@variable(model, N, Int)
@variable(model, 0 <= I <= 1.0)
@variable(model, min_feature_width <= trace_width <= coil_side_max / 2)
@variable(model, min_feature_width <= gap_width <= coil_side_max / 2)

coil_width = trace_width + gap_width
A = (coil_side_max - N * coil_width)^2 # Average coil cross-section area
coil_length = 4 * (coil_side_max - N * coil_width) * N * pcb_layers
R = ρ * coil_length / (trace_width * trace_thickness) # coil resistance

@objective(model, Max,  N * I * A * pcb_layers) # - I^2 * R +

@constraint(model, I^2 * R <= resistive_loss_limit) # .. Watts - Requirement?
@constraint(model, N * coil_width <= 0.5 * coil_side_max)
@constraint(model, V == I * R)


optimize!(model)


println("====================================")
println(":::Optimization results:::")
println(termination_status(model))
println("Objective_value: ", objective_value(model))
N = value.(N)
I = value.(I)
trace_width = value.(trace_width)
gap_width = value.(gap_width)

println("Coils per layer: ", N)
println("Total coils: ", N * pcb_layers)
println("I: ", I)
println("Trace width (m): ", trace_width)
println("Gap width (m): ", gap_width)


coil_width = trace_width + gap_width
A = (coil_side_max - N * coil_width)^2
coil_length = 4 * (coil_side_max - N * coil_width) * N * pcb_layers
R = ρ * coil_length / (trace_width * trace_thickness)

println("Total coils width (m): ", N * coil_width)
println("Total coil length (m): ", 4 * (coil_side_max - N * coil_width) * N * pcb_layers)
println("Cross section area (m^2): ", A)
println("Magnetic moment (A m^2): ", N * I * A * pcb_layers)
println("Resistance (Ohm): ", R)
println("Resistive Loss (W): ", I^2 * R)