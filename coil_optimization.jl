using LinearAlgebra

using JuMP, Juniper, Ipopt
ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
optimizer = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt)

# Power limit 
power_limit = 1.0

# Voltage
V = 5.0

# Copper Resistivity 
ρ = 1.724 * 10^-8

# Outermost Coil Size
coil_y = 0.09570 # mm - Values from XY board
coil_x = 0.07975

# Maximum width of Coils on PCB
max_coil_width = 0.0135 # mm - Value from XY board

# PCB trace manufacturing constraint
min_trace_width = 0.0001524 # 1oz thickness - 6 mils =  0.1524 mm
gap_width       = 0.0001524 # PCB recommended capabilities
#min_feature_width = 0.0002 # 2oz thickness - 0.2 mm

# PCB trace thickness
k = 0.75 # effective copper thickness (between 0.75 and 0.8 for PCBWay)
trace_thickness =  35 * 10^-6 * k # 1oz copper - 35um = 1.4 mils
#trace_thickness = 70 * 10^-6 # 2oz copper - 70um = 2.8 mils
 
# PCB layers
pcb_layers = 2


model = Model(optimizer)


@variable(model, N, Int)
@variable(model, 0 <= I <= 1.0)
@variable(model, min_trace_width <= trace_width <= max_coil_width)

coil_width = trace_width + gap_width

A = ((coil_x * coil_y) - (coil_x * coil_width * (N - 1) / 2) - (coil_y * coil_width * (N - 1) / 2) + (coil_width * coil_width * (2 * N * N - 3 * N + 1) / 6)) * N * pcb_layers

# sum_i=1^N (x - w * (i-1)) * (y - w * (i-1))
# sum_i=1^N (xy - xw * (i-1) - yw * (i-1) + w^2 * (i-1)^2)
# xyN - xw * (N * (N+1) / 2 - N) - yw * (N * (N+1) / 2 - N) + sum_i=1^N (w^2 * (i^2 - 2i + 1))
# xyN - xwN (N - 1) / 2 - ywN (N - 1) / 2 + w^2(N - N(N+1) + N(N+1)(2N+1) / 6)
# N * (xy - xw(N-1)/2 - yw(N-1)/2 + w^2(1 - N - 1 + (2N^2 + N + 2N + 1) / 6))
# N * (xy - xw(N-1)/2 - yw(N-1)/2 + w^2(2N^2 - 3N + 1) / 6)

coil_length = 2 * (coil_x + coil_y - (coil_width * (N - 1))) * N * pcb_layers

# sum_i=1^N 2 * ((x - w * (i-1)) + (y - w * (i-1)))
# sum_i=1^N 2 * (x + y - 2 * w * (i-1))
# sum_i=1^N (2x + 2y - 4w * (i-1))
# 2Nx + 2Ny - 4w * (sum_i=1^N i - N)
# 2N (x + y) - 4w * (N * (N+1) / 2 - N)
# 2N (x + y) - 4Nw * ((N+1) / 2 - 1)
# 2N (x + y - 2w * (N - 1) / 2)
# 2N (x + y - w * (N - 1))

R = ρ * coil_length / (trace_width * trace_thickness) # coil resistance

@objective(model, Max, I * A) # - I^2 * R +

@constraint(model, I^2 * R <= power_limit) # .. Watts - Requirement?
@constraint(model, N * coil_width <= max_coil_width)
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
println("I (A): ", I)
println("Trace width (m): ", trace_width)
println("Gap width (m): ", gap_width)


coil_width = trace_width + gap_width

A = ((coil_x * coil_y) - (coil_x * coil_width * (N - 1) / 2) - (coil_y * coil_width * (N - 1) / 2) + (coil_width * coil_width * (2 * N * N - 3 * N + 1) / 6)) * N * pcb_layers
coil_length = 2 * (coil_x + coil_y - (coil_width * (N - 1))) * N * pcb_layers

R = ρ * coil_length / (trace_width * trace_thickness)

println("Total coils width (m): ", N * coil_width)
println("Total coil length (m): ", coil_length)
println("Sum Cross section area (m^2): ", A)
println("Magnetic moment (A m^2): ", I * A)
println("Resistance (Ohm): ", R)
println("Resistive Loss (W): ", I^2 * R)