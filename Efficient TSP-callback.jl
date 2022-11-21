using Graphs
using JuMP
using GLPK
using BenchmarkTools
using Distances
using Plots
using Gurobi
using Test

struct Delivery
    x::Float64
    y::Float64
end

function random_delivery()
    x = rand(0.0:0.1:100.0)
    y = rand(0.0:0.1:100.0)
    return Delivery(x, y)
end

struct TSProblem
    deliveries::Array{Delivery}
end

function random_instance(n_deliveries)
    deliveries = [random_delivery() for _=1:n_deliveries]
    problem = TSProblem(deliveries)
    return problem
end

function dist(del1::Delivery, del2::Delivery)
    return sqrt((del1.x - del2.x)^2 + (del1.y - del2.y)^2)
end

function calc_travelmatrix(deliveries::Array{Delivery})
    tm = zeros(Float64, length(deliveries), length(deliveries))
    for i = 1:length(deliveries)
        for j = 1:length(deliveries)
            tm[i, j] = dist(deliveries[i], deliveries[j])
        end
    end
    return tm
end

test_m = [0 1 0 0; 
          1 0 0 0;
          0 0 0 1;
          0 0 1 0;]
"""
given an incidence matrix of a tsp route return nodes of the shortest cycle
"""
function shortest_subtour(matrix::Matrix{Int64})
    g = Graphs.SimpleDiGraph(matrix)
    cycles = []
    max_cycle_len = Graphs.nv(g)
    for node in Graphs.vertices(g)
        push!(cycles, Graphs.neighborhood(g, node, max_cycle_len))
    end
    cycles = filter(x -> max_cycle_len > length(x) > 1, cycles)
    if isempty(cycles)
        return []
    end
    return sort(cycles, by=length)[1]
end
#= @test shortest_subtour(test_m) == [1, 2]

test_m_empty = [0 0 0 0; 
          0 0 0 0;
          0 0 0 0;
          0 0 0 0;]
@test shortest_subtour(test_m_empty) == []

test_m_no_subtours = 
         [0 1 0 0; 
          0 0 1 0;
          0 0 0 1;
          1 0 0 0;]
@test shortest_subtour(test_m_no_subtours) == []

test_subtours = 
    [0 0 1 0 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     0 1 0 0 0 0;
     0 0 0 1 0 0;
     1 0 0 0 0 0]
@test shortest_subtour(test_subtours) == [1, 3, 6]

test_subtours2 = 
   [0 0 1 0 0; 
    1 0 0 0 0;
    0 0 0 0 1;
    0 1 0 0 0;
    0 0 0 1 0]
@test shortest_subtour(test_subtours2) == [] =#









function display_solution(problem, route)
    x_pos = [c.x for c in problem.deliveries]
    y_pos = [c.y for c in problem.deliveries]
    plot_result = scatter(x_pos, y_pos, shape = :circle, markersize = 6)
    for i in 1:length(problem.deliveries)
        for j in 1:length(problem.deliveries)
            val = route[i, j]
            if val > 0
                del1 = problem.deliveries[i]
                del2 = problem.deliveries[j]
                plot!([del1.x, del2.x], [del1.y, del2.y], legend = false)
            end
        end
    end
    return plot_result
end;






#= function solve_tsp(deliveries::Int64, solver, show_viz=false)
    problem=random_instance(deliveries)
    travelmatrix = calc_travelmatrix(problem.deliveries)
    model = Model(solver)
    set_silent(model)
    # route is an adjence matrix representing a route traveled
    route=@variable(model, route[1:length(problem.deliveries), 1:length(problem.deliveries)], Bin)

    # ensure all events are planned
    @constraint(model, [i = 1:length(problem.deliveries)], sum(route[i, :]) == 1.0)
    @constraint(model, [c = 1:length(problem.deliveries)], sum(route[:, c]) == 1.0)
    # disallow traveling to itself
    @constraint(model, [j = 1:length(problem.deliveries)], route[j, j] == 0)

    traveltime = travelmatrix.* route 
    @objective(model, Min, sum(traveltime))

    function callback(cb_data)
        status = callback_node_status(cb_data, model)
        if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
            return
        end
        x_val = callback_value.(cb_data, route)
        x_val = round.(x_val)
        x_val = Int64.(x_val)
        cycle = shortest_subtour(x_val)
        sub_inds = [(i, j) for (i, j) in Iterators.product(cycle, cycle) if i != j]
        if length(sub_inds) > 0
            con = @build_constraint(sum(route[i, j] for (i,j) in sub_inds) <= length(cycle) -1 )
            MOI.submit(model, MOI.LazyConstraint(cb_data), con)
        end
    end
    MOI.set(model, MOI.LazyConstraintCallback(), callback)
    optimize!(model)
    route_val = JuMP.value.(route)
    if show_viz
        println(termination_status(model))
        return display_solution(problem, route_val)
    end
end

@time solve_tsp(30, GLPK.Optimizer, true) =#







###########
deliveries = 30;
solver = GLPK.Optimizer
problem=random_instance(deliveries)
travelmatrix = calc_travelmatrix(problem.deliveries)
cvrp = Model(GLPK.Optimizer)
#set_silent(model)
# route is an adjence matrix representing a route traveled
route=@variable(cvrp, route[1:length(problem.deliveries), 1:length(problem.deliveries)], Bin)

# ensure all events are planned
@constraint(cvrp, [i = 1:length(problem.deliveries)], sum(route[i, :]) == 1.0)
@constraint(cvrp, [c = 1:length(problem.deliveries)], sum(route[:, c]) == 1.0)
# disallow traveling to itself
@constraint(cvrp, [j = 1:length(problem.deliveries)], route[j, j] == 0)

traveltime = travelmatrix.* route 
@objective(cvrp, Min, sum(traveltime))

function callback(cb_data)
    status = callback_node_status(cb_data, cvrp)
    if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
        return
    end
    x_val = callback_value.(cb_data, route)
    #println("1-x_val : $(x_val)")
    x_val = round.(x_val)
    #println("2-x_val : $(x_val)")
    x_val = Int64.(x_val)
    #println("3-x_val : $(x_val)")
    cycle = shortest_subtour(x_val)
    #println("cycle : $(cycle)")
    sub_inds = [(i, j) for (i, j) in Iterators.product(cycle, cycle) if i != j]
    if length(sub_inds) > 0
        con = @build_constraint(sum(route[i, j] for (i,j) in sub_inds) <= length(cycle) -1 )
        MOI.submit(cvrp, MOI.LazyConstraint(cb_data), con)
        println("Adding cut: $(con)")
    end
end
MOI.set(cvrp, MOI.LazyConstraintCallback(), callback)
optimize!(cvrp)
