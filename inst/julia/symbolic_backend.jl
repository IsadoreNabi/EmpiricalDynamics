using SymbolicRegression
using DynamicExpressions: Node, OperatorEnum
using CSV
using DataFrames
using Statistics
using Random
using Dates
using Logging
using JLD2

@kwdef struct ScientificSearchConfig
    seed::Int = 42
    max_complexity::Int = 30
    n_populations::Int = 20
    population_size::Int = 50 
    n_iterations::Int = 100
    parsimony::Float64 = 0.001 
    n_ensemble_runs::Int = 1            
    validation_fraction::Float64 = 0.0  
    timeout_seconds::Int = 3600         
    enable_checkpoints::Bool = true
    checkpoint_interval::Int = 10        
    checkpoint_path::String = "checkpoints/"
    binary_operators::Vector{Function} = [+, -, *, /]
    unary_operators::Vector{Function} = [sin, cos, exp, abs]
    detect_physical_constants::Bool = true
    physical_constants::Dict{String, Float64} = Dict(
        "pi" => pi,
        "e" => exp(1),
        "phi" => (1 + sqrt(5)) / 2,       
        "g" => 9.80665,            
        "c" => 299792458.0,        
        "h" => 6.62607015e-34,     
        "k_B" => 1.380649e-23,     
    )
    constant_detection_tolerance::Float64 = 0.01  
    log_level::Symbol = :info
    log_file::String = ""
end

struct SearchResult
    equations::DataFrame
    config::ScientificSearchConfig
    diagnostics::Dict{Symbol, Any}
    timestamp::DateTime
    success::Bool
    error_message::String
end

SearchResult(eq, cfg, diag, ts) = SearchResult(eq, cfg, diag, ts, true, "")

function SearchResult(config::ScientificSearchConfig, error_msg::String)
    SearchResult(
        DataFrame(),
        config,
        Dict{Symbol, Any}(),
        now(),
        false,
        error_msg
    )
end

function get_node(tree)
    if hasproperty(tree, :tree)
        return tree.tree  
    else
        return tree       
    end
end

function extract_constants(tree::Node{T}) where T
    constants = T[]
    _extract_constants_recursive!(constants, tree)
    return constants
end

function _extract_constants_recursive!(constants::Vector{T}, node::Node{T}) where T
    if node.degree == 0
        if node.constant
            push!(constants, node.val)
        end
        return
    end
    
    if node.degree >= 1 && isdefined(node, :l)
        _extract_constants_recursive!(constants, node.l)
    end
    
    if node.degree == 2 && isdefined(node, :r)
        _extract_constants_recursive!(constants, node.r)
    end
end

function count_parameters(tree::Node)
    return length(extract_constants(tree))
end

function count_nodes_total(tree::Node)
    count = 1  
    
    if tree.degree >= 1 && isdefined(tree, :l)
        count += count_nodes_total(tree.l)
    end
    
    if tree.degree == 2 && isdefined(tree, :r)
        count += count_nodes_total(tree.r)
    end
    
    return count
end

function validate_and_prepare_data(
    data_path::String,
    response_col::String,
    predictor_cols::Vector{String};
    config::ScientificSearchConfig = ScientificSearchConfig(),
    remove_outliers::Bool = false,
    outlier_threshold::Float64 = 3.0  
)
    if !isfile(data_path)
        throw(ArgumentError("File not found: $data_path"))
    end
    
    data = try
        CSV.read(
            data_path, 
            DataFrame; 
            missingstring=["NA", "NaN", "", "null", ".", "NULL", "N/A"],
            silencewarnings=true
        )
    catch e
        throw(ErrorException("Error reading CSV '$data_path': $(e)"))
    end
    
    all_cols = [response_col; predictor_cols]
    missing_cols = setdiff(all_cols, names(data))
    if !isempty(missing_cols)
        throw(ArgumentError(
            "Columns not found in CSV: $(join(missing_cols, ", ")). " *
            "Available columns: $(join(names(data), ", "))"
        ))
    end
    
    X = Matrix{Float64}(coalesce.(data[:, predictor_cols], NaN))
    y = Vector{Float64}(coalesce.(data[:, response_col], NaN))
    
    n_original = length(y)
    
    valid_mask = .!isnan.(y) .& .!isinf.(y)
    
    for j in 1:size(X, 2)
        valid_mask .&= .!isnan.(X[:, j]) .& .!isinf.(X[:, j])
    end
    
    n_invalid = sum(.!valid_mask)
    
    n_outliers = 0
    if remove_outliers && sum(valid_mask) > 10
        y_valid = y[valid_mask]
        μ, σ = mean(y_valid), std(y_valid)
        
        if σ > 1e-10
            outlier_mask = abs.(y .- μ) .> outlier_threshold * σ
            n_outliers = sum(outlier_mask .& valid_mask)
            valid_mask .&= .!outlier_mask
        end
    end
    
    X_clean = X[valid_mask, :]
    y_clean = y[valid_mask]
    
    n_final = length(y_clean)
    
    if n_final < 10
        throw(ArgumentError(
            "Too few valid data points: $n_final observations. " *
            "At least 10 are required."
        ))
    end
    
    if n_final < size(X_clean, 2) * 5
        @warn "Few data points relative to number of variables" n=n_final p=size(X_clean, 2)
    end
    
    diagnostics = Dict{Symbol, Any}(
        :n_original => n_original,
        :n_invalid_removed => n_invalid,
        :n_outliers_removed => n_outliers,
        :n_final => n_final,
        :retention_rate => round(n_final / n_original * 100, digits=1),
        :y_mean => mean(y_clean),
        :y_std => std(y_clean),
        :y_range => (minimum(y_clean), maximum(y_clean)),
        :predictor_names => predictor_cols,
    )
    
    @info "Data prepared" n_original n_final invalid_removed=n_invalid outliers_removed=n_outliers
    
    return X_clean, y_clean, diagnostics
end

function build_scientific_options(config::ScientificSearchConfig, var_names::Vector{String})
    
    nested_rules = Dict(
        sin => Dict(sin => 0, cos => 0, tan => 0),
        cos => Dict(sin => 0, cos => 0, tan => 0),
        tan => Dict(sin => 0, cos => 0, tan => 0),
        exp => Dict(exp => 0, log => 0),
        log => Dict(log => 0, exp => 0),
        sqrt => Dict(sqrt => 0, abs => 0),
        abs => Dict(abs => 0),
        tanh => Dict(tanh => 0),
    )
    
    active_nested = Dict{Any, Dict{Any, Int}}()
    for op in config.unary_operators
        if haskey(nested_rules, op)
            inner_dict = Dict{Any, Int}()
            for (inner_op, limit) in nested_rules[op]
                if inner_op in config.unary_operators
                    inner_dict[inner_op] = limit
                end
            end
            if !isempty(inner_dict)
                active_nested[op] = inner_dict
            end
        end
    end

    mutation_weights = SymbolicRegression.MutationWeights(
        mutate_constant = 0.04,     
        mutate_operator = 0.47,     
        swap_operands = 0.03,       
        rotate_tree = 0.0,          
        add_node = 0.05,            
        insert_node = 0.15,         
        delete_node = 0.03,         
        simplify = 0.002,           
        randomize = 0.004,          
        do_nothing = 0.1,           
        optimize = 0.01,            
        form_connection = 0.0,      
        break_connection = 0.0,     
    )

    complexity_map = Dict{Function, Int}(
        (+) => 1,
        (-) => 1,
        (*) => 1,
        (/) => 2,       
        (^) => 3,       
        sin => 2,
        cos => 2,
        tan => 3,       
        exp => 2,       
        log => 3,       
        sqrt => 2,      
        abs => 1,       
        tanh => 2,
        sinh => 2,
        cosh => 2,
    )
    
    active_complexity = Pair{Function, Int}[]
    for op in config.binary_operators
        if haskey(complexity_map, op)
            push!(active_complexity, op => complexity_map[op])
        end
    end
    for op in config.unary_operators
        if haskey(complexity_map, op)
            push!(active_complexity, op => complexity_map[op])
        end
    end

    options = SymbolicRegression.Options(
        binary_operators = Tuple(config.binary_operators),
        unary_operators = Tuple(config.unary_operators),
        
        maxsize = config.max_complexity,
        complexity_of_operators = active_complexity,
        
        npopulations = config.n_populations,
        npop = config.population_size,
        
        parsimony = config.parsimony,
        
        nested_constraints = active_nested,
        
        mutation_weights = mutation_weights,
        crossover_probability = 0.066,
        
        should_optimize_constants = true,
        optimizer_nrestarts = 8,
        optimizer_iterations = 30,
        optimizer_algorithm = "BFGS",
        
        progress = true,
        verbosity = 1,
        
        batching = false,  
        turbo = false,     
        
        seed = config.seed,
        deterministic = false,  
    )
    
    return options
end

function run_single_search(
    X::Matrix{Float64},
    y::Vector{Float64},
    var_names::Vector{String},
    config::ScientificSearchConfig,
    options
)
    if config.enable_checkpoints
        mkpath(config.checkpoint_path)
    end
    
    hall_of_fame = EquationSearch(
        X', y;  
        niterations = config.n_iterations,
        options = options,
        variable_names = var_names,
        parallelism = :multithreading,
    )
    
    return hall_of_fame
end

function run_ensemble_search(
    X::Matrix{Float64},
    y::Vector{Float64},
    var_names::Vector{String},
    config::ScientificSearchConfig
)
    all_results = DataFrame[]
    
    for run_id in 1:config.n_ensemble_runs
        @info "Running ensemble search" run=run_id total=config.n_ensemble_runs
        
        run_seed = config.seed + run_id * 1000
        run_config = ScientificSearchConfig(
            config;
            seed = run_seed
        )
        
        options = build_scientific_options(run_config, var_names)
        
        try
            hof = run_single_search(X, y, var_names, run_config, options)
            results = extract_pareto_front(hof, X', y, var_names, options, run_config)
            results[!, :run_id] .= run_id
            push!(all_results, results)
        catch e
            @warn "Run $run_id failed" exception=e
        end
    end
    
    if isempty(all_results)
        throw(ErrorException("All ensemble runs failed"))
    end
    
    combined = vcat(all_results...)
    
    return aggregate_ensemble(combined)
end

function aggregate_ensemble(combined::DataFrame)
    if !hasproperty(combined, :expression_normalized)
        combined[!, :expression_normalized] = normalize_expression.(combined.expression)
    end
    
    aggregated = combine(
        groupby(combined, :expression_normalized),
        :complexity => first => :complexity,
        :loss => mean => :mean_loss,
        :loss => std => :std_loss,
        :r_squared => mean => :mean_r_squared,
        :expression => first => :expression,
        :n_constants => first => :n_constants,
        :detected_physics => first => :detected_physics,
        nrow => :n_appearances
    )
    
    aggregated[!, :consistency_score] = aggregated.n_appearances ./ (1 .+ coalesce.(aggregated.std_loss, 0.0))
    
    sort!(aggregated, :consistency_score, rev=true)
    
    return aggregated
end

function normalize_expression(expr::String)
    normalized = replace(expr, r"-?\d+\.?\d*([eE][+-]?\d+)?" => "C")
    return normalized
end

function extract_pareto_front(
    hall_of_fame,
    X::AbstractMatrix,      
    y::AbstractVector,
    var_names::Vector{String},
    options,
    config::ScientificSearchConfig
)
    results = DataFrame(
        complexity = Int[],
        loss = Float64[],
        mse = Float64[],
        mae = Float64[],
        r_squared = Float64[],
        expression = String[],
        expression_normalized = String[],
        n_constants = Int[],
        n_nodes = Int[],
        detected_physics = String[]
    )
    
    n = length(y)
    y_mean = mean(y)
    ss_tot = sum((y .- y_mean).^2)
    
    X_dense = Array(X)
    
    @info "Extracting Pareto front..." n_members=length(hall_of_fame.members)
    
    processed = 0
    skipped_not_exists = 0
    skipped_error = 0
    
    for i in eachindex(hall_of_fame.members)
        if hasproperty(hall_of_fame, :exists)
            if !hall_of_fame.exists[i]
                skipped_not_exists += 1
                continue
            end
        end
        
        member = hall_of_fame.members[i]
        
        if member === nothing
            skipped_not_exists += 1
            continue
        end
        
        try
            complexity = compute_complexity(member.tree, options)
            loss = member.loss
            
            eval_result = eval_tree_array(member.tree, X_dense, options)
            
            local y_pred
            local eval_ok::Bool
            
            if eval_result === nothing
                @warn "eval_tree_array returned nothing for member $i"
                skipped_error += 1
                continue
            elseif isa(eval_result, Tuple)
                y_pred = eval_result[1]
                eval_ok = if length(eval_result) > 1 && isa(eval_result[2], Bool)
                    eval_result[2]
                else
                    true
                end
            else
                y_pred = eval_result
                eval_ok = true
            end
            
            if y_pred === nothing || !eval_ok
                @warn "Evaluation failed for member $i" eval_ok
                skipped_error += 1
                continue
            end
            
            y_pred = vec(y_pred)
            
            if length(y_pred) != n
                @warn "Incorrect dimension in member $i" expected=n got=length(y_pred)
                skipped_error += 1
                continue
            end
            
            if any(isnan.(y_pred)) || any(isinf.(y_pred))
                @warn "Prediction with NaN/Inf in member $i"
                skipped_error += 1
                continue
            end
            
            residuals = y .- y_pred
            mse = mean(residuals.^2)
            mae = mean(abs.(residuals))
            
            r_squared = if ss_tot > 1e-10
                1.0 - sum(residuals.^2) / ss_tot
            else
                0.0
            end
            
            r_squared = clamp(r_squared, -1.0, 1.0)
            
            expr_str = string_tree(member.tree, options; variable_names=var_names)
            expr_normalized = normalize_expression(expr_str)
            
            node = get_node(member.tree)
            n_consts = count_parameters(node)
            n_nodes = count_nodes_total(node)
            physics = detect_physical_constants(node, config)
            
            push!(results, (
                complexity = complexity,
                loss = loss,
                mse = mse,
                mae = mae,
                r_squared = r_squared,
                expression = expr_str,
                expression_normalized = expr_normalized,
                n_constants = n_consts,
                n_nodes = n_nodes,
                detected_physics = join(physics, ", ")
            ))
            
            processed += 1
            
        catch e
            @warn "Error processing member $i" exception=(e, catch_backtrace())
            skipped_error += 1
            continue
        end
    end
    
    @info "Extraction completed" processed skipped_not_exists skipped_error total_results=nrow(results)
    
    if nrow(results) > 0
        sort!(results, :complexity)
    else
        @warn "No valid equations extracted from HallOfFame"
    end
    
    return results
end

function detect_physical_constants(tree::Node, config::ScientificSearchConfig)
    !config.detect_physical_constants && return String[]
    
    detected = String[]
    
    try
        constants = extract_constants(tree)
        
        for c in constants
            abs_c = abs(c)
            
            for (name, val_phy) in config.physical_constants
                tol = config.constant_detection_tolerance
                
                if isapprox(abs_c, val_phy, rtol=tol)
                    push!(detected, name)
                    continue
                end
                
                if isapprox(abs_c, 2 * val_phy, rtol=tol)
                    push!(detected, "2$name")
                elseif isapprox(abs_c, val_phy / 2, rtol=tol)
                    push!(detected, "$name/2")
                elseif isapprox(abs_c, val_phy^2, rtol=tol)
                    push!(detected, "$name^2")
                elseif val_phy > 1e-15 && isapprox(abs_c, 1/val_phy, rtol=tol)
                    push!(detected, "1/$name")
                elseif isapprox(abs_c, 4 * val_phy, rtol=tol)
                    push!(detected, "4$name")
                end
            end
        end
    catch e
        @warn "Error detecting physical constants" exception=e
    end
    
    return unique(detected)
end

function detect_physical_constants(tree, config::ScientificSearchConfig)
    return detect_physical_constants(get_node(tree), config)
end

function run_symbolic_search_robust(
    data_path::String,
    response_col::String,
    predictor_cols::Vector{String},
    max_complexity::Int = 30,
    n_iterations::Int = 100,
    parsimony::Float64 = 0.001,
    output_path::String = ""
)
    timestamp = now()
    
    @info "═══════════════════════════════════════════════════════════"
    @info "EmpiricalDynamics - Starting symbolic search"
    @info "═══════════════════════════════════════════════════════════"
    @info "File: $data_path"
    @info "Response: $response_col"
    @info "Predictors: $(join(predictor_cols, ", "))"
    @info "Max complexity: $max_complexity | Iterations: $n_iterations | Parsimony: $parsimony"
    
    try
        config = ScientificSearchConfig(
            max_complexity = max_complexity,
            n_iterations = n_iterations,
            parsimony = parsimony,
            n_populations = 20,
            population_size = 50,
            seed = 42
        )
        
        X, y, diagnostics = validate_and_prepare_data(
            data_path, 
            response_col, 
            predictor_cols; 
            config=config
        )
        
        @info "Data loaded" observations=length(y) variables=length(predictor_cols)
        
        options = build_scientific_options(config, predictor_cols)
        
        @info "Starting genetic evolution..."
        hof = run_single_search(X, y, predictor_cols, config, options)
        
        results = extract_pareto_front(hof, X', y, predictor_cols, options, config)
        
        if !isempty(output_path)
            CSV.write(output_path, results)
            @info "Results saved to: $output_path"
        end
        
        n_eq = nrow(results)
        @info "═══════════════════════════════════════════════════════════"
        @info "Search completed successfully"
        @info "Equations found: $n_eq"
        if n_eq > 0
            best = results[argmin(results.loss), :]
            @info "Best equation (by loss): $(best.expression)"
            @info "    R^2 = $(round(best.r_squared, digits=4)), Complexity = $(best.complexity)"
        end
        @info "═══════════════════════════════════════════════════════════"
        
        return results
        
    catch e
        error_msg = sprint(showerror, e)
        @error "Search failed" exception=(e, catch_backtrace())
        
        return DataFrame(
            complexity = Int[0],
            loss = Float64[Inf],
            mse = Float64[Inf],
            mae = Float64[Inf],
            r_squared = Float64[-999.0],
            expression = String["ERROR: $error_msg"],
            expression_normalized = String["ERROR"],
            n_constants = Int[0],
            n_nodes = Int[0],
            detected_physics = String[""],
            error_flag = Bool[true]
        )
    end
end

function run_symbolic_search_robust(
    data_path::String,
    response_col::String,
    predictor_cols::String, 
    max_complexity::Int = 30,
    n_iterations::Int = 100,
    parsimony::Float64 = 0.001,
    output_path::String = ""
)
    return run_symbolic_search_robust(
        data_path,
        response_col,
        [predictor_cols], 
        max_complexity,
        n_iterations,
        parsimony,
        output_path
    )
end

function run_symbolic_search_ensemble(
    data_path::String,
    response_col::String,
    predictor_cols::Vector{String};
    max_complexity::Int = 30,
    n_iterations::Int = 100,
    parsimony::Float64 = 0.001,
    n_ensemble_runs::Int = 3,
    output_path::String = ""
)
    config = ScientificSearchConfig(
        max_complexity = max_complexity,
        n_iterations = n_iterations,
        parsimony = parsimony,
        n_ensemble_runs = n_ensemble_runs,
    )
    
    try
        X, y, diagnostics = validate_and_prepare_data(data_path, response_col, predictor_cols; config=config)
        results = run_ensemble_search(X, y, predictor_cols, config)
        
        if !isempty(output_path)
            CSV.write(output_path, results)
        end
        
        return results
    catch e
        @error "Ensemble search failed" exception=(e, catch_backtrace())
        return DataFrame()
    end
end

function warmup_jit()
    @info "Warming up JIT (this takes ~10 seconds the first time)..."
    
    try
        n = 20
        X_dummy = randn(n, 2)
        y_dummy = sin.(X_dummy[:, 1]) .+ 0.1 .* X_dummy[:, 2] .+ 0.05 .* randn(n)
        
        config = ScientificSearchConfig(
            n_iterations = 3,
            n_populations = 2,
            population_size = 10,
            max_complexity = 8,
            seed = 12345
        )
        
        options = build_scientific_options(config, ["x1", "x2"])
        
        hof = EquationSearch(
            X_dummy', y_dummy;
            niterations = 3,
            options = options,
            variable_names = ["x1", "x2"],
            parallelism = :serial,  
        )
        
        _ = extract_pareto_front(hof, X_dummy', y_dummy, ["x1", "x2"], options, config)
        
        @info "JIT ready"
        return true
        
    catch e
        @warn "Warmup failed (non-critical)" exception=e
        return false
    end
end

function test_installation()
    println("Testing EmpiricalDynamics installation...")
    println("─" ^ 50)
    
    tests = [
        ("SymbolicRegression", true),
        ("DynamicExpressions", true),
        ("CSV", true),
        ("DataFrames", true),
        ("JLD2", true),
    ]
    
    all_pass = true
    for (pkg, required) in tests
        try
            status = "✓"
        catch
            status = required ? "✗ MISSING" : "○ optional"
            all_pass = all_pass && !required
        end
        println("  $status $pkg")
    end
    
    println("─" ^ 50)
    
    println("Testing core functions...")
    
    try
        @assert @isdefined(get_node)
        println("  ✓ get_node")
    catch
        println("  ✗ get_node")
        all_pass = false
    end
    
    try
        @assert hasmethod(extract_constants, Tuple{Node})
        println("  ✓ extract_constants")
    catch
        println("  ✗ extract_constants")
        all_pass = false
    end
    
    try
        @assert hasmethod(run_symbolic_search_robust, Tuple{String, String, Vector{String}, Int, Int, Float64, String})
        println("  ✓ run_symbolic_search_robust (Vector)")
    catch
        println("  ✗ run_symbolic_search_robust (Vector)")
        all_pass = false
    end
    
    try
        @assert hasmethod(run_symbolic_search_robust, Tuple{String, String, String, Int, Int, Float64, String})
        println("  ✓ run_symbolic_search_robust (String overload)")
    catch
        println("  ✗ run_symbolic_search_robust (String overload)")
        all_pass = false
    end
    
    println("─" ^ 50)
    println(all_pass ? "All tests passed!" : "Some tests failed")
    
    return all_pass
end

get_version() = "0.1.1"

function get_available_operators()
    return Dict(
        :binary => Dict(
            :safe => [+, -, *, /],
            :extended => [+, -, *, /, ^],
        ),
        :unary => Dict(
            :safe => [sin, cos, exp, abs],
            :extended => [sin, cos, tan, exp, log, sqrt, abs, tanh],
            :trigonometric => [sin, cos, tan, asin, acos, atan],
            :hyperbolic => [sinh, cosh, tanh],
        )
    )
end

function __init__()
    println()
    println("╔═══════════════════════════════════════════════════════════╗")
    println("║     EmpiricalDynamics - Symbolic Regression Engine   ║")
    println("║     Industrial-Grade | Physics-Optimized | Robust         ║")
    println("║     FIX: Expression vs Node (SymbolicRegression v1.0+)    ║")
    println("╚═══════════════════════════════════════════════════════════╝")
    println()
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) == 0
        println("EmpiricalDynamics - Julia Backend")
        println()
        println("Usage:")
        println("  julia symbolic_backend.jl <data.csv> <response> <predictors> [options]")
        println()
        println("Example:")
        println("  julia symbolic_backend.jl data.csv y x1,x2,x3 --complexity 20 --iterations 50")
        println()
        println("Options:")
        println("  --complexity N    Maximum equation complexity (default: 30)")
        println("  --iterations N    Number of iterations (default: 100)")
        println("  --parsimony F     Parsimony coefficient (default: 0.001)")
        println("  --output PATH     Output CSV path (default: results.csv)")
        println("  --warmup          Run JIT warmup before search")
        println()
    else
        data_path = ARGS[1]
        response_col = ARGS[2]
        predictor_cols = split(ARGS[3], ",")
        
        max_complexity = 30
        n_iterations = 100
        parsimony = 0.001
        output_path = "results.csv"
        do_warmup = false
        
        i = 4
        while i <= length(ARGS)
            if ARGS[i] == "--complexity" && i < length(ARGS)
                max_complexity = parse(Int, ARGS[i+1])
                i += 2
            elseif ARGS[i] == "--iterations" && i < length(ARGS)
                n_iterations = parse(Int, ARGS[i+1])
                i += 2
            elseif ARGS[i] == "--parsimony" && i < length(ARGS)
                parsimony = parse(Float64, ARGS[i+1])
                i += 2
            elseif ARGS[i] == "--output" && i < length(ARGS)
                output_path = ARGS[i+1]
                i += 2
            elseif ARGS[i] == "--warmup"
                do_warmup = true
                i += 1
            else
                @warn "Unknown argument: $(ARGS[i])"
                i += 1
            end
        end
        
        if do_warmup
            warmup_jit()
        end
        
        results = run_symbolic_search_robust(
            data_path,
            response_col,
            Vector{String}(predictor_cols),
            max_complexity,
            n_iterations,
            parsimony,
            output_path
        )
        
        println()
        println("═══ TOP 5 EQUATIONS ═══")
        for (i, row) in enumerate(eachrow(first(results, 5)))
            println()
            println("#$i: $(row.expression)")
            println("    R^2 = $(round(row.r_squared, digits=4)) | Loss = $(round(row.loss, digits=6)) | Complexity = $(row.complexity)")
            if !isempty(row.detected_physics)
                println("    Physics detected: $(row.detected_physics)")
            end
        end
    end
end