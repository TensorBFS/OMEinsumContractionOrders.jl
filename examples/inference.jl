using TensorInference, OMEinsumContractionOrders

problems = dataset_from_artifact("uai2014")["MAR"]
problem_set_name = "relational"
tamaki = [251, 3, 8, 101, 11]
for (id, problem) in problems[problem_set_name]
    @info "Testing: $(problem_set_name)_$id"
    # optimizer = TreeSA(ntrials=1, niters=5, βs=0.1:0.1:100)
    optimizer = Treewidth(; alg=ND(MF(), METISND(; ufactor=150); limit=200, level=6))

    tn = TensorNetworkModel(read_model(problem); optimizer, evidence=read_evidence(problem))
    ref_sol = read_solution(problem)
    evidence = read_evidence(problem)

    # does not optimize over open vertices
    @info "Contraction complexity: $(contraction_complexity(tn)), tamaki tw = $(tamaki[id])"
end