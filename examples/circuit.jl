module YaoQASMReader
using Yao
using Yao.EasyBuild: FSimGate
const GATE_TENSORS = Dict{Symbol, AbstractBlock}()

GATE_TENSORS[:i] = I2
GATE_TENSORS[:x] = X
GATE_TENSORS[:y] = Y
GATE_TENSORS[:z] = Z

@const_gate X12 = [1 -1im; -1im 1]/sqrt(2)
GATE_TENSORS[:x_1_2] = X12
@const_gate Y12 = [1 -1; 1 1]/sqrt(2)
GATE_TENSORS[:y_1_2] = Y12
@const_gate HZ12 = [1 1+1im; 1-1im 1]/2
GATE_TENSORS[:hz_1_2] = HZ12

GATE_TENSORS[:h] = H
GATE_TENSORS[:s] = ConstGate.S
GATE_TENSORS[:t] = ConstGate.T

GATE_TENSORS[:cx] = control(2,1,2=>X)
GATE_TENSORS[:cz] = control(2,1,2=>Z)

rz(theta) = Rz(theta)
fsim(theta, phi) = FSimGate(theta, phi)

"""
    yaocircuit_from_qasm(qasm_str::String)

Function to create the yao circuit described in the given qasm file.
"""
function yaocircuit_from_qasm(qasm_file::String)
    qasm_string = open(f -> read(f, String), qasm_file)
    qasm_lines = split(qasm_string, "\n")
    num_qubits = parse(Int, qasm_lines[1])
    qubit_map = Dict{Int, Int}(-1 => 0)
    c = chain(num_qubits)

    for line in qasm_lines[2:end]
        words = filter(x->length(x)>0, split(line, ['(', ')', ' ', ',']))
        if length(words) > 0
            targets, tensor = parse_words(words, qubit_map)
            push!(c, put(num_qubits, (targets...,)=>tensor))
        end
    end
    c
end


"""
Map the given target qubits to indices according to qubit_map

If the target qubit isn't mapped to anything, increment the total number of qubits and map to total number.
"""
function transform_targets!(targets, qubit_map)
    for i = 1:length(targets)
        if haskey(qubit_map, targets[i])
            targets[i] = qubit_map[targets[i]]
        else
            qubit_map[-1] += 1
            qubit_map[targets[i]] = qubit_map[-1]
            targets[i] = qubit_map[-1]
        end
    end
    targets
end

function parse_words(words, qubit_map)
    if words[2] == "fsim"
        theta = parse(Float64, words[3])
        phi = parse(Float64, words[4])
        tensor = fsim(theta, phi)
        targets = parse.(Int, words[5:6])

    elseif words[2] == "rz"
        theta = parse(Float64, words[3])
        tensor = rz(theta)
        targets = parse.(Int, words[4:end])

    else
        gate_symbol = Symbol(words[2])
        targets = parse.(Int, words[3:end])
        tensor = GATE_TENSORS[gate_symbol]
    end

    transform_targets!(targets, qubit_map)
    targets, tensor
end
end


using Yao
# Load functions to read qflex qasm files.
using .YaoQASMReader: yaocircuit_from_qasm

# circuit source: https://github.com/brenjohn/Contraction-Order-Bench/tree/main/data/circuits
cirq_name = "sycamore_53_20_0"
@info("running circuit: $(cirq_name)")

# Create the TensorNetworkCircuit object for the circuit
qasm_file = joinpath(@__DIR__, cirq_name * ".txt")
c = yaocircuit_from_qasm(qasm_file)
using Random; Random.seed!(2)
time_elapsed = @elapsed net = yao2einsum(c, initial_state=Dict(zip(1:nqubits(c), zeros(Int,nqubits(c)))), final_state=Dict(zip(1:nqubits(c), zeros(Int,nqubits(c)))), optimizer=TreeSA(ntrials=10, niters=20, sc_target=52))
@info "contraction complexity: $(contraction_complexity(net)), time cost: $(time_elapsed)s"