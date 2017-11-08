using JLD

include("src/Tightbinding.jl")
include("src/DaghoferModel.jl")

include("src/PairingModel.jl")
include("src/NaivePairingModel.jl")


const daghofer_parameter = DaghoferModel.DaghoferParameter(
        0.02,  0.06,  0.03, -0.01,  0.20,
        0.30, -0.20,  0.10,  0.40,  0.20)



dm = DaghoferModel.daghoferModel(Complex128, daghofer_parameter)

h = Tightbinding.make_realspace_dense(dm, (16,16); periodic=true)
@show real(h)
