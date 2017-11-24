#using Base.Test
using ArgParse
using JLD

include("src/Tightbinding.jl")
include("src/OnebandModel.jl")

include("src/PairingModel.jl")
include("src/NaivePairingModel.jl")

include("oneband.jl")

function parse_commandline()
  argsetting = ArgParseSettings("Oneband One Dimension",
                                version="Version 0.1",
                                add_version=true)

  @add_arg_table argsetting begin
    "--nx"
      arg_type = Int
      required = true
    "--ny"
      arg_type = Int
      required = true
    "--xi0"
      arg_type = Float64
      nargs = '+'
    "--xi1"
      arg_type = Float64
      nargs = '+'
    "--pairing_d"
      arg_type = Float64
      nargs = '+'
    "--pairing_s"
      arg_type = Float64
      nargs = '+'
    "--pairing_p"
      arg_type = Float64
      nargs = '+'
    "--outfile", "-o"
      arg_type = String
      required = true
      help = "output filename"
  end
  parsed_args = parse_args(ARGS, argsetting)
  return parsed_args
end


function main()
  args = parse_commandline()

  const n_nambu = 2
  const n_basis = 2
  const n_orbital = 1

  result = Dict()

  let 
    nx = args["nx"]
    ny = args["ny"]

    count = length(Base.product(args["pairing_d"], args["pairing_s"], args["pairing_p"], args["xi0"], args["xi1"]))
    for (idx, (Δd, Δs, Δp, ξ₀, ξ₁)) in enumerate(Base.product(args["pairing_d"], args["pairing_s"], args["pairing_p"], args["xi0"], args["xi1"]))
      r = compute_mixedspace(nx, ny, Δd, Δs, Δp, ξ₀, ξ₁)
      n = length(r["eigenvalues"])
      ψ = reshape(r["eigenvectors"], (n_nambu, n_orbital, n_basis, nx, n))
      ρ = reshape( sum( (abs(ψ[1,:,:,:, :]).^2), (1, 2)), (nx, n))
      result[(Δd, Δs, Δp, ξ₀, ξ₁, "mixed")] = Dict("eigenvalues" => r["eigenvalues"],
                                                   "spectralweights" => ρ,
                                                   "momentums" => r["momentums"])
    end
  end

  save(args["outfile"],
       "type", :ldos,
       "parameters", args,
       "result", result)
end

main()

