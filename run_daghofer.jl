#using Base.Test
using ArgParse
using JLD

include("src/Tightbinding.jl")
include("src/DaghoferModel.jl")

include("src/PairingModel.jl")
include("src/NaivePairingModel.jl")

include("daghofer.jl")

function parse_commandline()
  argsetting = ArgParseSettings("Daghofer One Dimension",
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
      required = true
    "--xi1"
      arg_type = Float64
      required = true
    "--pairing_d"
      arg_type = Float64
      required = true
    "--pairing_s"
      arg_type = Float64
      required = true
    "--pairing_p"
      arg_type = Float64
      required = true
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

  result = Dict{String, Any}()
  result["momentum"] = compute_momentumspace(args)
  result["mixed"] = compute_mixedspace(args)

  save(args["outfile"],
       "parameters", args,
       "mixed", result["mixed"],
       "momentum", result["momentum"])
end

main()

