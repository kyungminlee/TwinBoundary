

#using Base.Test
using ArgParse
using JLD
using ProgressMeter

include("src/Tightbinding.jl")
include("src/DaghoferModel.jl")

include("src/PairingModel.jl")
include("src/NaivePairingModel.jl")

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




"""
Check whether 1D and 0D(momentumspace) calculations yield the same result.
"""
function main()

  args = parse_commandline()

  const n_nambu = 2
  const n_basis = 2
  const n_orbital = 3

  const nx, ny = args["nx"], args["ny"]

  const Δd = args["pairing_d"]
  const Δs = args["pairing_s"]
  const Δp = args["pairing_p"]

  const ξ₀ = args["xi0"]
  const ξ₁ = args["xi1"]

  const kxs = linspace(0, 2π, nx+1)[1:end-1]
  const kys = linspace(0, 2π, ny+1)[1:end-1]

  const daghofer_parameter = DaghoferModel.DaghoferParameter(
          0.02,  0.06,  0.03, -0.01,  0.20,
          0.30, -0.20,  0.10,  0.40,  0.20)

  const daghofer_model = DaghoferModel.daghoferModel(Complex128, daghofer_parameter)
  const double_daghofer_model = Tightbinding.double_unitcell(daghofer_model)


  result = Dict{String, Any}()

  let
    kxky_bookkeep = []
    all_eigenvalues_momentum = []
    all_eigenvectors_momentum = []

    double_daghofer_momentumspace = Tightbinding.make_momentumspace(double_daghofer_model)
    pairing_momentumspace_d = make_pairing_momentumspace(Δd, formfactors["d1"])
    pairing_momentumspace_s = make_pairing_momentumspace(Δs, formfactors["s1"])
    @showprogress for (ky, kx) in Base.product(kys, kxs)
      hk = zeros(Complex128, (n_nambu, n_orbital, n_basis, n_nambu, n_orbital, n_basis))

      let
        tk = double_daghofer_momentumspace(kx, ky)
        Δk = pairing_momentumspace_d(kx, ky) + pairing_momentumspace_s(kx, ky)
        hk[1, :, :, 1, :, :] = reshape( tk, (n_orbital, n_basis, n_orbital, n_basis))
        hk[2, :, :, 2, :, :] = reshape(-tk, (n_orbital, n_basis, n_orbital, n_basis))
        for i_orbital in 1:n_orbital
          hk[1, i_orbital, :, 2, i_orbital, :] = reshape(Δk, (n_basis, n_basis))
          hk[2, i_orbital, :, 1, i_orbital, :] = reshape(Δk, (n_basis, n_basis))
        end
      end

      hamiltonian_matrix = reshape(hk, (n_nambu * n_orbital * n_basis, n_nambu * n_orbital * n_basis))
      hamiltonian_matrix = 0.5 * (hamiltonian_matrix + hamiltonian_matrix')
      eivals, eivecs = eig(Hermitian(hamiltonian_matrix))
      push!(kxky_bookkeep, [(kx, ky) for e in eivals])
      push!(all_eigenvalues_momentum, eivals)
      push!(all_eigenvectors_momentum, eivecs)
    end

    all_eigenvalues_momentum = vcat(all_eigenvalues_momentum...)
    all_eigenvectors_momentum = hcat(all_eigenvectors_momentum...)
    kxky_bookkeep = vcat(kxky_bookkeep...)

    let 
      idxperm = sortperm(all_eigenvalues_momentum)
      all_eigenvalues_momentum = all_eigenvalues_momentum[idxperm]
      all_eigenvectors_momentum = all_eigenvectors_momentum[:, idxperm]
      kxky_bookkeep = kxky_bookkeep[idxperm]
    end

    result["momentum"] = Dict{String, Any}("eigenvalues" => all_eigenvalues_momentum,
                                           "eigenvectors" => all_eigenvectors_momentum,
                                           "momentums" => kxky_bookkeep)
  end


  let
    ky_bookkeep = []
    all_eigenvalues_mixed = []
    all_eigenvectors_mixed = []

    double_daghofer_mixedspace = Tightbinding.make_mixedspace(double_daghofer_model, nx; periodic=true)

    pairing_orderparameter = [(        Δs * tanh((x - 1 - nx)/(2*ξ₀))
                               + 1im * Δp * tanh((x - 1 - nx)/(2*ξ₁))) for x in 1:2*nx]
    #pairing_orderparameter = reshape(pairing_orderparameter, (2, nx))
    pairing_orderparameter = ones(Complex128, (2, nx)) * Δs

    pairing_mixedspace_uniform = make_pairing_mixedspace(Δd * ones(Complex128, (2,nx)),
                                                         formfactors["d1"];
                                                         periodic=true)

    pairing_mixedspace = make_pairing_mixedspace(pairing_orderparameter,
                                                 formfactors["s1"];
                                                 periodic=true)

    @showprogress for ky in kys
      hamiltonian = zeros(Complex128, (n_nambu, n_orbital, n_basis, nx,
                                       n_nambu, n_orbital, n_basis, nx))
      let
        tky  = double_daghofer_mixedspace( ky)
        Δky = pairing_mixedspace_uniform(ky) + pairing_mixedspace(ky)

        hamiltonian[1, :, :, :, 1, :, :, :] = reshape( tky, (n_orbital, n_basis, nx,
                                                             n_orbital, n_basis, nx))
        hamiltonian[2, :, :, :, 2, :, :, :] = reshape(-tky, (n_orbital, n_basis, nx,
                                                             n_orbital, n_basis, nx))
        for i_orbital in 1:n_orbital
          hamiltonian[1, i_orbital, :, :, 2, i_orbital, :, :] = reshape(Δky, (n_basis, nx, n_basis, nx))
          hamiltonian[2, i_orbital, :, :, 1, i_orbital, :, :] = reshape(Δky, (n_basis, nx, n_basis, nx))
        end
      end

      hamiltonian_matrix = reshape(hamiltonian,
                                   (n_nambu * n_orbital * n_basis * nx,
                                    n_nambu * n_orbital * n_basis * nx))
      eivals, eivecs = eig(Hermitian(hamiltonian_matrix))
      push!(all_eigenvalues_mixed, eivals)
      push!(all_eigenvectors_mixed, eivecs)
      push!(ky_bookkeep, ky * ones(Float64, size(eivals)))
    end

    all_eigenvalues_mixed = vcat(all_eigenvalues_mixed...)
    all_eigenvectors_mixed = hcat(all_eigenvectors_mixed...)
    ky_bookkeep = vcat(ky_bookkeep...)

    let 
      idxperm = sortperm(all_eigenvalues_mixed)
      all_eigenvalues_mixed = all_eigenvalues_mixed[idxperm]
      all_eigenvectors_mixed = all_eigenvectors_mixed[:, idxperm]
      ky_bookkeep = ky_bookkeep[idxperm]
    end

    result["mixed"] = Dict{String, Any}("eigenvalues" => all_eigenvalues_mixed,
                                        "eigenvectors" => all_eigenvectors_mixed,
                                        "momentums" => ky_bookkeep)
  end

  save(args["outfile"],
       "parameters", args,
       "mixed", result["mixed"],
       "momentum", result["momentum"])
end

main()

