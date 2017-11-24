using Base.Test
using PyCall
using ProgressMeter
using PyPlot

include("spinful_daghofer.jl")

## Test Uniform
#=
let
  const λ = 0.1
  const Δd = 0.2
  const Δs = 0.3

  const n_nambu = 2
  const n_basis = 2
  const n_orbital = 3
  const n_spin = 2

  const UP = 1
  const DN = 2

  const nx = 3
  const ny = 4

  const kxs = linspace(0, 2π, nx+1)[1:end-1]
  const kys = linspace(0, 2π, ny+1)[1:end-1]

  const daghofer_parameter = DaghoferModel.DaghoferParameter(
          0.02,  0.06,  0.03, -0.01,  0.20,
          0.30, -0.20,  0.10,  0.40,  0.20)

  const spinful_daghofer_model = DaghoferModel.spinfulDaghoferModel(Complex128, daghofer_parameter, λ)
  const double_spinful_daghofer_model = Tightbinding.double_unitcell(spinful_daghofer_model)

  # 1/3 Momentum Space
  double_spinful_daghofer_momentumspace = Tightbinding.make_momentumspace(double_spinful_daghofer_model)
  eigenvalues_momentumspace = []
  eigenvectors_momentumspace = []
  for kx in kxs, ky in kys
    h = double_spinful_daghofer_momentumspace(kx, ky)
    eigenvalues, eigenvectors = eig(Hermitian(0.5 * (h + h')))
    push!(eigenvalues_momentumspace, eigenvalues)
    #push!(eigenvectors_momentumspace, eigenvectors)
  end
  eigenvalues_momentumspace = vcat(eigenvalues_momentumspace...)
  sort!(eigenvalues_momentumspace)

  # 2/3 Mixed Space (x, ky)
  double_spinful_daghofer_mixedspace = Tightbinding.make_mixedspace(double_spinful_daghofer_model, nx)
  eigenvalues_mixedspace = []
  eigenvectors_mixedspace = []
  for ky in kys
    h = double_spinful_daghofer_mixedspace(ky)
    eigenvalues, eigenvectors = eig(Hermitian(0.5 * (h + h')))
    push!(eigenvalues_mixedspace, eigenvalues)
    #push!(eigenvectors_mixedspace, eigenvectors)
  end
  eigenvalues_mixedspace = vcat(eigenvalues_mixedspace...)
  sort!(eigenvalues_mixedspace)

  # 1/3 Real Space
  eigenvalues_realspace, eigenvectors_realspace = let
    hamiltonian_realspace = Tightbinding.make_realspace_dense(double_spinful_daghofer_model, (nx, ny))
    eig(Hermitian(0.5 * (hamiltonian_realspace + hamiltonian_realspace')  ))
  end

  plot(eigenvalues_momentumspace .+ 0)
  plot(eigenvalues_mixedspace .+ 1)
  plot(eigenvalues_realspace .+ 2)
  println("BEFORE SHOW")
  show()
  println("AFTER SHOW")
end
=#




#=
Comparing realspace, mixedspace, and momentumspace pairings only.
formfactor is d-wave.
=#
@testset "Pairing d-wave" begin
  nx, ny = 4, 3

  const kxs = linspace(0, 2π, nx+1)[1:end-1]
  const kys = linspace(0, 2π, ny+1)[1:end-1]

  pairing_realspace = make_pairing_realspace_dense(ones(Complex128, (2, nx, ny)), FORMFACTORS["d1"]; periodic=true)
  pairing_mixedspace = make_pairing_mixedspace(ones(Complex128, (2, nx)), FORMFACTORS["d1"]; periodic=true)
  pairing_momentumspace = make_pairing_momentumspace(1.0, FORMFACTORS["d1"])
  eigenvalues_realspace = eigvals(pairing_realspace)
  eigenvalues_mixedspace = []
  eigenvalues_momentumspace = []
  for ky in kys
    foo = pairing_mixedspace(ky)
    push!(eigenvalues_mixedspace, eigvals(foo))
    for kx in kxs
      bar = pairing_momentumspace(kx, ky)
      push!(eigenvalues_momentumspace, eigvals(bar))
    end
  end
  eigenvalues_mixedspace = vcat(eigenvalues_mixedspace...)
  eigenvalues_momentumspace = vcat(eigenvalues_momentumspace...)

  sort!(eigenvalues_realspace)
  sort!(eigenvalues_mixedspace)
  sort!(eigenvalues_momentumspace)

  @test isapprox(eigenvalues_momentumspace, eigenvalues_mixedspace)
  @test isapprox(eigenvalues_momentumspace, eigenvalues_realspace)

  figure()
  plot(eigenvalues_momentumspace .+ 0, ".-", alpha=0.5, label="momentumspace", linewidth=2)
  plot(eigenvalues_mixedspace .+ 0, ".-", alpha=0.5, label="mixedspace", linewidth=2)
  plot(eigenvalues_realspace .+ 0, ".-", alpha=0.5, label="realspace", linewidth=2)
  title("Pairing Eigenvalues. d-wave pairing")
  legend()
end

#=
Comparing realspace, mixedspace, and momentumspace pairings only.
formfactor is s-wave.
=#
@testset "Pairing s-wave" begin
  nx, ny = 4, 3

  const kxs = linspace(0, 2π, nx+1)[1:end-1]
  const kys = linspace(0, 2π, ny+1)[1:end-1]

  pairing_realspace = make_pairing_realspace_dense(ones(Complex128, (2, nx, ny)), FORMFACTORS["s1"]; periodic=true)
  pairing_mixedspace = make_pairing_mixedspace(ones(Complex128, (2, nx)), FORMFACTORS["s1"]; periodic=true)
  pairing_momentumspace = make_pairing_momentumspace(1.0, FORMFACTORS["s1"])
  eigenvalues_realspace = eigvals(pairing_realspace)
  eigenvalues_mixedspace = []
  eigenvalues_momentumspace = []
  for ky in kys
    foo = pairing_mixedspace(ky)
    push!(eigenvalues_mixedspace, eigvals(foo))
    for kx in kxs
      bar = pairing_momentumspace(kx, ky)
      push!(eigenvalues_momentumspace, eigvals(bar))
    end
  end
  eigenvalues_mixedspace = vcat(eigenvalues_mixedspace...)
  eigenvalues_momentumspace = vcat(eigenvalues_momentumspace...)

  sort!(eigenvalues_realspace)
  sort!(eigenvalues_mixedspace)
  sort!(eigenvalues_momentumspace)

  @test isapprox(eigenvalues_momentumspace, eigenvalues_mixedspace)
  @test isapprox(eigenvalues_momentumspace, eigenvalues_realspace)

  figure()
  plot(eigenvalues_momentumspace .+ 0, ".-", alpha=0.5, label="momentumspace", linewidth=2)
  plot(eigenvalues_mixedspace .+ 0, ".-", alpha=0.5, label="mixedspace", linewidth=2)
  plot(eigenvalues_realspace .+ 0, ".-", alpha=0.5, label="realspace", linewidth=2)
  title("Pairing Eigenvalues. s-wave pairing")
  legend()
end






# TEST INCLUDING PAIRING

@testset "Full D" begin
  # no param confirmed
  const λ = 0.1 # confirmed
  const Δd = 0.2
  const Δs = 0.0
  const Δp = 0.0
  const ξ₀ = 1E-8
  const ξ₁ = 1E-8

  const n_nambu = 2
  const n_basis = 2
  const n_orbital = 3
  const n_spin = 2

  const UP = 1
  const DN = 2

  const nx = 3
  const ny = 4

  const kxs = linspace(0, 2π, nx+1)[1:end-1]
  const kys = linspace(0, 2π, ny+1)[1:end-1]

  const daghofer_parameter = DaghoferModel.DaghoferParameter(
          0.02,  0.06,  0.03, -0.01,  0.20,
          0.30, -0.20,  0.10,  0.40,  0.20)

  result_momentumspace = compute_momentumspace(nx, ny, λ, Δd, Δs, Δp, ξ₀, ξ₁)
  result_mixedspace = compute_mixedspace(nx, ny, λ, Δd, Δs, Δp, ξ₀, ξ₁)
  result_realspace = compute_realspace_dense(nx, ny, λ, Δd, Δs, Δp, ξ₀, ξ₁)

  eigenvalues_momentumspace = vcat(result_momentumspace["eigenvalues"]...)
  eigenvalues_mixedspace = vcat(result_mixedspace["eigenvalues"]...)
  eigenvalues_realspace = result_realspace["eigenvalues"]

  @test isapprox(eigenvalues_momentumspace, eigenvalues_mixedspace)
  @test isapprox(eigenvalues_momentumspace, eigenvalues_realspace)

  figure()
  plot(eigenvalues_momentumspace .+ 0, ".-", alpha=0.5, label="momentumspace", linewidth=2)
  plot(eigenvalues_mixedspace .+ 0, ".-", alpha=0.5, label="mixedspace", linewidth=2)
  plot(eigenvalues_realspace .+ 0, ".-", alpha=0.5, label="realspace", linewidth=2)
  title("Full Hamiltonian with uniform d-wave only.")
  legend()
end


@testset "Full D+S+P" begin
  # no param confirmed
  const λ = 0.0 # confirmed
  const Δd = 0.3
  const Δs = 0.2
  const Δp = 1.0
  const ξ₀ = 1.0
  const ξ₁ = 2.0

  const n_nambu = 2
  const n_basis = 2
  const n_orbital = 3
  const n_spin = 2

  const UP = 1
  const DN = 2

  const nx = 3
  const ny = 4

  const kxs = linspace(0, 2π, nx+1)[1:end-1]
  const kys = linspace(0, 2π, ny+1)[1:end-1]

  const daghofer_parameter = DaghoferModel.DaghoferParameter(
          0.02,  0.06,  0.03, -0.01,  0.20,
          0.30, -0.20,  0.10,  0.40,  0.20)

  result_momentumspace = compute_momentumspace(nx, ny, λ, Δd, Δs, Δp, ξ₀, ξ₁)
  result_mixedspace = compute_mixedspace(nx, ny, λ, Δd, Δs, Δp, ξ₀, ξ₁)
  result_realspace = compute_realspace_dense(nx, ny, λ, Δd, Δs, Δp, ξ₀, ξ₁)

  eigenvalues_momentumspace = vcat(result_momentumspace["eigenvalues"]...)
  eigenvalues_mixedspace = vcat(result_mixedspace["eigenvalues"]...)
  eigenvalues_realspace = result_realspace["eigenvalues"]

  # domain wall breaks translation symmetry
  @test ! isapprox(eigenvalues_momentumspace, eigenvalues_realspace)
  # but Ty is still preserved
  @test isapprox(eigenvalues_mixedspace, eigenvalues_realspace)

  figure()
  plot(eigenvalues_momentumspace .+ 0, ".-", alpha=0.5, label="momentumspace", linewidth=2)
  plot(eigenvalues_mixedspace .+ 0, ".-", alpha=0.5, label="mixedspace", linewidth=2)
  plot(eigenvalues_realspace .+ 0, ".-", alpha=0.5, label="realspace", linewidth=2)
  title("D+S+P")
  legend()
end

show()
