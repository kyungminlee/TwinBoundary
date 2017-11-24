#
# Compute Real Space Spectrum With an Edge
#
#
using Base.Test
using PyPlot

include("spinful_daghofer.jl")

@testset "nambufy" begin
  # no param confirmed
  const λ = 0.1 # confirmed
  const Δd = 0.0
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
  const spinful_daghofer_model = DaghoferModel.spinfulDaghoferModel(Complex128, daghofer_parameter, λ)
  const double_spinful_daghofer_model = Tightbinding.double_unitcell(spinful_daghofer_model)
  const nambu_double_spinful_daghofer_model = Tightbinding.nambufy(double_spinful_daghofer_model)
  nambu_dense = Tightbinding.make_realspace_dense(nambu_double_spinful_daghofer_model, (nx, ny))
  nambu_sparse = Tightbinding.make_realspace_sparse(nambu_double_spinful_daghofer_model, (nx, ny))


  direct_dense = let
    hamiltonian = zeros(Complex128, (n_nambu, n_orbital* n_spin* n_basis* nx* ny,
                                     n_nambu, n_orbital* n_spin* n_basis* nx* ny))
    hopping = Tightbinding.make_realspace_dense(double_spinful_daghofer_model, (nx, ny))
    hamiltonian[1, :, 1, :] =  hopping
    hamiltonian[2, :, 2, :] = -transpose(hopping)
    reshape(hamiltonian, (n_nambu * n_orbital* n_spin* n_basis* nx* ny,
                          n_nambu * n_orbital* n_spin* n_basis* nx* ny))
  end

  @test all(isapprox.(direct_dense, nambu_dense))
  @test all(isapprox.(direct_dense, full(nambu_sparse)))

  eigenvalues_dense = eigvals(Hermitian(nambu_dense))
  #eigenvalues_direct_dense = eigvals(Hermitian(direct_dense))
  eigenvalues_sparse, eigenvectors_sparse = eigs(Hermitian(nambu_sparse); nev=10, which=:SM)

  eigenvalues_sparse = real(eigenvalues_sparse)

  e1d = sort(eigenvalues_dense[eigenvalues_dense .>= 0.0])
  e2d = -sort(-eigenvalues_dense[eigenvalues_dense .< 0.0])

  e1s = sort(eigenvalues_sparse[eigenvalues_sparse .>= 0.0])
  e2s = -sort(-eigenvalues_sparse[eigenvalues_sparse .< 0.0])

  figure()
  plot(e1d, ".-", alpha=0.5, label="Dense", linewidth=2)
  plot(e1s, ".-", alpha=0.5, label="Sparse", linewidth=2)
  plot(e2d, ".-", alpha=0.5, label="Dense", linewidth=2)
  plot(e2s, ".-", alpha=0.5, label="Sparse", linewidth=2)

  xlim(-0.5, max(length(e1s), length(e2s)) + 1.5)
  ylim(e2s[end]-0.1, e1s[end]+0.1)
  title("Normal State in Nambu Space (Purely)")
  legend()
end

@testset "NormalNambu" begin
  # no param confirmed
  const λ = 0.1 # confirmed
  const Δd = 0.0
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

  result_realspace_dense = compute_realspace_dense(nx, ny, λ, Δd, Δs, Δp, ξ₀, ξ₁, periodic=false)
  result_realspace_sparse = compute_realspace_sparse(nx, ny, λ, Δd, Δs, Δp, ξ₀, ξ₁, periodic=false)

  eigenvalues_realspace_dense = result_realspace_dense["eigenvalues"]
  eigenvalues_realspace_sparse = result_realspace_sparse["eigenvalues"]

  e1d = sort( eigenvalues_realspace_dense[ eigenvalues_realspace_dense .>= 0.0] )
  e1s = sort( eigenvalues_realspace_sparse[ eigenvalues_realspace_sparse .>= 0.0] )

  e2d = -sort(-eigenvalues_realspace_dense[ eigenvalues_realspace_dense .< 0.0] )
  e2s = -sort(-eigenvalues_realspace_sparse[ eigenvalues_realspace_sparse .< 0.0] )

  figure()
  plot(e1d, ".-", alpha=0.5, label="Dense", linewidth=2)
  plot(e1s, ".-", alpha=0.5, label="Sparse", linewidth=2)
  plot(e2d, ".-", alpha=0.5, label="Dense", linewidth=2)
  plot(e2s, ".-", alpha=0.5, label="Sparse", linewidth=2)

  xlim(-0.5, max(length(e1s), length(e2s)) + 1.5)
  ylim(e2s[end]-0.1, e1s[end]+0.1)
  title("Normal State in Nambu Space")
  legend()
end

@testset "D-wave" begin
  # no param confirmed
  const λ = 0.1 # confirmed
  const Δd = 0.05
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

  result_realspace_dense = compute_realspace_dense(nx, ny, λ, Δd, Δs, Δp, ξ₀, ξ₁, periodic=false)
  result_realspace_sparse = compute_realspace_sparse(nx, ny, λ, Δd, Δs, Δp, ξ₀, ξ₁, periodic=false)

  eigenvalues_realspace_dense = result_realspace_dense["eigenvalues"]
  eigenvalues_realspace_sparse = result_realspace_sparse["eigenvalues"]

  e1d = sort( eigenvalues_realspace_dense[ eigenvalues_realspace_dense .>= 0.0] )
  e1s = sort( eigenvalues_realspace_sparse[ eigenvalues_realspace_sparse .>= 0.0] )

  e2d = -sort(-eigenvalues_realspace_dense[ eigenvalues_realspace_dense .< 0.0] )
  e2s = -sort(-eigenvalues_realspace_sparse[ eigenvalues_realspace_sparse .< 0.0] )

  figure()
  plot(e1d, ".-", alpha=0.5, label="Dense", linewidth=2)
  plot(e1s, ".-", alpha=0.5, label="Sparse", linewidth=2)
  plot(e2d, ".-", alpha=0.5, label="Dense", linewidth=2)
  plot(e2s, ".-", alpha=0.5, label="Sparse", linewidth=2)

  xlim(-0.5, max(length(e1s), length(e2s)) + 1.5)
  ylim(e2s[end]-0.1, e1s[end]+0.1)
  title("D-wave Only")
  legend()
end


@testset "Full" begin
  # no param confirmed
  const λ = 0.1 # confirmed
  const Δd = 0.05
  const Δs = 0.025
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

  result_realspace_dense = compute_realspace_dense(nx, ny, λ, Δd, Δs, Δp, ξ₀, ξ₁, periodic=false)
  result_realspace_sparse = compute_realspace_sparse(nx, ny, λ, Δd, Δs, Δp, ξ₀, ξ₁, periodic=false)

  eigenvalues_realspace_dense = result_realspace_dense["eigenvalues"]
  eigenvalues_realspace_sparse = result_realspace_sparse["eigenvalues"]

  e1d = sort( eigenvalues_realspace_dense[ eigenvalues_realspace_dense .>= 0.0] )
  e1s = sort( eigenvalues_realspace_sparse[ eigenvalues_realspace_sparse .>= 0.0] )

  e2d = -sort(-eigenvalues_realspace_dense[ eigenvalues_realspace_dense .< 0.0] )
  e2s = -sort(-eigenvalues_realspace_sparse[ eigenvalues_realspace_sparse .< 0.0] )

  figure()
  plot(e1d, ".-", alpha=0.5, label="Dense", linewidth=2)
  plot(e1s, ".-", alpha=0.5, label="Sparse", linewidth=2)
  plot(e2d, ".-", alpha=0.5, label="Dense", linewidth=2)
  plot(e2s, ".-", alpha=0.5, label="Sparse", linewidth=2)

  xlim(-0.5, max(length(e1s), length(e2s)) + 1.5)
  ylim(e2s[end]-0.1, e1s[end]+0.1)
  title("Full")
  legend()
end

show()

# Tests Successfull 2017-11-24
