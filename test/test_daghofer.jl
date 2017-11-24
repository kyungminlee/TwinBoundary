using Base.Test

include("../src/Tightbinding.jl")
include("../src/BandStructures/DaghoferModel.jl")

include("../src/DoublePlaquetteModel/PairingModel.jl")
include("../src/DoublePlaquetteModel/NaivePairingModel.jl")

daghofer_parameter = DaghoferModel.DaghoferParameter(
    0.02,
    0.06,
    0.03,
    -0.01,
    0.2,
    0.3,
    -0.2,
    0.1,
    0.4,
    0.2)

zero_daghofer_parameter = DaghoferModel.DaghoferParameter(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

const n_nambu = 2
const n_basis = 2
const n_orbital = 3

const nx, ny = 3, 5
const kxs = linspace(0, 2π, nx+1)[1:end-1]
const kys = linspace(0, 2π, ny+1)[1:end-1]

const daghofer_model = DaghoferModel.daghoferModel(Complex128, daghofer_parameter)
const double_daghofer_model = Tightbinding.double_unitcell(daghofer_model)

const Δ₀ = 0.5 + 0.2im
const Δ1d = Δ₀ * ones(Complex128, (2, nx))
const Δ2d = Δ₀ * ones(Complex128, (2, nx, ny))

all_eigenvalues_real = Float64[]

@testset "2D" begin
  double_daghofer_realspace_dense = Tightbinding.make_realspace_dense(double_daghofer_model, (nx, ny); periodic=true)

  hamiltonian = zeros(Complex128, (n_nambu, n_orbital, n_basis, nx, ny, n_nambu, n_orbital, n_basis, nx, ny))

  let t = double_daghofer_realspace_dense
    hamiltonian[1, :, :, :, :, 1, :, :, :, :] = reshape( t, (n_orbital, n_basis, nx, ny, n_orbital, n_basis, nx, ny))
    hamiltonian[2, :, :, :, :, 2, :, :, :, :] = reshape(-t, (n_orbital, n_basis, nx, ny, n_orbital, n_basis, nx, ny))
  end

  Δx = make_pairing_realspace_dense(Δ2d, FORMFACTORS["s1"]; periodic=true)
  # test if Δx is Hermitian
  @test isapprox(Δx, Δx'; atol=1E-8)

  let
    Δx_naive = NaivePairingModel.make_pairing_realspace_dense(Δ₀, nx, ny; periodic=true)
    # test if it matches with naive version
    @test isapprox(Δx, Δx_naive; atol=1E-8)
  end

  for i_orbital in 1:n_orbital
    hamiltonian[1, i_orbital, :, :, :, 2, i_orbital, :, :, :] = reshape(Δx, (n_basis, nx, ny, n_basis, nx, ny))
    hamiltonian[2, i_orbital, :, :, :, 1, i_orbital, :, :, :] = reshape(Δx, (n_basis, nx, ny, n_basis, nx, ny))
  end

  hamiltonian_matrix = reshape(hamiltonian,
                              (n_nambu * n_orbital * n_basis * nx * ny,
                                n_nambu * n_orbital * n_basis * nx * ny))

  eivals = eigvals(Hermitian(hamiltonian_matrix))
  append!(all_eigenvalues_real, eivals)
end

all_eigenvalues_mixed = Float64[]
@testset "1D" begin
  double_daghofer_mixedspace = Tightbinding.make_mixedspace(double_daghofer_model, nx; periodic=true)

  fullsize = (n_orbital, n_basis, nx, n_orbital, n_basis, nx)
  matsize = (n_orbital * n_basis * nx, n_orbital * n_basis * nx)

  pairing_mixedspace = make_pairing_mixedspace(Δ1d, FORMFACTORS["s1"]; periodic=true)
  pairing_mixedspace_naive = NaivePairingModel.make_pairing_mixedspace(Δ₀, nx)

  for ky in kys
    hamiltonian = zeros(Complex128, (n_nambu, n_orbital, n_basis, nx, n_nambu, n_orbital, n_basis, nx))

    Δky  = pairing_mixedspace(ky)
    @test isapprox(Δky, Δky'; atol=1E-8)

    let
      Δky_naive = pairing_mixedspace_naive(ky)
      @test isapprox(Δky, Δky_naive; atol=1E-8)
    end
    let
      tky  = double_daghofer_mixedspace( ky)
      let
        tmky = double_daghofer_mixedspace(-ky)
        @test isapprox(tky, tmky.'; atol=1E-8)
      end
      hamiltonian[1, :, :, :, 1, :, :, :] = reshape( tky, fullsize)
      hamiltonian[2, :, :, :, 2, :, :, :] = reshape(-tky, fullsize)
    end

    for i_orbital in 1:n_orbital
      hamiltonian[1, i_orbital, :, :, 2, i_orbital, :, :] = reshape(Δky, (n_basis, nx, n_basis, nx))
      hamiltonian[2, i_orbital, :, :, 1, i_orbital, :, :] = reshape(Δky, (n_basis, nx, n_basis, nx))
    end

    hamiltonian_matrix = reshape(hamiltonian, (n_nambu * n_orbital * n_basis * nx, n_nambu * n_orbital * n_basis * nx))

    eivals = eigvals(Hermitian(hamiltonian_matrix))
    append!(all_eigenvalues_mixed, eivals)
  end
end

all_eigenvalues_momentum = Float64[]
@testset "0D" begin
  double_daghofer_momentumspace = Tightbinding.make_momentumspace(double_daghofer_model)

  pairing_momentumspace = make_pairing_momentumspace(Δ₀, FORMFACTORS["s1"])
  pairing_momentumspace_naive = NaivePairingModel.make_pairing_momentumspace(Δ₀)

  for kx in kxs, ky in kys
    hk = zeros(Complex128, (n_nambu, n_orbital, n_basis, n_nambu, n_orbital, n_basis))

    let
      tk = double_daghofer_momentumspace(kx, ky)
      let
        tmk = double_daghofer_momentumspace(-kx, -ky)
        @test isapprox(tk, tmk.'; atol=1E-8)
      end
      hk[1, :, :, 1, :, :] = reshape( tk, (n_orbital, n_basis, n_orbital, n_basis))
      hk[2, :, :, 2, :, :] = reshape(-tk, (n_orbital, n_basis, n_orbital, n_basis))
    end

    let
      Δk = pairing_momentumspace(kx, ky)
      #@test isapprox(Δk, Δk')
      let
        Δk_naive = pairing_momentumspace_naive(kx, ky)
        @test isapprox(Δk, Δk_naive; atol=1E-8)
      end

      for i_orbital in 1:n_orbital
        hk[1, i_orbital, :, 2, i_orbital, :] = reshape(Δk, (n_basis, n_basis))
        hk[2, i_orbital, :, 1, i_orbital, :] = reshape(Δk, (n_basis, n_basis))
      end
    end

    hamiltonian_matrix = reshape(hk, (n_nambu * n_orbital * n_basis, n_nambu * n_orbital * n_basis))
    @test isapprox(hamiltonian_matrix, hamiltonian_matrix'; atol=1E-8)

    hamiltonian_matrix = 0.5 * (hamiltonian_matrix + hamiltonian_matrix')
    eivals = eigvals(Hermitian(hamiltonian_matrix))
    append!(all_eigenvalues_momentum, eivals)
  end
end

@testset "eigenvalues" begin
  sort!(all_eigenvalues_real)
  sort!(all_eigenvalues_mixed)
  sort!(all_eigenvalues_momentum)

  @test isapprox(all_eigenvalues_real, all_eigenvalues_mixed)
  @test isapprox(all_eigenvalues_real, all_eigenvalues_momentum)
end
