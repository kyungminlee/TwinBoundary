using Base.Test



include("src/Tightbinding.jl")
include("src/DaghoferModel.jl")

include("src/PairingModel.jl")
include("src/NaivePairingModel.jl")

include("spinful_daghofer.jl")


@testset "Test Nambufy" begin
  const n_nambu = 2
  const n_basis = 2
  const n_orbital = 3
  const n_spin = 2

  const UP = 1
  const DN = 2
  const λ = 0.1

  const daghofer_parameter = DaghoferModel.DaghoferParameter(
            0.02,  0.06,  0.03, -0.01,  0.20,
            0.30, -0.20,  0.10,  0.40,  0.20)
  const spinful_daghofer_model = DaghoferModel.spinfulDaghoferModel(Complex128, daghofer_parameter, λ)
  const double_spinful_daghofer_model = Tightbinding.double_unitcell(spinful_daghofer_model)
  const nambu_double_spinful_daghofer_model = Tightbinding.nambufy(double_spinful_daghofer_model)

  const nx = 4
  const ny = 3

  const kxs = linspace(0, 2π, nx+1)[1:end-1]
  const kys = linspace(0, 2π, ny+1)[1:end-1]

  @testset "MomentumSpace" begin
    hopping_direct = Tightbinding.make_momentumspace(double_spinful_daghofer_model)
    hopping_nambufy = Tightbinding.make_momentumspace(nambu_double_spinful_daghofer_model)
    for kx in kxs, ky in kys
      hamiltonian_direct = zeros(Complex128, (n_nambu, n_orbital * n_spin * n_basis,
                                              n_nambu, n_orbital * n_spin * n_basis))
      hamiltonian_direct[1, :, 1, :] = hopping_direct(kx, ky)
      hamiltonian_direct[2, :, 2, :] = -transpose(hopping_direct(-kx, -ky))
      hamiltonian_direct = reshape(hamiltonian_direct, (n_nambu * n_orbital * n_spin * n_basis,
                                                        n_nambu * n_orbital * n_spin * n_basis))
      hamiltonian_nambufy = hopping_nambufy(kx, ky)
      @test isapprox( maximum( abs.(hamiltonian_direct - hamiltonian_nambufy) ), 0.0)
    end
  end

  @testset "MixedSpace" begin
    hopping_direct = Tightbinding.make_mixedspace(double_spinful_daghofer_model, nx)
    hopping_nambufy = Tightbinding.make_mixedspace(nambu_double_spinful_daghofer_model, nx)
    for ky in kys
      hamiltonian_direct = zeros(Complex128, (n_nambu, n_orbital * n_spin * n_basis * nx,
                                              n_nambu, n_orbital * n_spin * n_basis * nx))
      hamiltonian_direct[1, :, 1, :] = hopping_direct(ky)
      hamiltonian_direct[2, :, 2, :] = -transpose(hopping_direct(-ky))
      hamiltonian_direct = reshape(hamiltonian_direct, (n_nambu * n_orbital * n_spin * n_basis * nx,
                                                        n_nambu * n_orbital * n_spin * n_basis * nx))
      hamiltonian_nambufy = hopping_nambufy(ky)
      @test isapprox( maximum( abs.(hamiltonian_direct - hamiltonian_nambufy) ), 0.0)
    end
  end

  @testset "Realspace" begin
    hamiltonian_direct = zeros(Complex128, (n_nambu, n_orbital* n_spin* n_basis* nx* ny,
                                            n_nambu, n_orbital* n_spin* n_basis* nx* ny))
    let
      hopping = Tightbinding.make_realspace_dense(double_spinful_daghofer_model, (nx, ny); periodic=true)
      hamiltonian_direct[1, :, 1, :] =  hopping
      hamiltonian_direct[2, :, 2, :] = -transpose(hopping)
    end
    hamiltonian_direct = reshape(hamiltonian_direct,
                                 (n_nambu * n_orbital* n_spin* n_basis* nx* ny,
                                  n_nambu * n_orbital* n_spin* n_basis* nx* ny))
    hamiltonian_nambufy = Tightbinding.make_realspace_dense(nambu_double_spinful_daghofer_model, (nx, ny); periodic=true)
    @test isapprox( maximum( abs.(hamiltonian_direct - hamiltonian_nambufy) ), 0.0)
    @test ! isapprox(maximum( abs.(hamiltonian_direct)), 0.0)
  end

end
