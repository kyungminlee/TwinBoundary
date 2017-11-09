using ProgressMeter
using PyCall
@pyimport numpy.linalg as npl

#function compute_momentumspace(param ::Dict)
function compute_momentumspace(nx ::Int, ny ::Int,
  λ ::Real,
  Δd ::Real, Δs ::Real, Δp ::Real, ξ₀ ::Real, ξ₁ ::Real)
  const n_nambu = 2
  const n_basis = 2
  const n_orbital = 3
  const n_spin = 2

  const UP = 1
  const DN = 2

  const kxs = linspace(0, 2π, nx+1)[1:end-1]
  const kys = linspace(0, 2π, ny+1)[1:end-1]

  const daghofer_parameter = DaghoferModel.DaghoferParameter(
          0.02,  0.06,  0.03, -0.01,  0.20,
          0.30, -0.20,  0.10,  0.40,  0.20)

  const spinful_daghofer_model = DaghoferModel.spinfulDaghoferModel(Complex128, daghofer_parameter, λ)
  const double_spinful_daghofer_model = Tightbinding.double_unitcell(spinful_daghofer_model)

  kxky_bookkeep = []
  all_eigenvalues_momentum = []
  all_eigenvectors_momentum = []

  double_spinful_daghofer_momentumspace = Tightbinding.make_momentumspace(double_spinful_daghofer_model)
  pairing_momentumspace_d = make_pairing_momentumspace(Δd, formfactors["d1"])
  pairing_momentumspace_s = make_pairing_momentumspace(Δs, formfactors["s1"])
  @showprogress for (ky, kx) in Base.product(kys, kxs)
    hk = zeros(Complex128, (n_nambu, n_orbital, n_spin, n_basis,
                            n_nambu, n_orbital, n_spin, n_basis))

    let
      tk = double_spinful_daghofer_momentumspace(kx, ky)
      Δk = pairing_momentumspace_d(kx, ky) + pairing_momentumspace_s(kx, ky)
      hk[1, :, :, :, 1, :, :, :] = reshape( tk, (n_orbital, n_spin, n_basis,
                                           n_orbital, n_spin, n_basis))
      hk[2, :, :, :, 2, :, :, :] = reshape(-tk, (n_orbital, n_spin, n_basis,
                                           n_orbital, n_spin, n_basis))

      # Assuming singlet (because I am using Delta k instead of Delta -k)
      for i_orbital in 1:n_orbital
        hk[1, i_orbital, UP, :, 2, i_orbital, DN, :] = reshape(Δk, (n_basis, n_basis))
        hk[1, i_orbital, DN, :, 2, i_orbital, UP, :] = -reshape(Δk, (n_basis, n_basis))
        hk[2, i_orbital, DN, :, 1, i_orbital, UP, :] = -reshape(Δk, (n_basis, n_basis))
        hk[2, i_orbital, UP, :, 1, i_orbital, DN, :] = reshape(Δk, (n_basis, n_basis))
      end
    end

    hamiltonian_matrix = reshape(hk, (n_nambu * n_orbital * n_spin * n_basis,
                                      n_nambu * n_orbital * n_spin * n_basis))
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

  return Dict{String, Any}("eigenvalues" => all_eigenvalues_momentum,
                           "eigenvectors" => all_eigenvectors_momentum,
                           "momentums" => kxky_bookkeep)
end


#function compute_mixedspace(param ::Dict)
function compute_mixedspace(nx ::Int, ny ::Int,
  λ ::Real,
  Δd ::Real, Δs ::Real, Δp ::Real, ξ₀ ::Real, ξ₁ ::Real)

  const n_nambu = 2
  const n_basis = 2
  const n_orbital = 3
  const n_spin = 2

  const UP = 1
  const DN = 2

  const kxs = linspace(0, 2π, nx+1)[1:end-1]
  const kys = linspace(0, 2π, ny+1)[1:end-1]

  const daghofer_parameter = DaghoferModel.DaghoferParameter(
          0.02,  0.06,  0.03, -0.01,  0.20,
          0.30, -0.20,  0.10,  0.40,  0.20)

  const spinful_daghofer_model = DaghoferModel.spinfulDaghoferModel(Complex128, daghofer_parameter, λ)
  const double_spinful_daghofer_model = Tightbinding.double_unitcell(spinful_daghofer_model)

  ky_bookkeep = []
  all_eigenvalues_mixed = []
  all_eigenvectors_mixed = []

  double_spinful_daghofer_mixedspace = Tightbinding.make_mixedspace(double_spinful_daghofer_model, nx; periodic=false)

  ψs = reshape( [Δs * tanh((x - 0.5 - nx)/(sqrt(2.0) * ξ₀)) for x in 1:2*nx], (2, nx) )
  ψd = reshape( [Δd for x in 1:2*nx], (2,nx))
  dz = reshape( [1im * Δp * tanh((x - 0.5 - nx)/(sqrt(2.0) * ξ₁)) for x in 1:2*nx], (2, nx) )

  pairing_mixedspace_uniform = make_pairing_mixedspace(ψd,
                                                       FORMFACTORS["d1"];
                                                       periodic=false)
  #
  # Δ = | 0  ψ| + |-dx+idy   dz  | = (ψ + ⃗d ⋅ ⃗σ) (i σ₂)
  #     |-ψ  0|   |  dz    dx+idy|
  #
  # Produces Δ(↑↓)
  pairing_mixedspace = make_pairing_mixedspace(ψs + dz,
                                               FORMFACTORS["s1"];
                                               periodic=false)
  @showprogress for ky in kys
    hamiltonian = zeros(Complex128, (n_nambu, n_orbital, n_spin, n_basis, nx,
                                     n_nambu, n_orbital, n_spin, n_basis, nx))
    let
      tky  = double_spinful_daghofer_mixedspace( ky)
      tmky = double_spinful_daghofer_mixedspace(-ky)
      Δky  = pairing_mixedspace_uniform( ky) + pairing_mixedspace( ky)
      Δmky = pairing_mixedspace_uniform(-ky) + pairing_mixedspace(-ky)

      hamiltonian[1, :, :, :, :, 1, :, :, :, :] = reshape( tky, (n_orbital, n_spin, n_basis, nx,
                                                                 n_orbital, n_spin, n_basis, nx))
      hamiltonian[2, :, :, :, :, 2, :, :, :, :] = reshape(-transpose(tmky),
                                                                (n_orbital, n_spin, n_basis, nx,
                                                                 n_orbital, n_spin, n_basis, nx))
      for i_orbital in 1:n_orbital
        hamiltonian[1, i_orbital, UP, :, :, 2, i_orbital, DN, :, :] =  reshape(Δky, (n_basis, nx, n_basis, nx))
        hamiltonian[1, i_orbital, DN, :, :, 2, i_orbital, UP, :, :] = -reshape(transpose(Δmky), (n_basis, nx, n_basis, nx))

        hamiltonian[2, i_orbital, UP, :, :, 1, i_orbital, DN, :, :] = -reshape(conj(Δmky), (n_basis, nx, n_basis, nx))
        hamiltonian[2, i_orbital, DN, :, :, 1, i_orbital, UP, :, :] = reshape(conj(transpose(Δky)), (n_basis, nx, n_basis, nx))
      end
    end
    #@show size(hamiltonian)
    #@show n_nambu, n_orbital, n_basis, nx

    hamiltonian_matrix = reshape(hamiltonian,
                                 (n_nambu * n_orbital * n_spin * n_basis * nx,
                                  n_nambu * n_orbital * n_spin * n_basis * nx))
    #eivals, eivecs = eig(Hermitian(hamiltonian_matrix))
    eivals, eivecs = npl.eigh(hamiltonian_matrix)
    push!(all_eigenvalues_mixed, eivals)
    push!(all_eigenvectors_mixed, eivecs)
    push!(ky_bookkeep, ky * ones(Float64, size(eivals)))
  end

  all_eigenvalues_mixed  = vcat(all_eigenvalues_mixed...)
  all_eigenvectors_mixed = hcat(all_eigenvectors_mixed...)
  ky_bookkeep = vcat(ky_bookkeep...)

  let
    idxperm = sortperm(all_eigenvalues_mixed)
    all_eigenvalues_mixed = all_eigenvalues_mixed[idxperm]
    all_eigenvectors_mixed = all_eigenvectors_mixed[:, idxperm]
    ky_bookkeep = ky_bookkeep[idxperm]
  end

  return Dict{String, Any}("eigenvalues" => all_eigenvalues_mixed,
                           "eigenvectors" => all_eigenvectors_mixed,
                           "momentums" => ky_bookkeep)
end
