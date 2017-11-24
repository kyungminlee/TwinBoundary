using ProgressMeter

function compute_momentumspace(nx ::Int, ny ::Int, Δd ::Real, Δs ::Real, Δp ::Real, ξ₀ ::Real, ξ₁ ::Real)
  const n_nambu = 2
  const n_basis = 2
  const n_orbital = 1

  const kxs = linspace(0, 2π, nx+1)[1:end-1]
  const kys = linspace(0, 2π, ny+1)[1:end-1]

  const oneband_parameter = OnebandModel.OnebandParameter(
    -1.0, 0.0, -3.0)

  const oneband_model = OnebandModel.onebandModel(Complex128, oneband_parameter)
  const double_oneband_model = Tightbinding.double_unitcell(oneband_model)

  kxky_bookkeep = []
  all_eigenvalues_momentum = []
  all_eigenvectors_momentum = []

  double_oneband_momentumspace = Tightbinding.make_momentumspace(double_oneband_model)
  pairing_momentumspace_d = make_pairing_momentumspace(Δd, formfactors["d1"])
  pairing_momentumspace_s = make_pairing_momentumspace(Δs, formfactors["s1"])
  @showprogress for (ky, kx) in Base.product(kys, kxs)
    hk = zeros(Complex128, (n_nambu, n_orbital, n_basis, n_nambu, n_orbital, n_basis))

    let
      tk = double_oneband_momentumspace(kx, ky)
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

  return Dict{String, Any}("eigenvalues" => all_eigenvalues_momentum,
                           "eigenvectors" => all_eigenvectors_momentum,
                           "momentums" => kxky_bookkeep)
end


function compute_mixedspace(nx ::Int, ny ::Int, Δd ::Real, Δs ::Real, Δp ::Real, ξ₀ ::Real, ξ₁ ::Real)
  const n_nambu = 2
  const n_basis = 2
  const n_orbital = 1

  const kxs = linspace(0, 2π, nx+1)[1:end-1]
  const kys = linspace(0, 2π, ny+1)[1:end-1]

  const oneband_parameter = OnebandModel.OnebandParameter(
    -1.0, 0.0, -3.0)

  const oneband_model = OnebandModel.onebandModel(Complex128, oneband_parameter)
  const double_oneband_model = Tightbinding.double_unitcell(oneband_model)

  ky_bookkeep = []
  all_eigenvalues_mixed = []
  all_eigenvectors_mixed = []

  double_oneband_mixedspace = Tightbinding.make_mixedspace(double_oneband_model, nx; periodic=false)

  pairing_orderparameter = [(        Δs * tanh((x - 1 - nx)/(sqrt(0.5) * ξ₀))
                             + 1im * Δp * tanh((x - 1 - nx)/(sqrt(0.5) * ξ₁))) for x in 1:2*nx]
  pairing_orderparameter = reshape(pairing_orderparameter, (2, nx))

  pairing_mixedspace_uniform = make_pairing_mixedspace(Δd * ones(Complex128, (2, nx)),
                                                       formfactors["d1"];
                                                       periodic=false)

  pairing_mixedspace = make_pairing_mixedspace(pairing_orderparameter,
                                               formfactors["s1"];
                                               periodic=false)

  @showprogress for ky in kys
    hamiltonian = zeros(Complex128, (n_nambu, n_orbital, n_basis, nx,
                                     n_nambu, n_orbital, n_basis, nx))
    let
      tky  = double_oneband_mixedspace( ky)
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

  return Dict{String, Any}("eigenvalues" => all_eigenvalues_mixed,
                           "eigenvectors" => all_eigenvectors_mixed,
                           "momentums" => ky_bookkeep)
end

