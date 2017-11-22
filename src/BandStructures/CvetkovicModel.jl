module CvetkovicModel

using Tightbinding
export CvetkovicParameter

type CvetkovicParameter
  ϵ_Γ ::Float64
  ϵ_1 ::Float64
  ϵ_2 ::Float64
  t_Γ ::Float64
  t_1 ::Float64
  t_2 ::Float64
  a_1 ::Float64
  a_3 ::Float64
  b   ::Float64
  c   ::Float64
  v   ::Float64
  p1  ::Float64
  p2  ::Float64
end

function hoppingList(param ::CvetkovicParameter)
  ϵ_Γ = param.ϵ_Γ
  t_Γ = param.t_Γ
  b   = param.b
  c   = param.c

  b_4 = 0.25 * b

  return [
    (( 0,  0), "A", "A", ϵ_Γ + 4 * t_Γ),
    (( 0,  0), "B", "B", ϵ_Γ + 4 * t_Γ),

    (( 1,  0), "A", "A", -t_Γ),    ((-1,  0), "A", "A", -t_Γ),
    (( 0,  1), "A", "A", -t_Γ),    (( 0, -1), "A", "A", -t_Γ),

    (( 1,  0), "B", "B", -t_Γ),    ((-1,  0), "B", "B", -t_Γ),
    (( 0,  1), "B", "B", -t_Γ),    (( 0, -1), "B", "B", -t_Γ),

    (( 1,  1), "A", "A", -b_4),    ((-1,  1), "A", "A",  b_4),
    (( 1, -1), "A", "A",  b_4),    ((-1, -1), "A", "A", -b_4),

    (( 1,  1), "B", "B",  b_4),    ((-1,  1), "B", "B", -b_4),
    (( 1, -1), "B", "B", -b_4),    ((-1, -1), "B", "B",  b_4),

    (( 1,  0), "A", "B", -c  ),    ((-1,  0), "A", "B", -c  ),
    (( 0,  1), "A", "B",  c  ),    (( 0, -1), "A", "B",  c  ),

    (( 1,  0), "B", "A", -c  ),    ((-1,  0), "B", "A", -c  ),
    (( 0,  1), "B", "A",  c  ),    (( 0, -1), "B", "A",  c  )
  ]
end


function spinfulHoppingList(param ::CvetkovicParameter;
  λ_Γ::Float64 = 40.0)
  ϵ_Γ = param.ϵ_Γ
  t_Γ = param.t_Γ
  b   = param.b
  c   = param.c

  b_4 = 0.25 * b

  h1 = [
    (( 0,  0), "A", "A", ϵ_Γ + 4 * t_Γ),
    (( 0,  0), "B", "B", ϵ_Γ + 4 * t_Γ),

    (( 1,  0), "A", "A", -t_Γ),    ((-1,  0), "A", "A", -t_Γ),
    (( 0,  1), "A", "A", -t_Γ),    (( 0, -1), "A", "A", -t_Γ),

    (( 1,  0), "B", "B", -t_Γ),    ((-1,  0), "B", "B", -t_Γ),
    (( 0,  1), "B", "B", -t_Γ),    (( 0, -1), "B", "B", -t_Γ),

    (( 1,  1), "A", "A", -b_4),    ((-1,  1), "A", "A",  b_4),
    (( 1, -1), "A", "A",  b_4),    ((-1, -1), "A", "A", -b_4),

    (( 1,  1), "B", "B",  b_4),    ((-1,  1), "B", "B", -b_4),
    (( 1, -1), "B", "B", -b_4),    ((-1, -1), "B", "B",  b_4),

    (( 1,  0), "A", "B", -c  ),    ((-1,  0), "A", "B", -c  ),
    (( 0,  1), "A", "B",  c  ),    (( 0, -1), "A", "B",  c  ),

    (( 1,  0), "B", "A", -c  ),    ((-1,  0), "B", "A", -c  ),
    (( 0,  1), "B", "A",  c  ),    (( 0, -1), "B", "A",  c  )
  ]

  # SOC First
  h2 = [
    ((0,0), "Aup", "Bup", -0.5im * λ_Γ),
    ((0,0), "Bup", "Aup",  0.5im * λ_Γ),
    ((0,0), "Adn", "Bdn",  0.5im * λ_Γ),
    ((0,0), "Bdn", "Adn", -0.5im * λ_Γ)
  ]
  for spin in ["up", "dn"], (dis, ro, co, val) in h1
    push!(h2, (dis, join([ro, spin]), join([co, spin]), val))
  end

  return h2
end


function cvetkovicModel(S::Type, param::CvetkovicParameter)
  @assert( S <: Number )
  tb_model = TightbindingModel{2, S}([
    ("A", [0.0, 0.0]),
    ("B", [0.0, 0.0])
  ])

  hops = hoppingList(param)
  for (dis, ro, co, v) in hops
    add_hopping!(tb_model, dis, ro, co, v)
  end
  return tb_model
end

function cvetkovicSpinfulModel(S::Type, param::CvetkovicParameter;
  λ_Γ::Float64 = 40.0)
  @assert( S <: Number )
  tb_model = TightbindingModel{2, S}([
    ("Aup", [0.0, 0.0]),
    ("Adn", [0.0, 0.0]),
    ("Bup", [0.0, 0.0]),
    ("Bdn", [0.0, 0.0])
  ])

  hops = spinfulHoppingList(param; λ_Γ=λ_Γ)
  for (dis, ro, co, v) in hops
    add_hopping!(tb_model, dis, ro, co, v)
  end
  return tb_model
end


end
