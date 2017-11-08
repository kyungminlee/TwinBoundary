module OnebandModel

using Tightbinding

export OnebandParameter

type OnebandParameter
  t1 ::Float64
  t2 ::Float64
  μ  ::Float64
end

function hoppingList(param ::OnebandParameter)
  t1 = param.t1
  t2 = param.t2
  μ  = param.μ

  hops = [
    (( 0, 0), "s", "s", μ),
    (( 1, 0), "s", "s", t1), ((-1, 0), "s", "s", t1), 
    (( 0, 1), "s", "s", t1), (( 0,-1), "s", "s", t1),
    (( 1, 1), "s", "s", t2), ((-1, 1), "s", "s", t2), 
    (( 1,-1), "s", "s", t2), ((-1,-1), "s", "s", t2),
  ]
  return hops
end

function onebandModel(S ::Type, param ::OnebandParameter)
  @assert(S<: Number)
  tb_model = TightbindingModel{2, S}([("s", [0.0, 0.0])])

  hops = hoppingList(param)
  for (dis, ro, co, v) in hops
    add_hopping!(tb_model, dis, ro, co, v)
  end
  return tb_model
end


end