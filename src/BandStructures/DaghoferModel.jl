module DaghoferModel

using Formatting
using Tightbinding

export DaghoferParameter

type DaghoferParameter
  t1  ::Float64
  t2  ::Float64
  t3  ::Float64
  t4  ::Float64
  t5  ::Float64
  t6  ::Float64
  t7  ::Float64
  t8  ::Float64
  Δxy ::Float64
  μ   ::Float64
end

function hoppingList(param ::DaghoferParameter)
    t1 = param.t1
    t2 = param.t2
    t3 = param.t3
    t4 = param.t4
    t5 = param.t5
    t6 = param.t6
    t7 = param.t7
    t8 = param.t8
    Δxy = param.Δxy
    μ   = param.μ
    hl = [
        (( 1, 0), "xz", "xz", t1),    ((-1, 0), "xz", "xz", t1),
        (( 0, 1), "xz", "xz", t2),    (( 0,-1), "xz", "xz", t2),
        (( 1, 1), "xz", "xz", t3),    (( 1,-1), "xz", "xz", t3),
        ((-1, 1), "xz", "xz", t3),    ((-1,-1), "xz", "xz", t3),
        (( 0, 0), "xz", "xz", -μ),

        (( 1, 0), "yz", "yz", t2),    ((-1, 0), "yz", "yz", t2),
        (( 0, 1), "yz", "yz", t1),    (( 0,-1), "yz", "yz", t1),
        (( 1, 1), "yz", "yz", t3),    (( 1,-1), "yz", "yz", t3),
        ((-1, 1), "yz", "yz", t3),    ((-1,-1), "yz", "yz", t3),
        (( 0, 0), "yz", "yz", -μ),

        (( 1, 0), "xy", "xy", t5),    ((-1, 0), "xy", "xy", t5),
        (( 0, 1), "xy", "xy", t5),    (( 0,-1), "xy", "xy", t5),
        (( 1, 1), "xy", "xy", t6),    (( 1,-1), "xy", "xy", t6),
        ((-1, 1), "xy", "xy", t6),    ((-1,-1), "xy", "xy", t6),
        (( 0, 0), "xy", "xy", Δxy-μ),

        (( 1, 1), "xz", "yz",-t4),    (( 1,-1), "xz", "yz", t4),
        ((-1, 1), "xz", "yz", t4),    ((-1,-1), "xz", "yz",-t4),

        (( 1, 0), "xz", "xy",-t7),    ((-1, 0), "xz", "xy", t7),

        (( 0, 1), "yz", "xy",-t7),    (( 0,-1), "yz", "xy", t7),

        (( 1, 1), "xz", "xy",-t8),    ((-1, 1), "xz", "xy", t8),
        (( 1,-1), "xz", "xy",-t8),    ((-1,-1), "xz", "xy", t8),

        (( 1, 1), "yz", "xy",-t8),    ((-1, 1), "yz", "xy",-t8),
        (( 1,-1), "yz", "xy", t8),    ((-1,-1), "yz", "xy", t8),
    ]
    hl2 = copy(hl)
    for ((dx, dy), ro, co, v) in hl
        if ro != co
            push!(hl2, ((-dx, -dy), co, ro, conj(v)))
        end
    end

    hl2
end

function daghoferModel(S :: Type, param :: DaghoferParameter)
    @assert( S <: Number )
    tb_model = TightbindingModel{2, S}([
        ("xz", [0.0, 0.0]),
        ("yz", [0.0, 0.0]) ,
        ("xy", [0.0, 0.0])])

    hops = hoppingList(param)
    for (dis, ro, co, v) in hops
        add_hopping!(tb_model, dis, ro, co, v)
    end
    return tb_model
end

function spinfulDaghoferModel(S ::Type,
                              param ::DaghoferParameter,
                              soc::Float64)
    @assert( S <: Number )
    tb_model = TightbindingModel{2, S}([
        ("xz_up", [0.0, 0.0]),
        ("yz_up", [0.0, 0.0]),
        ("xy_up", [0.0, 0.0]),
        ("xz_dn", [0.0, 0.0]),
        ("yz_dn", [0.0, 0.0]),
        ("xy_dn", [0.0, 0.0]),
    ])

    # spin-independent hoppings
    hops = hoppingList(param)
    for sp in ["up", "dn"], (dis, ro, co, v) in hops
        ros = format("{}_{}", ro, sp)
        cos = format("{}_{}", co, sp)
        add_hopping!(tb_model, dis, ros, cos, v)
    end
    # spin-orbit coupling
    add_hopping!(tb_model, (0,0), "xz_up", "yz_up", +im * soc)
    add_hopping!(tb_model, (0,0), "yz_up", "xz_up", -im * soc)
    add_hopping!(tb_model, (0,0), "xz_dn", "yz_dn", -im * soc)
    add_hopping!(tb_model, (0,0), "yz_dn", "xz_dn", +im * soc)

    return tb_model
end


end
