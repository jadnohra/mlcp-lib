module Lu
#
  importall Shared_types
  using Shared_funcs
  using Conv
#
  :Dense
  :Sparse
#
  function MatrixTypesT(T) return Union(Matrix{T}, SparseMatrixCSC{T}) end
  MatrixTypes = Union(Matrix, SparseMatrixCSC)
#
  type Problem{T}
    conv::Conv.Converter
    params::Params
    r::Int
    c::Int
    A::MatrixTypesT(T)
  #
    Problem() = new()
  end
#
type Solution{T}
  solved::Bool
  L::MatrixTypesT(T)
  U::MatrixTypesT(T)
  p::Vector{Int}
  q::Vector{Int}
  rs::Vector{T}
  #
  function Solution(params::Params) x = new(); x.solved = false; return x; end
end
#
  function is_dense(prob::Problem) return issparse(prob.A) == false end
  function is_dense(sol::Solution) return issparse(sol.L) == false end
  function get_form(prob::Problem) return (is_dense(prob) ? :Dense : :Sparse) end
  function get_t(prob::Problem) return prob.conv.t end
  function get_r(prob::Problem) return prob.r end
  function get_c(prob::Problem) return prob.c end
  function comp_density(prob::Problem) return Shared_funcs.comp_density(prob.A) end
  function construct_Problem(params::Params, A::MatrixTypes)
    conv = Conv.converter(params["type"])
    prob = Problem{conv.t}()
    prob.conv = conv; prob.params = params;
    return prob
  end
  function fill_problem{T}(params::Params, prob::Problem{T}, A::MatrixTypes)
    prob.r, prob.c = size(A)
    prob.A = Conv.matrix(prob.conv, A)
  end
  function create_Problem(params::Params, A::MatrixTypes)
    prob = construct_Problem(params, A)
    fill_problem(params, prob, A)
    return prob
  end
  function construct_solution(type_t::DataType, params::Params)
    return Solution{type_t}(params)
  end
  function calc_solution_distance(prob::Problem, sol::Solution)
    if (sol.solved == false) return Inf end
    #println(dense(sol.L)); println(dense(sol.U));
    if (all(isfinite(sol.L)) == false) return Inf end
    if (all(isfinite(sol.U)) == false) return Inf end
    if (istril(sol.L) == false) return Inf end
    if (istriu(sol.U) == false) return Inf end
    nA = prob.A
    nA = isdefined(sol, :rs) ? scale(sol.rs, nA) : nA
    nA = nA[isdefined(sol, :p) ? sol.p : range(1,end),isdefined(sol, :q) ? sol.q : range(1,end)]
    normp = 1 #2-norm not yet implemented for sparse matrices
    #println(sol.p); println(sol.q); println(sol.rs);
    #println(dense(nA)); println(dense(sol.L*sol.U));
    return norm(nA - sol.L*sol.U, normp)
  end
end
