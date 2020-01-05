module Lp_glpk
  using Lp
  # Waiting for https://github.com/JuliaLang/julia/pull/6884
  importall JuMP
  importall GLPKMathProgInterface

  function install()
    Pkg.add("JuMP")
    Pkg.add("GLPK")
    Pkg.add("GLPKMathProgInterface")
  end

  type Dat{T}
    prob::Lp.Cf0_problem{T}

    Dat() = new()
  end

  function construct_dat(T::DataType) return Dat{T}() end

  function fill_dat{T}(prob::Lp.Cf0_problem, dat::Dat{T})
    dat.prob = prob
  end

  function solve_dat{T}(dat::Dat{T}, sol::Lp.Solution{T})
    m = Model(solver = GLPKSolverLP(method=:Exact))

    lp_prob = dat.prob
    @defVar(m, lp_x[1:lp_prob.n] >= 0 )
    for ri = 1:lp_prob.m
      aff = AffExpr(lp_x[1:lp_prob.n], reshape(lp_prob.A[ri,1:lp_prob.n], lp_prob.n), 0.0)
      @addConstraint(m, aff <= lp_prob.b[ri])
    end
    obj = AffExpr(lp_x[1:lp_prob.n], lp_prob.c, 0.0)
    setObjective(m, :Min, obj)

    @time status = JuMP.solve(m)

    status_dict = { :Optimal => :Optimal, :Unbounded => :Unbounded, :Infeasible => :Infeasible, :UserLimit => :Maxit, :Error => :Error, :NotSolved => :Created }
    sol.status = status_dict[status]
    sol.solved = (sol.status == :Optimal)
    sol.z = getObjectiveValue(m)
    if (sol.solved)
      sol.x = Array(lp_prob.conv.t, lp_prob.n)
      for i=1:lp_prob.n
        sol.x[i] = JuMP.getValue(lp_x[i])
      end
    end
  end

end
