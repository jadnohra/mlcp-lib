
module Lp
	using Shared_funcs
	using Dcd
	using Conv

	:Created
	:Infeasible
	:Optimal
	:Unbounded
	:Maxit
	:Error

	:Cf0
	:Cf2b

	type SingleVar_translation{T}
		index::Int
		op1_mul::T
		op2_add::T
		SingleVar_translation() = (x = new(); index = -1; op1_mul = one(T); op2_add = zero(T); return x;)
	end

	type Cf0_problem{T} #Minimize
		conv::Conv.Converter
		n::Int
		m::Int
		c::Vector{T}
		A::Matrix{T}
		b::Vector{T}
		params::Params
		z_transl::SingleVar_translation{T}

		Cf0_problem() = new()
	end

	type Cf2b_problem{T} #Minimize
		cf0::Cf0_problem{T}
		l::Vector{T}
		h::Vector{T}

		Cf2b_problem() = new()
	end

	type Solution{T}
		solved::Bool
		status::Symbol
		z::T
		x::Vector{T}
		iter::Int

		dcd::Dcd.Session

		function Solution(params::Params)
			x = new()
			x.solved = false
			x.status = :Created
			x.dcd = Dcd.Session(get(params, "Dcd", false))
			x.iter = 0
			return x
		end
	end

	function get_form(prob::Cf0_problem) return :CF0 end
	function get_form(prob::Cf2b_problem) return :CF2b end
	function get_t(prob::Cf0_problem) return prob.conv.t end
	function get_t(prob::Cf2b_problem) return get_t(prob.cf0) end
	function get_n(prob::Cf0_problem) return prob.n end
	function get_m(prob::Cf0_problem) return prob.m end
	function get_n(prob::Cf2b_problem) return get_n(prob.cf0) end
	function get_m(prob::Cf2b_problem) return get_m(prob.cf0) end

	function comp_density(prob::Cf0_problem) return Shared_funcs.comp_density(prob.A[1:prob.m, 1:prob.n]) end
	function comp_density(prob::Cf2b_problem) return comp_density(prob.cf0) end

	function construct_solution(type_t::DataType, params::Params)
		return Solution{type_t}(params)
	end

	function transl_single{T}(transl::SingleVar_translation{T}, v::T)
		return (v*transl.op1_mul)+transl.op2_add
	end

	function translate_solution{T}(prob::Cf0_problem{T}, raw_sol::Solution{T})
		sol = construct_solution(prob.conv.t, prob.params)
		sol.solved = raw_sol.solved
		sol.status = raw_sol.status
		sol.dcd = raw_sol.dcd
		sol.iter = raw_sol.iter
		if (sol.solved)
			sol.z = transl_single(prob.z_transl, raw_sol.z)
			sol.x = raw_sol.x
		else
			sol.z = zero(T)
			sol.x = Array(T, 0)
		end

		return sol
	end

	function translate_solution{T}(prob::Cf2b_problem{T}, raw_sol::Solution{T})
		return translate_solution(prob.cf0, raw_sol)
	end

	function construct_Cf0_problem(params::Params)
		conv = Conv.converter(params["type"])
		prob = Cf0_problem{conv.t}()
		prob.conv = conv
		prob.z_transl = SingleVar_translation{prob.conv.t}()
		prob.params = params
		return prob
	end

	function construct_Cf2b_problem(params::Params)
		cf0 = construct_Cf0_problem(params)
		prob = Cf2b_problem{cf0.conv.t}()
		prob.cf0 = cf0
		return prob
	end

	# min: z = cx, subj: Ax <= b, x >= 0
	function fill_problem{T}(params::Params, prob::Cf0_problem{T}, c::Vector, A::Matrix, b::Vector)
		prob.n = length(c)
		prob.m = length(b)
		prob.c = vcat( Conv.vector(prob.conv, c), zeros(T, prob.m) )
		prob.A = hcat( Conv.matrix(prob.conv, A), Conv.matrix(prob.conv, eye(prob.m)) )
		prob.b = Conv.vector(prob.conv, b)
	end

	# min: z = cx, subj: Ax <= b, x >= 0
	function create_min_Cf0_problem(params::Params, c::Vector, A::Matrix, b::Vector)
		prob = construct_Cf0_problem(params)
		fill_problem(params, prob, c, A, b)
		return prob
	end

	# max: z = cx, subj: Ax <= b, x >= 0
	function create_max_Cf0_problem(params::Params, c::Vector, A::Matrix, b::Vector)
		mone = Conv.conv(params["type"], -1)
		prob = create_min_Cf0_problem(params, mone*(c), A, b)
		prob.z_transl.op1_mul = mone
		return prob
	end

	# min: z = cx, subj: zl <= Ax <= zh, xl <= x <= xh
	function fill_problem{T}(params::Params, prob::Cf2b_problem{T}, c::Vector, A::Matrix,
			zl::Vector, zh::Vector, xl::Vector, xh::Vector)
		fill_problem(params, prob.cf0, c, A, Array(T, 0))
		prob.l = vcat( Conv.vector(prob.cf0.conv, xl), Conv.vector(prob.cf0.conv, zl) )
		prob.h = vcat( Conv.vector(prob.cf0.conv, xh), Conv.vector(prob.cf0.conv, zh) )
	end

	# min: z = cx, subj: zl <= Ax <= zh, xl <= x <= xh
	function create_min_Cf2b_problem(params::Params, c::Vector, A::Matrix,
			zl::Vector, zh::Vector, xl::Vector, xh::Vector)
		prob = construct_Cf2b_problem(params)
		fill_problem(params, prob, c, A, zlh, xlh)
		return prob
	end

	# max: z = cx, subj: zl <= Ax <= zh, xl <= x <= xh
	function create_max_Cf2b_problem(params::Params, c::Vector, A::Matrix,
			zl::Vector, zh::Vector, xl::Vector, xh::Vector)
		mone = Conv.conv(params["type"], -1)
		prob = create_min_Cf2b_problem(params, mone*(c), A, zlh, xlh)
		prob.cf0.z_transl.op1_mul = mone
		return prob
	end

	function dcd_var_info(n, i) return i <= n ? ("x", i) : ("w", i-n); end
	function dcd_var(n, i) s,i = dcd_var_info(n, i); return @sprintf "%s%d" s i; end
	function dcd_iR(n, m, iB)
		iR = Int[]
		for i = 1:m+n
			if findfirst(iB, i) == 0 push!(iR, i) end
		end
		return iR
	end

	function dcd_pivots_impl(sol, typ)

		function dcd_print_basis(n, iB)
			for i = 1:length(iB)
				vs,vi = dcd_var_info(n, iB[i])
				print(vs, vi)
			end
			@printf "\n"
		end

		function dcd_print_pivot_type(n, iB, iR, typ)
			if (typ == 2) dcd_print_basis(n, iB)
			elseif (typ == 3) dcd_print_basis(n, iR)
			elseif (typ == 4) println(sort(iB))
			elseif (typ == 5) println(sort(iR))
			end
		end

		function dcd_print_pivot(i, iter, n, iB, iR, typ)
			if (haskey(iter, "pivot"))
				c,r = iter["pivot"]
				if (typ == 1) println( i, ". ",  dcd_var(n, iR[c]), ",", dcd_var(n, iB[r]) ) end
				iR[c], iB[r] = iB[r], iR[c]
				if (typ != 1) dcd_print_pivot_type(n, iB, iR, typ) end
			end
		end

		if isempty(sol.dcd.iters) return; end
		n = sol.prob.n
		iB = deepcopy(sol.dcd.iters[1]["iB"])
		iR = dcd_iR(sol.prob.n, sol.prob.m, iB)
		if (typ != 1) dcd_print_pivot_type(n, iB, iR, typ) end
		for i = 1:length(sol.dcd.iters)
			iter = sol.dcd.iters[i]
			dcd_print_pivot(i, iter, n, iB, iR, typ)
		end
	end

	function dcd_pivots(sol) dcd_pivots_impl(sol, 1) end
	function dcd_basis(sol) dcd_pivots_impl(sol, 2) end
	function dcd_nbasis(sol) dcd_pivots_impl(sol, 3) end
	function dcd_ibasis(sol) dcd_pivots_impl(sol, 4) end
	function dcd_inbasis(sol) dcd_pivots_impl(sol, 5) end
	function dcd_iters(sol) println(length(sol.dcd.iters)) end
	function dcd_key(sol, key)
		for i = 1:length(sol.dcd.iters)
			iter = sol.dcd.iters[i]
			if (haskey(iter, key))
				val = iter[key]
				println(i, ". ", val)
			end
		end
	end
end
