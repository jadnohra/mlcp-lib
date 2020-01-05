# István Maros. Computational Techniques of the Simplex Method. Kluwer Academic Publishers, Boston, 2003. [CTSM]
# Robert J. Vanderbei. Linear Programming: Foundations and Extensions. Springer, second edition, 2001. [LPFE]
# Suhl, Suhl. A fast LU update for linear programming. [LULP]

#=
	List of algorithms

	Linear Programming
		I. Solve, Revised, Canonical Simplex.
			1. revised simplex, basis inversion, no reuse
			2. revised simplex, basis factorization, no reuse
			3. revised simplex, basis factorization, reuse, no re-factorization
			4. revised simplex, basis factorization, reuse, re-factorization
			5. revised simplex, basis factorization, reuse, re-factorization, degen (CTSM p230), sing
			6. 5. with 'textbook' Phase-I
			7. 5. with Maros' Phase-I (p203)

		II. Solve, Revised, Sparse, Canonical Simplex.
			1. revised sparse simplex

		III. Generalized Simplex variants of I, then II (Presolve, Postsolve)
			1. generalized version of I.1, using Maros' CF-2 (p13-17), but without translation, calling it CF-2b. (Handwritten p29).
			2. 1. with variable transformation for less lohi calculations, the proper CF-2 formulation
			3. 2. with infeasibility tolerance, type-0 handling and Harris ratio. (CSTM p178)

	LU Factorization
		I. Suhl-suhl (CTSM p136), and original [LULP]
=#

module Jad

module lp_I_1
	using Dcd
	using Lp

	type Dat{T}
		prob::Lp.Cf0_problem
		sol::Lp.Solution
		maxit::Int
		n::Int
		m::Int
		status::Symbol

		# Using CTSM notation.
		# LPFE notation translation.
		# αq - ΔxB, β - xB*, q - j, i - p, t - θ

		iB::Array{Int, 1}
		iR::Array{Int, 1}
		B::Array{T, 2}
		Binv::Array{T, 2}
		β::Array{T, 1}
		π::Array{T, 1}
		cBT::Array{T, 2}
		dJ::Array{T, 1}
		αq::Array{T, 1}
		z::T
		zero::T

		dcd::Dcd.Session
		it::Int
		q::Int
		p::Int
		θ::T

		Dat() = new()
	end

	function construct_dat(T::DataType) return Dat{T}() end

	function fill_dat{T}(prob::Lp.Cf0_problem, dat::Dat{T})
		n = prob.n; m = prob.m;
		dat.maxit = get(prob.params, "maxit", 0)
		dat.prob = prob
		dat.n = n
		dat.m = m
		dat.iB = Array(Int, m)
		dat.iR = Array(Int, n)
		dat.B = Array(T, (m,m))
		#dat.Binv = Array(T, (m,m))
		dat.β = Array(T, m)
		#dat.π = Array(T, m+n)
		#dat.cBT = Array(T, m)
		dat.dJ = Array(T, n)
		#dat.αq = Array(T, n)
		dat.zero = zero(T)
	end

	function warn{T}(dat::Dat{T}, warning::Symbol)
		println("Warning:", warning)
	end

	function succeed{T}(dat::Dat{T})
		dat.status = :Ended
		sol = dat.sol
		sol.iter = dat.it
		sol.solved = true
		sol.status = :Optimal
		sol.z = dat.z
		sol.x = zeros(dat.prob.conv.t, dat.n)
		for i = 1:length(dat.iB)
			xi = dat.iB[i]
			if (xi <= dat.n) sol.x[xi] = dat.β[i]; end
		end
	end

	function fail{T}(dat::Dat{T}, status::Symbol)
		dat.status = :Ended
		sol = dat.sol
		sol.iter = dat.it
		sol.solved = false
		sol.status = status
	end

	function sel_cols{T}(src::Matrix{T}, dst::Matrix{T}, cols::Vector{Int})
		for i = 1:length(cols)
			c = cols[i]
			for r = 1:size(dst)[1]
				dst[r,i] = src[r,c]
			end
		end
	end

	function set_basis_logical{T}(dat::Dat{T})
		dat.iB = [i for i in dat.n+1:dat.n+dat.m]
		dat.iR = [i for i in 1:dat.n]
	end

	function comp_B_R{T}(dat::Dat{T}) sel_cols(dat.prob.A, dat.B, dat.iB) end
	function comp_Binv{T}(dat::Dat{T}) dat.Binv = inv(dat.B) end
	function comp_cBT{T}(dat::Dat{T}) dat.cBT = transpose(dat.prob.c[dat.iB]) end
	function init_β{T}(dat::Dat{T}) dat.β = dat.Binv * (dat.prob.b) end
	function update_β{T}(dat::Dat{T}, p::Int, θ::T) dat.β = dat.β - (θ*dat.αq); dat.β[p] = θ; end
	function init_z{T}(dat::Dat{T}) dat.z = dot(reshape(dat.cBT, length(dat.cBT)), dat.β); end
	function update_z{T}(dat::Dat{T}, q::Int, θ::T) dat.z = dat.z + θ * dat.dJ[q]; end
	function comp_π{T}(dat::Dat{T}) dat.π = reshape(dat.cBT * dat.Binv, length(dat.cBT)); end
	function calc_dj{T}(dat::Dat{T}, j::Int) i = dat.iR[j]; dj = dat.prob.c[i] - dot(dat.π, dat.prob.A[:,i]); return dj; end
	function comp_dJ{T}(dat::Dat{T}) for j = 1:dat.n dat.dJ[j] = calc_dj(dat, j) end; end
	function comp_αq{T}(dat::Dat{T}, q::Int) aq = dat.iR[q]; dat.αq = dat.Binv * (dat.prob.A[:,aq]); end
	function check_optimal_dJ{T}(dat::Dat{T}) return all( dj->(dj >= dat.zero), dat.dJ ); end
	function check_feasible_β{T}(dat::Dat{T}) return all( β->(β >= dat.zero), dat.β ); end

	function pivot_iB_iR{T}(dat::Dat{T}, q::Int, p::Int) dat.iR[q], dat.iB[p] = dat.iB[p], dat.iR[q] end

	function price_full_dantzig{T}(dat::Dat{T})
		# [CTSM].p187
		min_i = 0; min_dj = dat.zero;
		for i = 1:length(dat.dJ)
			dj = dat.dJ[i]
			if (dj < dat.zero && dj <= min_dj)
				if (dj == min_dj)
					warn(dat, :Degen)
				end
				min_i = i; min_dj = dj;
			end
		end
		return min_i
	end

	function chuzro{T}(dat::Dat{T})
		min_i = 0; min_ratio = dat.zero;
		for i = 1:length(dat.αq)
			if (dat.αq[i] > dat.zero)
				ratio = dat.β[i] / dat.αq[i]
				if (min_i == 0 || ratio <= min_ratio)
					if (ratio == min_ratio)
						warn(dat, :Degen)
					end
					min_i = i; min_ratio = ratio;
				end
			end
		end
		return (min_i, min_ratio)
	end


	function solve_init{T}(dat::Dat{T}, sol::Lp.Solution{T})
		dat.sol = sol
		dat.dcd = sol.dcd
		Dcd.@it(dcd)
		dat.it = 0
		dat.status = :Iterating
	end

	function solve_step0{T}(dat::Dat{T})
		set_basis_logical(dat)
		Dcd.@set(dcd, "iB", dat.iB)
		#Initializations
		comp_cBT(dat); comp_B_R(dat); comp_Binv(dat);
		init_β(dat); init_z(dat);
		Dcd.@set(dcd, "β0", dat.β)
	end

	function solve_step_phaseI{T}(dat::Dat{T})
		if (check_feasible_β(dat) == false)
			warn(degen, :PhaseI)
			fail(dat, sol, :Infeasible)
		end
	end

	function solve_iter{T}(dat::Dat{T})
		return (dat.status == :Iterating && (dat.maxit == 0 || dat.it < dat.maxit))
	end

	function solve_next{T}(dat::Dat{T})
		dat.it = dat.it + 1
	end

	function solve_step1{T}(dat::Dat{T})
		Dcd.@it(dcd)
		#Step 1
		Dcd.@set(dcd, "B", dat.B); Dcd.@set(dcd, "Binv", dat.Binv);
		Dcd.@set(dcd, "cBT", dat.cBT);
		comp_π(dat); Dcd.@set(dcd, "π", dat.π);
	end

	function solve_step_price{T}(dat::Dat{T})
		comp_dJ(dat); Dcd.@set(dcd, "dJ", dat.dJ);
		if check_optimal_dJ(dat) succeed(dat); return; end
		dat.q = price_full_dantzig(dat); Dcd.@set(dcd, "q", dat.q);
	end

	function solve_step3{T}(dat::Dat{T})
		comp_αq(dat, dat.q); Dcd.@set(dcd, "αq", dat.αq);
	end

	function solve_step_chuzro{T}(dat::Dat{T})
		dat.p, dat.θ = chuzro(dat); Dcd.@set(dcd, "p", dat.p); Dcd.@set(dcd, "θ", dat.θ);
		if (dat.p == 0) fail(dat, :Unbounded); return; end
		Dcd.@set(dcd, "pivot", (dat.q, dat.p))
	end

	function solve_step5{T}(dat::Dat{T})
		pivot_iB_iR(dat, dat.q, dat.p)
	end

	function solve_update{T}(dat::Dat{T})
		update_β(dat, dat.p, dat.θ); Dcd.@set(dcd, "β", dat.β);
		update_z(dat, dat.q, dat.θ); Dcd.@set(dcd, "z", dat.z);
		comp_cBT(dat); comp_B_R(dat); comp_Binv(dat);
	end

	function solve_dat{T}(dat::Dat{T}, sol::Lp.Solution{T})
		# [CTSM].p33,p30, but with naive basis reinversion instead of update.
		solve_init(dat, sol)
		solve_step0(dat)
		solve_step_phaseI(dat)

		while solve_iter(dat)
			solve_step1(dat)
			solve_step_price(dat);
			if (dat.status != :Iterating) break; end
			solve_step3(dat)
			solve_step_chuzro(dat); if (dat.status != :Iterating) break; end
			solve_step5(dat)
			solve_update(dat)
			solve_next(dat)
		end
		if (dat.status == :Iterating) fail(dat, :Maxit); end
	end

end

module lp_III_1
	using Dcd
	using Lp
	using Jad.lp_I_1

	type Dat{T}
		Base::Module
		prob::Lp.Cf2b_problem{T}
		bdat::lp_I_1.Dat{T}
		is_hi::Vector{Bool}

		p_hi::Bool

		Dat() = new()
	end

	function construct_dat(T::DataType) return Dat{T}() end

	function fill_dat{T}(prob::Lp.Cf2b_problem{T}, dat::Dat{T})
		dat.Base = lp_I_1
		dat.prob = prob
		dat.bdat = dat.Base.construct_dat(T)
		Base.fill_dat(prob.base, dat.bdat)
		dat.is_hi = zeros(Bool, dat.prob.n+dat.prob.m) # all variables at low.
	end

	function init_β{T}(dat::Dat{T}, xlo::Vector{T})
		R = dat.base.prob.A[:,dat.iR]
		dat.β = dat.Binv * (dat.prob.b - R*xlo)
	end

	function init_z{T}(dat::Dat{T}, xlo::Vector{T})
		cR = dat.prob.c[dat.iR]
		zR = dot(cR, xlo)
		zB = dot(reshape(dat.cBT, length(dat.cBT)), dat.β)
		dat.z = zB + zR
	end

	function solve_step0{T}(sdat::Dat{T})
		Base = dat.Base; dat = sdat.bdat;
		Base.set_basis_logical(dat)
		Dcd.@set(dcd, "iB", dat.iB)
		#Initializations
		Base.comp_cBT(dat); Base.comp_B_R(dat); Base.comp_Binv(dat);
		xlo = dat.prob.l[dat.iR]
		init_β(sdat, xlo); init_z(sdat, xlo);
		Dcd.@set(dcd, "β0", dat.β)
	end

	function price_full_dantzig{T}(sdat::Dat{T})
		# [CTSM].p187
		Base = sdat.Base; dat = sdat.bdat;
		min_i = 0; min_dj = dat.zero;
		for i = 1:length(dat.dJ)
			dj = dat.dJ[i]
			if (dj < dat.zero && dj <= min_dj && (sdat.is_hi[dat.iB[dj]] == false) )
				if (dj == min_dj)
					Base.warn(dat, :Degen)
				end
				min_i = i; min_dj = dj;
			end
		end
		return min_i
	end

	function solve_step_price{T}(sdat::Dat{T})
		Base = sdat.Base; dat = sdat.bdat;
		Base.comp_dJ(dat); Dcd.@set(dcd, "dJ", dat.dJ);
		if Base.check_optimal_dJ(dat) Base.succeed(dat); return; end
		dat.q = price_full_dantzig(sdat); Dcd.@set(dcd, "q", dat.q);
	end

	function chuzro{T}(dat::Dat{T})
		# Handwritten notes (p29)
		glob_i = 0; glob_max_θ = dat.zero; glob_hi = false;
		for i = 1:length(dat.αq)
			set_hi = false
			if (dat.αq[i] > dat.zero)
				lo = dat.prob.l[dat.iB[i]]
				max_θ = (dat.β[i] - lo) / dat.αq[i]
				set_hi = false
			elseif (dat.αq[i] < dat.zero)
				hi = dat.prob.h[dat.iB[i]]
				max_θ = (dat.β[i] - hi) / dat.αq[i]
				set_hi = true
			else
				continue
			end
			if (glob_i == 0 || max_θ >= glob_max_θ)
				if (max_θ == globa_max_θ) warn(dat, :Degen); end
				glob_i = i; glob_max_θ = max_θ; glob_hi = set_hi;
			end
		end

		return (glob_i, glob_max_θ, glob_hi)
	end

	function solve_step_chuzro{T}(dat::Dat{T})
		dat.p, dat.θ, dat.p_hi = chuzro(dat); Dcd.@set(dcd, "p", dat.p); Dcd.@set(dcd, "θ", dat.θ);
		if (dat.p == 0) fail(dat, :Unbounded); return; end
		Dcd.@set(dcd, "pivot", (dat.q, dat.p))
	end

	function solve_update{T}(dat::Dat{T})
		dat.Base.solve_update(dat.bdat)
		dat.is_hi[dat.p] = dat.p_hi
	end

	function solve_dat{T}(sdat::Dat{T}, sol::Lp.Solution{T})
		# lp_I_1, with generalized forumulation (w/o translation)
		Base = sdat.Base;	bdat = sdat.bdat;
		Base.solve_init(bdat, sol)
		solve_step0(dat)
		Base.solve_step_phaseI(bdat)

		while solve_iter(bdat)
			Base.solve_step1(bdat)
			solve_step_price(dat); if (bdat.status != :Iterating) break; end
			Base.solve_step3(bdat)
			solve_step_chuzro(dat); if (bdat.status != :Iterating) break; end
			Base.solve_step5(bdat)
			solve_update(dat)
			Base.solve_next(bdat)
		end
		if (dat.status == Iterating) Base.fail(bdat, :Maxit); end
	end

end

module lu_I
	#function factorize{}(B::Array{T, 2})
	#end
end

#=
module lp_I_2
	using Dcd
	using Lp

	type Dat{T}
		prob::Lp.Cf0_problem
		maxit::Int
		n::Int
		m::Int
		# Using CTSM notation.
		iB::Array{Int, 1}
		iR::Array{Int, 1}
		B::Array{T, 2}
		BL::Array{T, 2}
		BU::Array{T, 2}
		BP::Array{T, 1}
		β::Array{T, 1}
		π::Array{T, 1}
		cBT::Array{T, 2}
		dJ::Array{T, 1}
		αq::Array{T, 1}
		z::T
		zero::T

		# LPFE notation translation.
		# αq - ΔxB, β - xB*, q - j, i - p, t - θ

		Dat() = new()
	end

	function fill_dat{T}(prob::Lp.Cf0_problem, dat::Dat{T})
		n = prob.n; m = prob.m;
		dat.maxit = get(prob.params, "maxit", 0)
		dat.prob = prob
		dat.n = n
		dat.m = m
		dat.iB = Array(Int, m)
		dat.iR = Array(Int, n)
		dat.B = Array(T, (m,m))
		#dat.Binv = Array(T, (m,m))
		dat.β = Array(T, m)
		#dat.π = Array(T, m+n)
		#dat.cBT = Array(T, m)
		dat.dJ = Array(T, n)
		#dat.αq = Array(T, n)
		dat.zero = zero(T)
	end

	function sel_cols{T}(src::Matrix{T}, dst::Matrix{T}, cols::Vector{Int})
		for i = 1:length(cols)
			c = cols[i]
			for r = 1:size(dst)[1]
				dst[r,i] = src[r,c]
			end
		end
	end

	function set_basis_logical{T}(dat::Dat{T})
		dat.iB = [i for i in dat.n+1:dat.n+dat.m]
		dat.iR = [i for i in 1:dat.n]
	end

	function comp_B_R{T}(dat::Dat{T}) sel_cols(dat.prob.A, dat.B, dat.iB) end
	function comp_BLU{T}(dat::Dat{T}) dat.BL, dat.BU, dat.BP = lu(dat.B) end
	function comp_cBT{T}(dat::Dat{T}) dat.cBT = transpose(dat.prob.c[dat.iB]) end
	function init_β{T}(dat::Dat{T}) dat.β = dat.Binv * (dat.prob.b) end
	function update_β{T}(dat::Dat{T}, p::Int, θ::T) dat.β = dat.β - (θ*dat.αq); dat.β[p] = θ; end
	function init_z{T}(dat::Dat{T}) dat.z = dot(reshape(dat.cBT, length(dat.cBT)), dat.β); end
	function update_z{T}(dat::Dat{T}, q::Int, θ::T) dat.z = dat.z + θ * dat.dJ[q]; end
	function comp_π{T}(dat::Dat{T}) dat.π = reshape(dat.cBT * dat.Binv, length(dat.cBT)); end
	function calc_dj{T}(dat::Dat{T}, j::Int) i = dat.iR[j]; dj = dat.prob.c[i] - dot(dat.π, dat.prob.A[:,i]); return dj; end
	function comp_dJ{T}(dat::Dat{T}) for j = 1:dat.n dat.dJ[j] = calc_dj(dat, j) end; end
	function comp_αq{T}(dat::Dat{T}, q::Int) aq = dat.iR[q]; dat.αq = dat.Binv * (dat.prob.A[:,aq]); end
	function check_optimal_dJ{T}(dat::Dat{T}) return all( dj->(dj >= dat.zero), dat.dJ ); end
	function check_feasible_β{T}(dat::Dat{T}) return all( β->(β >= dat.zero), dat.β ); end

	function pivot_iB_iR{T}(dat::Dat{T}, q::Int, p::Int) dat.iR[q], dat.iB[p] = dat.iB[p], dat.iR[q] end

	function price_full_dantzig{T}(dat::Dat{T})
		# [CTSM].p187
		min_i = 0; min_dj = dat.zero;
		for i = 1:length(dat.dJ)
			dj = dat.dJ[i]
			if (dj < dat.zero && dj <= min_dj)
				if (dj == min_dj)
					println("Warning: degeneracy.")
				end
				min_i = i; min_dj = dj;
			end
		end
		return min_i
	end

	function chuzro{T}(dat::Dat{T})
		min_i = 0; min_ratio = dat.zero;
		for i = 1:length(dat.αq)
			if (dat.αq[i] > dat.zero)
				ratio = dat.β[i] / dat.αq[i]
				if (min_i == 0 || ratio <= min_ratio)
					if (ratio == min_ratio)
						println("Warning: degeneracy.")
					end
					min_i = i; min_ratio = ratio;
				end
			end
		end
		return (min_i, min_ratio)
	end

	function succeed{T}(dat::Dat{T})
		sol.solved = true
		sol.status = :Optimal
		sol.z = dat.z
		sol.x = eval(parse( "zeros($(dat.prob.type_s), $(dat.n))" ))
		for i = 1:length(dat.iB)
			xi = dat.iB[i]
			if (xi <= dat.n) sol.x[xi] = dat.β[i]; end
		end
	end

	function fail(dat, sol::Lp.Solution, status::Symbol) sol.solved = false; sol.status = status; end

	function solve_dat{T}(dat::Dat{T}, sol::Lp.Solution{T})
		# [CTSM].p33,p30, but with naive basis reinversion instead of update.
		dcd = sol.dcd
		Dcd.@it(dcd)
		it = 0
		#Step 0
		set_basis_logical(dat)
		Dcd.@set(dcd, "iB", dat.iB)
		#Initializations
		comp_cBT(dat); comp_B_R(dat); comp_Binv(dat);
		init_β(dat); init_z(dat);
		Dcd.@set(dcd, "β0", dat.β)
		#todo phaseI
		if (check_feasible_β(dat) == false) println("Warning: phaseI."); fail(dat, sol, :Infeasible); return; end

		while(dat.maxit == 0 || it < dat.maxit)
			Dcd.@it(dcd)
			#Step 1
			Dcd.@set(dcd, "B", dat.B); Dcd.@set(dcd, "Binv", dat.Binv);
			Dcd.@set(dcd, "cBT", dat.cBT);
			comp_π(dat); Dcd.@set(dcd, "π", dat.π);
			#Step 2
			comp_dJ(dat); Dcd.@set(dcd, "dJ", dat.dJ);
			if check_optimal_dJ(dat) succeed(dat, sol); return; end
			q = price_full_dantzig(dat); Dcd.@set(dcd, "q", q);
			#Step 3
			comp_αq(dat, q); Dcd.@set(dcd, "αq", dat.αq);
			#Step 4
			p,θ = chuzro(dat); Dcd.@set(dcd, "p", p); Dcd.@set(dcd, "θ", θ);
			if (p == 0) fail(dat, sol, :Unbounded); return; end
			Dcd.@set(dcd, "pivot", (q, p))
			#Step 5
			pivot_iB_iR(dat, q, p)
			#Updates
			update_β(dat, p, θ); Dcd.@set(dcd, "β", dat.β);
			update_z(dat, q, θ); Dcd.@set(dcd, "z", dat.z);
			comp_cBT(dat); comp_B_R(dat); comp_Binv(dat);
			it = it + 1
		end
		fail(dat, sol, :Maxit)
	end


end
=#

end #module Jad
