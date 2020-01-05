module Lu_julia_dense
#
	using Lu
#
		ArgName = "julia_dense"
		Descr = "Julia's default algorithm on a dense matrix."
#
	type Dat{T}
		prob::Lu.Problem{T}
		A::Array{T,2}
	#
		Dat() = new()
	end
#
	function construct_dat(T::DataType)
		return Dat{T}()
	end
	function fill_dat{T}(prob::Lu.Problem, dat::Dat{T})
		dat.prob = prob
		dat.A = full(prob.A)
	end
	function solve_dat{T}(dat::Dat{T}, sol::Lu.Solution{T})
		sol.L, sol.U, sol.p = lu(dat.A)
		sol.solved = true
	end
end
#
#
module Lu_julia_sparse
#
	using Lu
#
	ArgName = "julia_sparse"
	Descr = "Julia's default algorithm on a sparse matrix (uses umfpack)."
#
	type Dat{T}
		prob::Lu.Problem{T}
		A::SparseMatrixCSC{T}
	#
		Dat() = new()
	end
#
	function construct_dat(T::DataType)
		return Dat{T}()
	end
	function fill_dat{T}(prob::Lu.Problem, dat::Dat{T})
		dat.prob = prob
		dat.A = sparse(prob.A)
	end
	function solve_dat{T}(dat::Dat{T}, sol::Lu.Solution{T})
		# Note, lufact with sparse matrices seems to only be possible with a specific type (e.g Float64)
		# try: methods(lufact)
		# https://github.com/JuliaLang/julia/issues/4439 says: "I do believe that UMFPACK does not have a Float32/Complex64 interface"
		luF = lufact(dat.A); L,U,p,q,Rs = luF[:(:)]
		sol.L = L; sol.U = U; sol.p = p; sol.q = q; sol.rs = Rs; sol.solved = true;
	end
end
