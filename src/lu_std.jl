module Lu_elim_I
#
	using Lu
#
		ArgName = "elim_I"
		Descr = "Elimination without pivoting."
		Source = "[Trefethen] Numerical Linear Algebra (p151)."
		Detail = "O(n³)"
#
	type Dat{T}
		prob::Lu.Problem{T}
		U::Lu.MatrixTypesT(T)
		L::Matrix{T}
	#
		Dat() = new()
	end
#
	function construct_dat(T::DataType)
		return Dat{T}()
	end
	function fill_dat{T}(prob::Lu.Problem, dat::Dat{T})
		dat.prob = prob
		dat.U = deepcopy(prob.A); dat.L = eye(T, prob.r);
	end
	function solve_dat{T}(dat::Dat{T}, sol::Lu.Solution{T})
		U = dat.U; L = dat.L; m = dat.prob.c;
		for k in 1:m-1
			for j in k+1:m
				L[j,k] = U[j,k]/U[k,k]
				U[j,k:m] = U[j,k:m] - L[j,k] * U[k,k:m]
			end
		end
		sol.L = tril(L); sol.U = triu(U); sol.solved = all(isfinite(L)) && all(isfinite(U));
	end
end
