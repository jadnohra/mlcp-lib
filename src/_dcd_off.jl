# Data Crunching Debugger

module Dcd

	type Session
		enable::Bool
		iters::Array{Dict{String,Any}, 1}

		Session(enable = false) = ( x = new(); x.enable = enable; x.iters = Array(Dict{String,Any}, 0); return x; )
	end

	macro set(sess, key, val) end
	macro it(sess) end
end
