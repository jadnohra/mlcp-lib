module Shared_funcs
#
  # Waiting for https://github.com/JuliaLang/julia/pull/6884
  using GZip
#
  function format_percent(v)
    return strip(strip(@sprintf("%0.2f", 100.0 * v), '0'), '.') * "%"
  end
  function comp_density(A::AbstractMatrix)
    return countnz(A) / length(A)
  end
  function _read_matrix_mtx(filename, infoonly::Bool=false)
    function skewsymmetric!(M::AbstractMatrix)
        m,n = size(M)
        m == n || throw(DimensionMismatch())
        for i=1:n, j=1:n
            if M[i,j] != 0
                M[j,i] = -M[i,j]
            end
        end
        return M
    end
    function symmetric!(M::AbstractMatrix)
        m,n = size(M)
        m == n || throw(DimensionMismatch())
        for i=1:n, j=1:n
            if M[i,j] != 0
                M[j,i] = M[i,j]
            end
        end
        return M
    end
    function hermitian!(M::AbstractMatrix)
        m,n = size(M)
        m == n || throw(DimensionMismatch())
        for i=1:n, j=1:n
            if M[i,j] != 0
                M[j,i] = conj(M[i,j])
            end
        end
        return M
    end
    #println(pwd())
    mmfile = Nothing()
    if (endswith(filename, ".gz"))
      mmfile = GZip.open(filename,"r")
    else
      mmfile = open(filename,"r")
    end
    # Read first line
    firstline = chomp(readline(mmfile))
    tokens = split(firstline)
    if length(tokens) != 5
      throw(ParseError(string("Not enough words on first line: ", ll)))
    end
    if tokens[1] != "%%MatrixMarket"
      throw(ParseError(string("Not a valid MatrixMarket header:", ll)))
    end
    (head1, rep, field, symm) = map(lowercase, tokens[2:5])
    if head1 != "matrix"
        throw(ParseError("Unknown MatrixMarket data type: $head1 (only \"matrix\" is supported)"))
    end
    eltype = field == "real" ? Float64 :
             field == "complex" ? Complex128 :
             field == "pattern" ? Bool :
             throw(ParseError("Unsupported field $field (only real and complex are supported)"))

    symlabel = symm == "general" ? identity :
               symm == "symmetric" ? symmetric! :
               symm == "hermitian" ? hermitian! :
               symm == "skew-symmetric" ? skewsymmetric! :
               throw(ParseError("Unknown matrix symmetry: $symm (only general, symmetric, skew-symmetric and hermitian are supported)"))
    # Skip all comments and empty lines
    ll   = readline(mmfile)
    while length(chomp(ll))==0 || (length(ll) > 0 && ll[1] == '%')
        ll = readline(mmfile)
    end
    # Read matrix dimensions (and number of entries) from first non-comment line
    dd = map(parseint, split(ll))
    if length(dd) < (rep == "coordinate" ? 3 : 2)
        throw(ParseError(string("Could not read in matrix dimensions from line: ", ll)))
    end
    rows = dd[1]
    cols = dd[2]
    entries = (rep == "coordinate") ? dd[3] : (rows * cols)
    if infoonly
        return (rows, cols, entries, rep, field, symm)
    end
    if rep == "coordinate"
        rr = Array(Int, entries)
        cc = Array(Int, entries)
        xx = Array(eltype, entries)
        for i in 1:entries
            flds = split(readline(mmfile))
            rr[i] = parseint(flds[1])
            cc[i] = parsefloat(flds[2])
            if eltype == Complex128
                xx[i] = Complex128(parsefloat(flds[3]), parsefloat(flds[4]))
            elseif eltype == Float64
                xx[i] = parsefloat(flds[3])
            else
                xx[i] = true
            end
        end
        return symlabel(sparse(rr, cc, xx, rows, cols))
    end
    return symlabel(reshape([parsefloat(readline(mmfile)) for i in 1:entries], (rows,cols)))
  end
  function read_matrix_mtx(filename)
    return _read_matrix_mtx(filename)
  end
  function random_matrix(seed, scale, dense, n, m, show)
    function toindex(x,n) return clamp(int((n)* 0.5*(1.0+x)), 1, n) end
    rng = MersenneTwister(seed)
    B = A = scale * rand(rng, (m, n))
    if dense != 100
      #zeros = rand(rng, 1:(m*n), zero_n) # Needs next version of Julia
      zero_n = int(round((m*n) * (1.0 - dense / 100.0)))
      if (zero_n > 0)
        nz_i = [1:m*n]
        while length(nz_i) + zero_n > (m*n)
          splice!(nz_i, toindex(rand(rng), length(nz_i)))
        end
        B = zeros(Float64, (m,n))
        for i in nz_i B[i] = A[i] end
      end
    end
    if (show)
      println("A", B)
    end
    return B
  end
  function random_matrix(params)
    seed = int(get(params, "seed", int(time_ns()) % 32768))
    scale = int(get(params, "scale", 1.0))
    dense = int(get(params, "dense", 100 ))
    n = int(get(params, "n", 10))
    m = int(get(params, "m", n))
    show = get(params, "show", false)
    println("/seed:", seed)
    return random_matrix(seed, scale, dense, n, m, show)
  end
  function str_range(str, list = [])
  	function single_range(str)
      beg = int(split(str, "-")[1])
      nd = int(split(split(str, "-")[2],":")[1])
      stp = contains(str,":") ? int(split(str, ":")[2]) : 1
      cnt = int((nd-beg+1)/stp)
    	return [range(beg,stp,cnt)]
    end
    if str == "all" return 1:length(list) end
    rngs = Array(Int, 0)
    for x in split(str, ",")
      if length(x) != 0 rngs = vcat(rngs, contains(x, "-") ? single_range(x) : [int(x)]) end
    end
    return rngs
  end
  function print_choose(list, autochoose::Bool = false)
    if length(list) == 0 return [] end
    if length(list) == 1 && autochoose return list end
    cols = [:white, :yellow]
    for i in 1:length(list)
      print_with_color(cols[1+i%2], " $i. $(list[i]) \n")
    end
    print(" > ")
    return str_range(strip(readline()), list)
  end
end
