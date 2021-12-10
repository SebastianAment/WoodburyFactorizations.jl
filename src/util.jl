function Base.Matrix(F::BunchKaufman)
	U = F.U[:, F.ipiv]
	*(U, F.D, U')
end
Base.AbstractMatrix(F::BunchKaufman) = Matrix(F)
function LinearAlgebra.rdiv!(X::AbstractMatrix, F::BunchKaufman)
	# ldiv!(F, X')' # doesn't work since X' does not count as strided vector (requirement for ldiv!)
	X .= adjoint(F \ X')
end
LinearAlgebra.logdet(B::Bidiagonal) = logabsdet(B)[1]
function LinearAlgebra.logabsdet(B::Bidiagonal)
	d = @view B[diagind(B)]
	D = Diagonal(d)
	logabsdet(D)
end
