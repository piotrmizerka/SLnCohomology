# The functions here are adapted from https://git.rwth-aachen.de/jens.brandt/plesken-souvignier/-/tree/master/src

using LinearAlgebra
using GAP

function shortestVectors(gramMatrix::Array{Int64,2}, limitNormSquare::Int64 = -1)::Array{Array{Int64,1},1}
	@assert(size(gramMatrix,1) == size(gramMatrix,2))
	if limitNormSquare < 0
		limitNormSquare = maximum(Diagonal(gramMatrix))
	end
	gramMatrixGAP = GAP.julia_to_gap(map(GAP.julia_to_gap, gramMatrix))
	recGAP = GAP.Globals.ShortestVectors(gramMatrixGAP, limitNormSquare)
	vectors = GAP.gap_to_julia(Array{Array{Int64,1},1}, recGAP.vectors)
	return vcat(map( v -> [v,-v], vectors)...)
end

function scalarProduct(gramMatrix::Array{Int64,2}, vector1::Array{Int64,1}, vector2::Array{Int64,1} = vector1)::Int64
	@assert(size(gramMatrix,1) == length(vector1))
	@assert(size(gramMatrix,2) == length(vector2))
	return transpose(vector1) * gramMatrix * vector2
end

function hasVectorAccordingScalarProducts(gramMatrix::Array{Int64,2}, vector::Array{Int64,1}, vectorComparisonList::Array{Array{Int64,1},1}, scalarProductsComparisonList::Array{Int64,1})::Bool
	n = length(vectorComparisonList)
	@assert(size(gramMatrix,1) == size(gramMatrix,2))
	@assert(size(gramMatrix,1) == length(vector))
	for v in vectorComparisonList
		@assert(size(gramMatrix,1) == length(v))
	end
	@assert(length(scalarProductsComparisonList) == n)
	for i in 1:n
		if scalarProduct(gramMatrix, vector, vectorComparisonList[i]) != scalarProductsComparisonList[i]
			return false
		end
	end
	return true
end

function pleskenSouvignier_two_matrices(gramMatrix1::Array{Int64,2},gramMatrix2::Array{Int64,2}, shortVectors::Array{Array{Int64,1},1}, imagesOfBasis::Array{Array{Int64,1},1} = convert(Array{Array{Int64,1},1}, []))::Array{Array{Int64,2},1}
    n = length(imagesOfBasis) + 1
	if n > size(gramMatrix1,1)
		return [ reshape(hcat(imagesOfBasis...), (length(imagesOfBasis[1]), length(imagesOfBasis))) ]
	end
	imageCandidates = shortVectors
	# determine all vectors of the correct norm
	imageCandidates = filter( v -> scalarProduct(gramMatrix1, v) == gramMatrix2[n,n], imageCandidates)
	# determine all vectors of correct scalar products with given images
	imageCandidates = filter( v -> hasVectorAccordingScalarProducts(gramMatrix1, v, imagesOfBasis, gramMatrix2[n,1:(n-1)]), imageCandidates)
	res = []
	for possibleImage in imageCandidates
		nextImagesOfBasis = vcat(imagesOfBasis, [possibleImage])
		res = vcat(res, pleskenSouvignier_two_matrices(gramMatrix1, gramMatrix2, shortVectors, nextImagesOfBasis))
	end
	return res
end

