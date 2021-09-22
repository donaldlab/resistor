module MutSigProbs

using Printf
using Logging
using FASTX
using BioSequences
using LightGraphs, SimpleWeightedGraphs
using ArgParse

include("translation-tables.jl")
using .TranslationTables

export calculateMutationalProbabilities

FiveMerProb = @NamedTuple{fiveMer::LongDNASeq, prob::Float64}

function groupByAAType(d)
	outD = Dict{AminoAcid, Float64}()

	for (key, value) in d
		aa = CodonTable[DNACodon(key[2:4])]
		if haskey(outD, aa)
			outD[aa] += value
		else
			outD[aa] = value
		end
	end

	return outD
end

notBase(base) = setdiff(['A', 'C', 'G', 'T'], [convert(Char, base)])

"Returns a sequence of FiveMerProb for all five-mers that are reachable through mutating one of the middle three bases" 
function reachableInOneStep(probs, fiveMer)
	s = convert(String, fiveMer)
	first_bases = [FiveMerProb((LongDNASeq([s[1], j, s[3], s[4], s[5]]), probs[s[1:3]][String([s[1], j, s[3]])])) for j=notBase(fiveMer[2])]
	second_bases = [FiveMerProb((LongDNASeq([s[1], s[2], j, s[4], s[5]]), probs[s[2:4]][String([s[2], j, s[4]])])) for j=notBase(fiveMer[3])]
	third_bases = [FiveMerProb((LongDNASeq([s[1], s[2], s[3], j, s[5]]), probs[s[3:5]][String([s[3], j, s[5]])])) for j=notBase(fiveMer[4])]
	union(first_bases, second_bases, third_bases)
end

function descendents(graph, indexDict, fiveMers, fiveMer, depth)

	function descendents_i(graph, prob, fiveMer, depth, hist, probs)

		if depth < 1
			return []
		end

		idx = indexDict[fiveMer]
		state = []
		for neighbor in LightGraphs.neighbors(graph, idx)
			push!(state, (
						  fiveMers[neighbor], prob * graph.weights[neighbor, idx], 
						  join([convert(String, fm) for fm in vcat(hist, [fiveMers[neighbor]])], "->"),
						  "(" * join(vcat(probs, graph.weights[neighbor, idx]), " * ") * ")",
						  )
				  )
		end

		for neighbor in LightGraphs.neighbors(graph, idx)
			state = vcat(state, descendents_i(graph, prob * graph.weights[neighbor, idx], fiveMers[neighbor], depth - 1, vcat(hist, [fiveMers[neighbor]]), vcat(probs, graph.weights[neighbor, idx])))
		end

		return state
	end

	return descendents_i(graph, 1.0, fiveMer, depth, [fiveMer], [1.0])
end

"""
    calculateMutationalProbabilities(probs, sequence)

Given a probabilities dict and a sequence, return two dicts, first of probability of
AA to mutate to another AA, second of the actual fivemer at the location of the amino acid.

Dicts are of type (int, AminoAcid, AminoAcid) => float|string
where the int is the position of the AminoAcid, the first AminoAcid is the identity of the AminoAcid,
and the second AminoAcid is the mutation target AminoAcid.
"""
function calculateMutationalProbabilities(probs, sequence)
	bases = [DNA_A, DNA_C, DNA_G, DNA_T]
	fiveMers = [LongDNASeq([i, j, k, l, m]) for i=bases, j=bases, k=bases, l=bases, m=bases]
	d = Dict{LongDNASeq, Integer}()
	graph = SimpleWeightedDiGraph(4^5)

	# Setup a dictionary mapping fivemers to their index in the array,
	# The index is how nodes in the graph are accessed.
	for (index, value) in enumerate(fiveMers)
		d[value] = index
	end

	# Build the graph
	for (index, fivemer) in enumerate(fiveMers)
		for next_fivemer in reachableInOneStep(probs, fivemer)
			add_edge!(graph, index, d[next_fivemer.fiveMer], next_fivemer.prob)
			if next_fivemer.prob != graph.weights[d[next_fivemer.fiveMer], index]
				@error "ERROR!: ($(reachable.prob) != $(graph.weights[d[reachable.fiveMer], index]))"
			end
		end
	end

	allProbsDict = Dict()
	fiveMersDict = Dict()
	startResidue = 2

	# For each codon (exempting first and last codons)
	for i in range(4; length=(length(sequence) ÷ 3 - 2), step=3)
		codon = DNACodon(sequence[i:i+2])
		fivemer = sequence[i-1 : i+3]

		# Find all codons that can be reached with two mutations
		conversionTargets = descendents(graph, d, fiveMers, fivemer, 2)
		for (key, value, history, probs) in conversionTargets
			@debug @sprintf "%s\t%.6e\t%s\t\t%s\n" convert(String, key) value history probs
		end
		
		# Group all reached codons and probabilities by mutational amino acid type
		targetAAs = groupByAAType(conversionTargets)

		for aa in Set(values(CodonTable))
			prob = get(targetAAs, aa, 0)
			allProbsDict[(startResidue, CodonTable[codon], aa)] = prob
			fiveMersDict[(startResidue, CodonTable[codon], aa)] = fivemer
		end
		startResidue += 1
	end

	return (allProbsDict, fiveMersDict)
end

end # end module
