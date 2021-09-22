using FASTX
using JSON
using BioSequences
using LightGraphs, SimpleWeightedGraphs
using Printf
using XLSX
using ArgParse

include("translation-tables.jl")
using .TranslationTables

global probs
global graph

startResidue = 2
bases = [DNA_A, DNA_C, DNA_G, DNA_T]
fiveMers = [LongDNASeq([i, j, k, l, m]) for i=bases, j=bases, k=bases, l=bases, m=bases]
d = Dict{LongDNASeq, Integer}()
FiveMerProb = @NamedTuple{fiveMer::LongDNASeq, prob::Float64}

function parse_command_line()
	s = ArgParseSettings()
	@add_arg_table! s begin
		"--mut-prob"
			help = "the path to the JSON file with three-mer mutational probabilities"
			required = true
		"--fasta"
			help = "the path to the fasta file with the cDNA sequence"
			required = true
		"--identifier"
			help = "the identifier of the sequence in the fasta to extract"
			required = true
	end

	return parse_args(ARGS, s)
end

function parseJsonFile(f)
	open(f, "r") do reader
		JSON.parse(reader)
	end
end

function toAa(s)
	AminoAcid(AA3_1[titlecase(s)])
end

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

function getRecordSequence(f, identifier)
	open(FASTA.Reader, f) do reader
		for record in reader
			if FASTA.identifier(record) == identifier
				return FASTA.sequence(record)
			end
		end
	end
end

notBase(base) = setdiff(['A', 'C', 'G', 'T'], [convert(Char, base)])

function reachableInOneStep(fiveMer)
	s = convert(String, fiveMer)
	first_bases = [FiveMerProb((LongDNASeq([s[1], j, s[3], s[4], s[5]]), probs[s[1:3]][String([s[1], j, s[3]])])) for j=notBase(fiveMer[2])]
	second_bases = [FiveMerProb((LongDNASeq([s[1], s[2], j, s[4], s[5]]), probs[s[2:4]][String([s[2], j, s[4]])])) for j=notBase(fiveMer[3])]
	third_bases = [FiveMerProb((LongDNASeq([s[1], s[2], s[3], j, s[5]]), probs[s[3:5]][String([s[3], j, s[5]])])) for j=notBase(fiveMer[4])]
	union(first_bases, second_bases, third_bases)
end

function descendents(graph, fiveMer, depth)

	function descendents_i(graph, prob, fiveMer, depth, hist, probs)

		if depth < 1
			return []
		end

		idx = d[fiveMer]
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


function main()
	global probs
	global graph

	args = parse_command_line()

	# block 1
	probs = parseJsonFile(args["mut-prob"])
	graph = SimpleWeightedDiGraph(4^5)

	# block 2
	for tup in enumerate(fiveMers)
		d[tup[2]] = tup[1]
	end

	for tup in enumerate(fiveMers)
		i = tup[1]
		fmer = tup[2]
		for reachable in reachableInOneStep(fmer)
			add_edge!(graph, i, d[reachable.fiveMer], reachable.prob)
			if reachable.prob != graph.weights[d[reachable.fiveMer], i]
				@printf "ERROR!: (%f != %f)\n" reachable.prob graph.weights[d[reachable.fiveMer], i]
			end

		end
	end

	# block 3
	fastaIn = getRecordSequence(args["fasta"], args["identifier"])
	allProbsDict = Dict()
	fiveMersDict = Dict()

	for i in range(4; length=(length(fastaIn) ÷ 3 - 2), step=3)
		codon = DNACodon(fastaIn[i:i+2])
		fivemer = fastaIn[i-1 : i+3]
		conversionTargets = descendents(graph, fivemer, 2)
		for (key, value, history, probs) in conversionTargets
			@printf "%s\t%.6e\t%s\t\t%s\n" convert(String, key) value history probs
		end
		
		targetAAs = groupByAAType(conversionTargets)
		global startResidue

		for aa in Set(values(CodonTable))
			prob = get(targetAAs, aa, 0)
			allProbsDict[(startResidue, CodonTable[codon], aa)] = prob
			fiveMersDict[(startResidue, CodonTable[codon], aa)] = fivemer
		end
		startResidue += 1
	end

	excelFile = "/home/nsg/nas/duke/collaborations/teresa/braf-results/braf-inhibitor-resistance.xlsx"
	XLSX.openxlsx(excelFile, mode="rw") do xf
		for name in ["vemurafenib-resistance", "encorafenib-resistance", "PLX8394-resistance", "dabrafenib-resistance"]
			sheet = xf[name]

			for row in XLSX.eachrow(sheet)
				rowNum = XLSX.row_number(row)

				if rowNum == 1
					sheet[@sprintf "N%d" rowNum] = "sig_prob"
					sheet[@sprintf "O%s" rowNum] = "codon"
					continue
				end

				resCell = XLSX.getcell(row, 1)
				wtCell = XLSX.getcell(row, 2)
				mutCell = XLSX.getcell(row, 4)

				residue = XLSX.getdata(sheet, resCell)
				wt = toAa(XLSX.getdata(sheet, wtCell))
				mut = toAa(XLSX.getdata(sheet, mutCell))

				sheet[@sprintf "N%d" rowNum] = allProbsDict[(residue, wt, mut)]
				sheet[@sprintf "O%s" rowNum] = convert(String, fiveMersDict[(residue, wt, mut)])
			end
		end
	end
end

main()
