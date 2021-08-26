using FASTX
using JSON
using BioSequences
using LightGraphs, SimpleWeightedGraphs
using Printf
using XLSX

function parseJsonFile(f)
	open(f, "r") do reader
		JSON.parse(reader)
	end
end

mutProbFile = "./melanoma.json"
startResidue = 2
bases = [DNA_A, DNA_C, DNA_G, DNA_T]
fiveMers = [LongDNASeq([i, j, k, l, m]) for i=bases, j=bases, k=bases, l=bases, m=bases]
d = Dict{LongDNASeq, Integer}()
probs = parseJsonFile(mutProbFile)
graph = SimpleWeightedDiGraph(4^5)
FiveMerProb = @NamedTuple{fiveMer::LongDNASeq, prob::Float64}

CodonTable = Dict{DNACodon, AminoAcid}(
	DNACodon("ATT") => AA_I, 
	DNACodon("ATC") => AA_I, 
	DNACodon("ATA") => AA_I,
	DNACodon("CTT") => AA_L,
	DNACodon("CTC") => AA_L,
	DNACodon("CTA") => AA_L,
	DNACodon("CTG") => AA_L,
	DNACodon("TTA") => AA_L,
	DNACodon("TTG") => AA_L,
	DNACodon("GTT") => AA_V,
	DNACodon("GTC") => AA_V,
	DNACodon("GTA") => AA_V,
	DNACodon("GTG") => AA_V,
	DNACodon("TTT") => AA_F,
	DNACodon("TTC") => AA_F,
	DNACodon("ATG") => AA_M,
	DNACodon("TGT") => AA_C,
	DNACodon("TGC") => AA_C,
	DNACodon("GCT") => AA_A,
	DNACodon("GCC") => AA_A,
	DNACodon("GCA") => AA_A,
	DNACodon("GCG") => AA_A,
	DNACodon("GGT") => AA_G,
	DNACodon("GGC") => AA_G, 
	DNACodon("GGA") => AA_G,
	DNACodon("GGG") => AA_G,
	DNACodon("CCT") => AA_P,
	DNACodon("CCC") => AA_P,
	DNACodon("CCA") => AA_P,
	DNACodon("CCG") => AA_P,
	DNACodon("ACT") => AA_T,
	DNACodon("ACC") => AA_T,
	DNACodon("ACA") => AA_T,
	DNACodon("ACG") => AA_T,
	DNACodon("TCT") => AA_S,
	DNACodon("TCC") => AA_S,
	DNACodon("TCA") => AA_S, 
	DNACodon("TCG") => AA_S,
	DNACodon("AGT") => AA_S,
	DNACodon("AGC") => AA_S,
	DNACodon("TAT") => AA_Y,
	DNACodon("TAC") => AA_Y,
	DNACodon("TGG") => AA_W,
	DNACodon("CAA") => AA_Q,
	DNACodon("CAG") => AA_Q,
	DNACodon("AAT") => AA_N,
	DNACodon("AAC") => AA_N,
	DNACodon("CAT") => AA_H,
	DNACodon("CAC") => AA_H,
	DNACodon("GAA") => AA_E,
	DNACodon("GAG") => AA_E,
	DNACodon("GAT") => AA_D,
	DNACodon("GAC") => AA_D,
	DNACodon("AAA") => AA_K,
	DNACodon("AAG") => AA_K,
	DNACodon("CGT") => AA_R,
	DNACodon("CGC") => AA_R,
	DNACodon("CGA") => AA_R,
	DNACodon("CGG") => AA_R,
	DNACodon("AGA") => AA_R,
	DNACodon("AGG") => AA_R,
	DNACodon("TAA") => AA_Term,
	DNACodon("TAG") => AA_Term,
	DNACodon("TGA") => AA_Term
)

AA1_3 = Dict{Char, String}(
	'A' => "Ala",
	'R' => "Arg",
    'N' => "Asn",
    'D' => "Asp",
    'C' => "Cys",
    'Q' => "Gln",
    'E' => "Glu",
    'G' => "Gly",
    'H' => "His",
    'I' => "Ile",
    'L' => "Leu",
    'K' => "Lys",
    'M' => "Met",
    'F' => "Phe",
    'P' => "Pro",
    'S' => "Ser",
    'T' => "Thr",
    'W' => "Trp",
    'Y' => "Tyr",
    'V' => "Val"
)

AA3_1 = Dict{String, Char}(
	"Ala" => 'A',
	"Arg" => 'R',
    "Asn" => 'N',
    "Asp" => 'D',
    "Cys" => 'C',
    "Gln" => 'Q',
    "Glu" => 'E',
    "Gly" => 'G',
    "His" => 'H',
    "Hip" => 'H',
    "Hid" => 'H',
    "Hie" => 'H',
    "Ile" => 'I',
    "Leu" => 'L',
    "Lys" => 'K',
    "Met" => 'M',
    "Phe" => 'F',
    "Pro" => 'P',
    "Ser" => 'S',
    "Thr" => 'T',
    "Trp" => 'W',
    "Tyr" => 'Y',
    "Val" => 'V'
)

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
return

fastaIn = getRecordSequence("BRAF-cDNA.fasta", "BRAF")
# fastaIn = getRecordSequence("BRAF-cDNA-test.fasta", "BRAF")
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
