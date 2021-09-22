using Test
using BioSequences
import JSON

include("mut-sig-probs.jl")
using .MutSigProbs

melanoma_probs = """
{
    "GGA" : {
       "GAA" : 0.218572395,
       "GCA" : 0.000958632,
       "GTA" : 0.002185331
    },
    "GAT" : {
       "GCT" : 5.0796e-05,
       "GGT" : 0.000731155,
       "GTT" : 0.001277772
    },
    "GGG" : {
       "GAG" : 0.077262288,
       "GCG" : 0.000476074,
       "GTG" : 0.000736481
    },
    "AAT" : {
       "ACT" : 0.000203184,
       "AGT" : 0.002559217,
       "ATT" : 0.003134217
    }
}
""" |> JSON.parse

function calculate_g466e_mutational_probability()
	# G466E with melanoma cancer probabilities
	# Codon TGGAT
	#
	# 1-step: 
	#	GGA ->(0.218572395) GAA
	# 2-step:
	#	GGA ->(0.000731155) GGG ->(0.077262288) GAG = 5.649070818264e-5
	#	GGA ->(0.218572395) GAA ->(0.002559217) GAG = 0.000559374189014715
	init(melanoma_probs)
	return calculateFivemerMutationalProbabilities(LongDNASeq("TGGAT"))[AA_E]
end

@testset "Mutational Signature Probability Calculations" begin
	@test calculate_g466e_mutational_probability() == 0.21918825989719737
end
