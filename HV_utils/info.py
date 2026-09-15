#######################################################################################
# DEFINITIONS

# Metabolite (verbose to model) definitions
metDict = {
    'atp': 'atp_c',        # ATP,              ChEBI 15422
    'ctp': 'ctp_c',        # CTP,              ChEBI 17677
    'gtp': 'gtp_c',        # GTP,              ChEBI 15996
    'utp': 'utp_c',        # UTP,              ChEBI 15713
    'amp': 'amp_c',        
    'cmp': 'cmp_c',        
    'gmp': 'gmp_c',       
    'ump': 'ump_c',        
    'datp': 'datp_c',      
    'dctp': 'dctp_c',      
    'dgtp': 'dgtp_c',      
    'dttp': 'dttp_c',      
    'A': 'ala__L_c',        # Alaline,          ChEBI 16977
    'R': 'arg__L_c',        # Arginine,         ChEBI 16467
    'N': 'asn__L_c',        # Asparagine,       ChEBI 17196
    'D': 'asp__L_c',        # Aspartate,        ChEBI 17053
    'C': 'cys__L_c',        # Cysteine,         ChEBI 17561
    'Q': 'gln__L_c',        # Glutamine,        ChEBI 18050
    'E': 'glu__L_c',        # Glutamate,        ChEBI 16015
    'G': 'gly_c',          # Glycine,          ChEBI 15428
    'H': 'his__L_c',        # Histidine,        ChEBI 15971
    'I': 'ile__L_c',        # Isoleucine,       ChEBI 17191
    'L': 'leu__L_c',        # Leucine,          ChEBI 15603
    'K': 'lys__L_c',        # Lysine,           ChEBI 18019
    'M': 'met__L_c',        # Methionine,       ChEBI 16643
    'F': 'phe__L_c',        # Phenylalanine,    ChEBI 17295
    'P': 'pro__L_c',        # Proline,          ChEBI 17203
    'S': 'ser__L_c',        # Serine,           ChEBI 17115
    'T': 'thr__L_c',        # Threonine,        ChEBI 16857
    'W': 'trp__L_c',        # Tryptophan,       ChEBI 16828
    'Y': 'tyr__L_c',        # Tyrosine,         ChEBI 17895
    'V': 'val__L_c',        # Valine,           ChEBI 16414
    'h2o': 'h2o_c',        # H2O
    'adp': 'adp_c',        # ADP
    'Pi': 'pi_c',          # Phosphate
    'h': 'h_c',            # Hydrogen [Proton]
    'PPi': 'ppi_c',        # Pyrophosphate
}

# Nucleotides Dictionary
ntpsDict = {
    "atp": 507.181,
    "gtp": 483.15644,
    "ctp": 523.18062,
    "utp": 484.14116,
    "datp": 491.1816,
    "dgtp": 467.1569,
    "dctp": 482.1683,
    "dttp": 507.181,
    "damp": 331.068,
    "dgmp": 347.063,
    "dcmp": 307.057,
    "dtmp": 322.057,
}
# Amino Acids Dictionary
aaDict = {
    'A': 89.09322,      # Alaline,          ChEBI 16977
    'R': 174.201,       # Arginine,         ChEBI 16467
    'N': 132.118,       # Asparagine,       ChEBI 17196
    'D': 133.1027,      # Aspartate,        ChEBI 17053
    'C': 121.158,       # Cysteine,         ChEBI 17561
    'Q': 146.14458,     # Glutamine,        ChEBI 18050
    'E': 147.1293,      # Glutamate,        ChEBI 16015
    'G': 75.06664,      # Glycine,          ChEBI 15428
    'H': 155.15468,     # Histidine,        ChEBI 15971
    'I': 131.17296,     # Isoleucine,       ChEBI 17191
    'L': 131.17296,     # Leucine,          ChEBI 15603
    'K': 146.18764,     # Lysine,           ChEBI 18019
    'M': 149.21238,     # Methionine,       ChEBI 16643
    'F': 165.18918,     # Phenylalanine,    ChEBI 17295
    'P': 115.1305,      # Proline,          ChEBI 17203
    'S': 105.09262,     # Serine,           ChEBI 17115
    'T': 119.1192,      # Threonine,        ChEBI 16857
    'W': 204.22526,     # Tryptophan,       ChEBI 16828
    'Y': 181.18858,     # Tyrosine,         ChEBI 17895
    'V': 117.14638,     # Valine,           ChEBI 16414
}
# Misc. Dictionary
miscDict = {
    'PPi': 173.94332,   # Pyrophosphate,    ChEBI 18361
}
# Avogadro's Number
N_A         = 6.0221409e+23
# ATP requirement coefficients
# source: Queka, Dietmaira, Hanschob, Mart�neza, Borthb, Nielsen
# Journal of Biotechnology 184 (2014) 172�178
k_atp_protein = 4.3
k_atp_dna = 1.4
k_atp_rna = 0.4

k_ppi = 1

proteins_per_mrna = 5616 # median taken from https://doi.org/10.1371/journal.pone.0073943


