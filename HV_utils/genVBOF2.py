from cobra import Reaction, Metabolite
from copy import deepcopy
from collections import Counter
from HV_utils.info import ntpsDict, metDict, aaDict, N_A, k_atp_protein, k_atp_rna, k_ppi, proteins_per_mrna

#######################################################################################
# Function Setup

# NOTE: the keys in this dictionary correspond to the name of the product in the genbank file
virus_composition = {
    "HpV18": {              # HPV-18 (locus)                
                "Cg": 1,            # Copy number for viral genome
                "proteins": {
                            # Copy numbers for structural proteins L1 and L2 for HPV
                            # https://ictv.global/report/chapter/papillomaviridae/papillomaviridae
                            # HPV is a non-enveloped virus, so the capsid is made of L1 and L2 proteins. 
                            # The capsid is composed of 360 copies of L1 and 12 copies of L2.
                            "L1": 360, 
                            "L2": 12,}
            }
}

def multiply_counter(counter, n):
    """return a collection.Counter whose values have been multiplied by a
    scalar"""
    # the multiplicand does not need to be a collections.Counter, but at least
    # it needs to be a dict
    assert isinstance(counter, dict)
    assert float(n)
    return Counter({key: value * n for key, value in counter.items()})

def get_virus_names(virus_record):
    """Return a tuple containing a short name and a
    full name for the virus. This is an alternative to how it is done in Sean's
    genVBOF."""
    long_name = virus_record.description.rstrip(", complete genome")
    # try to get a short name for the virus by taking the prefix of the first CDS
    short_name = next(feature for feature in virus_record.features if feature.type == "CDS").qualifiers["locus_tag"][0].rstrip("gp1")
    return (short_name, long_name)

model_name_to_met_dict_file = {
    "Swainston2016 - Reconstruction of human metabolic network (Recon 2.2)": "met_dicts/recon2_met_dict.txt",
    "Genome-scale metabolic model for hepatocytes, iHepatocytes2322": "met_dicts/iHepatocytes2322_met_dict.txt",
    "macrophage_model": "met_dicts/mac_met_dict.txt",
    "recon 2_LarsNielsen model": "met_dicts/recon2_LarsNielsen_met_dict.txt",
    "GSM_human model": "met_dicts/GSM_human_met_dict.txt", "lung model": "met_dicts/lung_met_dict.txt"
}

# def load_metabolite_id_dict(model, model_name=None):
#     """Provide a dictionary mapping generic metabolic name (a string like "A",
#     "atp", "h2o", "h" or "PPi") to the metabolite object in the model.
    
#     The model_name parameter allows to indicate the name of the model in case
#     the model object has no name attribute. Do not set this parameter if the
#     model object already has a name."""

#     name = model_name or model.name
#     if name in model_name_to_met_dict_file:
#         met_dict_file = model_name_to_met_dict_file[name]
#     else:
#         raise NotImplementedError("This model is not covered: \"{}\"".format(name))

#     with open(met_dict_file, "r") as fh:
#         met_tuples = (line.split(", ")[:2] for line in fh.readlines()[2:])
#     met_dict = {key: model.metabolites.get_by_id(met_id) for key, met_id in met_tuples}
#     return met_dict

def reverse_complement(seq):
    rc = ""
    for nuc in seq[::-1]:
        if nuc == "A":
            rc += "T"
        elif nuc == "T":
            rc += "A"
        elif nuc == "C":
            rc += "G"
        elif nuc == "G":
            rc += "C"
        else:
            raise ValueError("Invalid nucleotide: %s" % nuc)
    return rc

#######################################################################################
# Function Definition
# genVOBF.py takes user-supplied viral genome file (.gb) and creates a
# biomass objective function that is characterises the genomic, proteomic and
# energy requirements for the production of virus particles

# Inputs:
# VirusGB           User-supplied GenBank file (NCBI) for desired virus

# Outputs:
# VBOF              Virus biomass objective function for desired virus

def genVBOF2(virus_record, model, model_name=None):
    """New version of the genVBOF function by Hadrien.
    Builds a Virus Biomass Objective Function (basically a virus biomass
    production reaction, from aminoacids and nucleotides) from a genbank
    file.
    
    Params:
    - virus_record: genbank record of a virus (output from Bio.SeqIO.parse)
    - model: a cobra metabolic model (cobra.core.model.Model)

    Returns:
    - virus biomass objective function (cobra.core.reaction.Reaction)
    """

    # VIRUS IDENTIFICATION
    taxonomy = " ".join([taxon.lower() for taxon in virus_record.annotations["taxonomy"]])
    if "papillomaviridae" not in taxonomy:
        raise NotImplementedError('Virus family is not supported: Unable to create VBOF. Consult _README')
    short_name, full_name = get_virus_names(virus_record)

    # AMINOACID COUNT
    all_cds = [feature for feature in virus_record.features if feature.type == "CDS"]
    # Check that our own virus_composition dict contain exactly the
    # proteins defined in the genbank file, no more, no less.
    protein_names_in_gb_file = {cds.qualifiers["gene"][0] for cds in all_cds}
    protein_names_in_our_data = {protein_name for protein_name in virus_composition[short_name]["proteins"]}
    assert protein_names_in_our_data.issubset(protein_names_in_gb_file)

    virus_aa_composition = Counter()
    virus_mrna_composition = Counter()
    # protein name -> number of atp involved in its peptide bonds formations
    # (accounting for the number of copies of protein)
    peptide_bond_formation = dict()
    for cds in all_cds:
        if cds.qualifiers["gene"][0] in protein_names_in_our_data:
            protein_name = cds.qualifiers["gene"][0]
            aa_sequence = cds.qualifiers["translation"][0]
            dna_sequence = virus_record.seq[cds.location.start : cds.location.end]
            if cds.location.strand == -1:
                dna_sequence = reverse_complement(dna_sequence)
            aa_count = Counter(aa_sequence)
            nc_count = Counter(dna_sequence)

            copies_per_virus = virus_composition[short_name]["proteins"][protein_name]

            virus_aa_composition += multiply_counter(aa_count, copies_per_virus)
            virus_mrna_composition += multiply_counter(nc_count, copies_per_virus)

            peptide_bond_formation[protein_name] = (
                len(aa_sequence) * k_atp_protein - k_atp_protein
                ) * copies_per_virus

    # [3] Precursor frequency
    # Genome                            [Nucleotides]
    Cg = virus_composition[short_name]["Cg"] # number of genome copies per virus
    virus_nucl_count = Counter(str(virus_record.seq))
    countA  = virus_nucl_count["A"]
    countC  = virus_nucl_count["C"]
    countG  = virus_nucl_count["G"]
    countT  = virus_nucl_count["T"]    # Base 'T' is pseudo for base 'U'
    antiA   = countT
    antiC   = countG
    antiG   = countC
    antiT   = countA

    # Note: mRNA is 1-stranded, so we don't track the "anti" parts
    # Also: we didn't bother to replace T->U in the string; we do it here
    countA_rna = virus_mrna_composition["A"] / proteins_per_mrna
    countC_rna = virus_mrna_composition["C"] / proteins_per_mrna
    countG_rna = virus_mrna_composition["G"] / proteins_per_mrna
    countU_rna = virus_mrna_composition["T"] / proteins_per_mrna

    # Count summation
    totNTPS     = (Cg * (countA + countC + countG + countT + antiA + antiC + antiG + antiT))
    totAA       = sum(count for count in virus_aa_composition.values())

    # [4] VBOF Calculations
    # Nucleotides
    # mol.dntps/mol.virus
    V_a = (Cg*(countA + antiA))
    V_c = (Cg*(countC + antiC))
    V_g = (Cg*(countG + antiG))
    V_t = (Cg*(countT + antiT))
    # g.dmps/mol.virus
    G_a = V_a * ntpsDict["damp"]
    G_c = V_c * ntpsDict["dcmp"]
    G_g = V_g * ntpsDict["dgmp"]
    G_t = V_t * ntpsDict["dtmp"]

    # mol.ntps (mrna) / mol.virus
    V_a_rna = countA_rna
    V_c_rna = countC_rna
    V_g_rna = countG_rna
    V_u_rna = countU_rna

    # Amino Acids
    # g.a/mol.virus
    G_aa = {aa: count * aaDict[aa] for aa, count in virus_aa_composition.items()}
    # Total genomic and proteomic molar mass
    M_v = (G_a + G_c + G_g + G_t) + sum(G_aa.values())

    # Stoichiometric coefficients
    # Nucleotides [mmol.ntps/g.virus] (for the genome)
    S_datp = 1000 * (V_a / M_v)
    S_dctp = 1000 * (V_c / M_v)
    S_dgtp = 1000 * (V_g / M_v)
    S_dttp = 1000 * (V_t / M_v)

    # mRNA is handeled differently because the virus doesn't pack it. We will
    # allow it to decay into NMPs at a rate governed by the protein to mrna
    # ratio (for mass conservation)
    S_atp_rna = 1000 * (V_a_rna / M_v)
    S_ctp_rna = 1000 * (V_c_rna / M_v)
    S_gtp_rna = 1000 * (V_g_rna / M_v)
    S_utp_rna = 1000 * (V_u_rna / M_v)
    S_amp_rna = 1000 * (V_a_rna / M_v)
    S_cmp_rna = 1000 * (V_c_rna / M_v)
    S_gmp_rna = 1000 * (V_g_rna / M_v)
    S_ump_rna = 1000 * (V_u_rna / M_v)

    # Amino acids [mmol.aa/g.virus]
    S_aa = {aa: 1000 * V_aa / M_v for aa, V_aa in virus_aa_composition.items()}

    # Energy requirements
    n_proteins = len(virus_composition[short_name]["proteins"])
    S_ppi = (
        S_atp_rna
        + S_ctp_rna
        + S_gtp_rna
        + S_utp_rna
        + S_datp
        + S_dctp
        + S_dgtp
        + S_dttp
        - (2 + n_proteins / proteins_per_mrna) * 1000 / M_v
    ) * k_ppi

    # Proteome: Peptide bond formation [ATP + H2O]
    # Note: ATP used in this process is denoated as ATPe/Ae [e = energy version]
    V_Ae = sum(peptide_bond_formation.values())
    S_Ae = 1000 * (V_Ae / M_v)
 
    # [5] VBOF Reaction formatting and output
    # Left-hand terms: Nucleotides
    S_ATP = S_Ae * -1 - S_atp_rna
    S_CTP = S_ctp_rna * -1
    S_GTP = S_gtp_rna * -1
    S_UTP = S_utp_rna * -1
    S_AMP = S_amp_rna
    S_CMP = S_cmp_rna
    S_GMP = S_gmp_rna
    S_UMP = S_ump_rna
    S_DATP = S_datp * -1
    S_DCTP = S_dctp * -1
    S_DGTP = S_dgtp * -1
    S_DTTP = S_dttp * -1
 
    # Left-hand terms: Amino Acids
    S_AAf = {aa: -coef for aa, coef in S_aa.items()}
    # Left-hand terms: Energy Requirements
    S_H2O   = S_Ae * -1
    # Right-hand terms: Energy Requirements
    S_ADP   = S_Ae
    S_Pi    = S_Ae
    S_H     = S_Ae
    S_PPi   = S_ppi

    reaction_name       = short_name + '_prodrxn_VN'
    virus_reaction      = Reaction(reaction_name)
    virus_reaction.name = full_name + ' production reaction'
    virus_reaction.subsystem                = 'Virus Production'
    virus_reaction.lower_bound              = 0
    virus_reaction.upper_bound              = 1000
    model.add_reactions([virus_reaction])
    virus_reaction.add_metabolites(({
        metDict['atp']: S_ATP,
        metDict['ctp']: S_CTP,
        metDict['gtp']: S_GTP,
        metDict['utp']: S_UTP,
        metDict["amp"]: S_AMP,
        metDict["cmp"]: S_CMP,
        metDict["gmp"]: S_GMP,
        metDict["ump"]: S_UMP,
        metDict["datp"]: S_DATP,
        metDict["dctp"]: S_DCTP,
        metDict["dgtp"]: S_DGTP,
        metDict["dttp"]: S_DTTP,
        metDict['A']: S_AAf['A'],
        metDict['R']: S_AAf['R'],
        metDict['N']: S_AAf['N'],
        metDict['D']: S_AAf['D'],
        metDict['C']: S_AAf['C'],
        metDict['Q']: S_AAf['Q'],
        metDict['E']: S_AAf['E'],
        metDict['G']: S_AAf['G'],
        metDict['H']: S_AAf['H'],
        metDict['I']: S_AAf['I'],
        metDict['L']: S_AAf['L'],
        metDict['K']: S_AAf['K'],
        metDict['M']: S_AAf['M'],
        metDict['F']: S_AAf['F'],
        metDict['P']: S_AAf['P'],
        metDict['S']: S_AAf['S'],
        metDict['T']: S_AAf['T'],
        metDict['W']: S_AAf['W'],
        metDict['Y']: S_AAf['Y'],
        metDict['V']: S_AAf['V'],
        metDict['h2o']: S_H2O,
        metDict['adp']: S_ADP,
        metDict['Pi']:  S_Pi,
        metDict['h']:   S_H,
        metDict['PPi']: S_PPi}))
    return virus_reaction

# Function Definition
# genHVM.py takes a user-supplied model file and creates an intergrated host-
# virus model, given a VBOF

# Inputs:
# Model             User-supplied model (cobra.io.core.model.Model instance)
# VBOF              Virus biomass objective function created by genVBOF.py

# Outputs:
# hvm               Integrated host-virus model

def genHVM(Model,VBOF):
    "Generate_HVM"

    # Integrate the VBOF Reaction
    HVM = deepcopy(Model)
    HVM.add_reaction(VBOF)

    # Outputs
    return HVM

