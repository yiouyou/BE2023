from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.teaching import *

#Core Includes
from pyrosetta.rosetta.core.select import residue_selector as selections
from pyrosetta.rosetta.protocols import antibody

_pdb = "/home/songz/PyRosetta.notebooks/notebooks/inputs/"

# Intitlialization
init('-use_input_sc -ignore_unrecognized_res \
     -ignore_zero_occupancy false -load_PDB_components false -no_fconfig')

# Import and copy pose
pose = pose_from_pdb(_pdb+"2r0l_1_1.pdb")
original_pose = pose.clone()

# AntibodyInfo
ab_info = antibody.AntibodyInfo(pose, antibody.AHO_Scheme, antibody.North)
print(ab_info)

# Basic AntibodyInfo Access
print("h1", ab_info.get_CDR_start(antibody.h1, pose))
print("h2", ab_info.get_CDR_end(antibody.h2, pose))
for i in range(1, 7):
    print(i, ab_info.get_CDR_name(antibody.CDRNameEnum(i)))
for cdr in ['L1', 'l1', 'L2', 'l2', 'L3', 'H1', 'H2', 'H3']:
    print(cdr, str(ab_info.get_CDR_name_enum(cdr)))
print(str(antibody.h3))
print(int(antibody.h3))

# AntibodyEnumManager
enum_manager = antibody.AntibodyEnumManager()
print(enum_manager.numbering_scheme_enum_to_string(antibody.AHO_Scheme))
print(enum_manager.cdr_definition_enum_to_string(antibody.North))
print(enum_manager.cdr_name_string_to_enum("H1"))
print(enum_manager.antibody_region_enum_to_string(antibody.framework_region))
for i in range(1, pose.size()+1):
    region = ab_info.get_region_of_residue(pose, i)
    if (region == antibody.cdr_region):
        print(i, enum_manager.cdr_name_enum_to_string(ab_info.get_CDRNameEnum_of_residue(pose, i)))
    else:
        print(i, enum_manager.antibody_region_enum_to_string(region))

# CDR Clusters
print(ab_info.get_CDR_length(antibody.l1))
print(ab_info.get_CDR_cluster(antibody.l1).cluster())
L1_cluster = ab_info.get_CDR_cluster(antibody.l1)
print(L1_cluster.normalized_distance_in_degrees())

# Conserved Inter-Domain Cysteine
rosetta_num = ab_info.get_landmark_resnum(pose, antibody.Kabat_Scheme, 'H', 22)
print(pose.pdb_info().pose2pdb(rosetta_num))
print(pose.residue(rosetta_num))
pre_cdr3_c = ab_info.get_landmark_resnum(pose, antibody.IMGT_Scheme, 'H', 104)
print(pose.pdb_info().pose2pdb(pre_cdr3_c))
print(pose.residue(pre_cdr3_c))

# Sequence
ab_seq = ab_info.get_antibody_sequence()
print(ab_seq)
L1_seq = ab_info.get_CDR_sequence_with_stem(antibody.l1, pose)
print("L1", L1_seq)
for i in range(1, 7):
    cdr = antibody.CDRNameEnum(i)
    print(cdr, ab_info.get_CDR_sequence_with_stem(cdr, pose))

# Function: get_cdr_loops()
h3_l3 = rosetta.utility.vector1_bool(6)
print(h3_l3)
h3_l3[antibody.h3] = True
h3_l3[antibody.l3] = True
h3_l3_loops = antibody.get_cdr_loops(ab_info, pose, h3_l3, 2)
print(h3_l3_loops)

# Function: select_epitope_residues()
epi_residues = antibody.select_epitope_residues(ab_info, pose, 8)
total=0
for i in range(1, len(epi_residues)+1):
    if epi_residues[i]:
        print(i)
        total+=1
print("Total Epitope Residues:", total)

# SasaMetric, TotalEnergyMetric, SelectedResiduesPyMOLMetric
epi_res_selector = selections.ReturnResidueSubsetSelector(epi_residues)
import pyrosetta.rosetta.core.simple_metrics.metrics as sm
sasa_metric = sm.SasaMetric(epi_res_selector)
print("\nSASA", sasa_metric.calculate(pose))
total_metric = sm.TotalEnergyMetric(epi_res_selector)
print("\nTOTAL RESIDUE ENERGY", total_metric.calculate(pose))
pymol_metric = sm.SelectedResiduesPyMOLMetric(epi_res_selector)
print("\nSELECTION", pymol_metric.calculate(pose))

# PerResidueSasaMetric
import pyrosetta.rosetta.core.simple_metrics.per_residue_metrics as residue_sm
import operator
res_sasa_metric = residue_sm.PerResidueSasaMetric()
res_sasa_metric.set_residue_selector(epi_res_selector)
per_res_sasa = res_sasa_metric.calculate(pose)
#print(per_res_sasa) 
for ele in sorted(per_res_sasa.items(), key=operator.itemgetter(1), reverse=False):
    print(ele)

# PerResidueEnergyMetric
res_energy_metric = residue_sm.PerResidueEnergyMetric()
res_energy_metric.set_residue_selector(epi_res_selector)
per_res_energy = res_sasa_metric.calculate(pose)
#print(per_res_sasa)
for ele in sorted(per_res_energy.items(), key=operator.itemgetter(1), reverse=False):
    print(ele[0], pose.pdb_info().pose2pdb(ele[0]), ele[1])

# CDRResidueSelector
from pyrosetta.rosetta.protocols.antibody.residue_selector import *
cdr_selector = CDRResidueSelector(ab_info)
cdr_selector.set_cdrs(h3_l3)
sele = cdr_selector.apply(pose)
for i in range(1, len(sele)):
    if sele[i]:
        print(i, pose.pdb_info().pose2pdb(i))

# AntibodyRegionSelector
region_selector = AntibodyRegionSelector(ab_info)
region_selector.set_region(antibody.cdr_region)
sele = region_selector.apply(pose)
for i in range(1, len(sele)):
    if sele[i]:
        print(i, pose.pdb_info().pose2pdb(i))


