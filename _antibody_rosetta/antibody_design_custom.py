import pyrosettacolabsetup
import pyrosetta

from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.teaching import *

from pyrosetta.rosetta.utility import *
from pyrosetta.rosetta.protocols.rosetta_scripts import *
from pyrosetta.rosetta.protocols.antibody import *
from pyrosetta.rosetta.protocols.antibody.design import *
from pyrosetta.rosetta.protocols.antibody.residue_selector import *
from pyrosetta.rosetta.protocols.antibody.task_operations import *
from pyrosetta.rosetta.protocols.antibody.constraints import *
from pyrosetta.rosetta.protocols.analysis.simple_metrics import RunSimpleMetricsMover
from pyrosetta.rosetta.core.simple_metrics.metrics import *
from pyrosetta.rosetta.core.simple_metrics.per_residue_metrics import *
from pyrosetta.rosetta.core.select.movemap import *
from pyrosetta.rosetta.core.pack.task import *

from typing import *
import pandas
from pathlib import Path
import json
import re
import os


def load_json_scorefile(file_path: Path, sort_by: str="dG_separated") -> pandas.DataFrame:
        """
        Read scorefile lines as a dataframe, sorted by total_score with Nan's correctly replaced.
        """
        local_lines = open(file_path, 'r').readlines()
        decoys=[]
        for line in local_lines:
                o = json.loads(line.replace("nan", "NaN"))
                # print o[self.decoy_field_name]
                # print repr(o)
                decoys.append(o)
        local_df = pandas.DataFrame.from_dict(decoys)
        local_df = local_df.infer_objects()
        # df.to_csv("debugging.csv", sep=",")
        local_df = local_df.sort_values(sort_by, ascending=True)
        return local_df

def drop_cluster_columns(local_df: pandas.DataFrame, keep_cdrs: List[str]=None) -> pandas.DataFrame:
        """
        Drop cluster columns that RAbD outputs to make it easier to work with the dataframe.
        """
        to_drop = []
        for column in local_df.columns:
            if re.search("cdr_cluster", column):
                skip=False
                if (keep_cdrs):
                    for cdr in keep_cdrs:
                        if re.search(cdr, column):
                            skip=True
                            break
                if not skip:
                    to_drop.append(column)
        return local_df.drop(columns=to_drop)


pyrosettacolabsetup.install_pyrosetta()
pyrosetta.init()
init('-no_fconfig @/home/songz/PyRosetta.notebooks/notebooks/inputs/rabd/common -check_cdr_chainbreaks false')

pose = pose_from_pdb("/home/songz/PyRosetta.notebooks/notebooks/inputs/rabd/my_ab.pdb")
original_pose = pose.clone()

# Custom Design Protocol in code
pose = original_pose.clone()
ab_info = AntibodyInfo(pose)

cart = create_score_function("ref2015_cart")

#############################
## Setup Residue Selectors ##
#############################
#1 <CDR name="light_cdrs" cdrs="L1,L2,L3" />
#2  <CDR name="L1" cdrs="L1"/>
#3  <CDR name="L2" cdrs="L2"/>
#4  <CDR name="L3" cdrs="L3"/>
#5  <AntibodyRegion name="antigen" region="antigen_region" />
#6  <AntibodyRegion name="framework" region="framework_region" />
#7  <AntibodyRegion name="cdrs" region="cdr_region" />

#1 <CDR name="light_cdrs"
light_sele = CDRResidueSelector()
light_list = vector1_protocols_antibody_CDRNameEnum()
[light_list.append(x) for x in [l1, l2, l3]]

light_sele.set_cdrs(light_list)

#2 <CDR name="L1"
l1_sele = CDRResidueSelector()
l1_sele.set_cdr(l1)

#3 <CDR name="L2"
l2_sele = CDRResidueSelector()
l2_sele.set_cdr(l2)

#4 <CDR name="L3"
l3_sele = CDRResidueSelector()
l3_sele.set_cdr(l3)

#5 <AntibodyRegion name="antigen"
antigen_sele = AntibodyRegionSelector()
antigen_sele.set_region(antigen_region)

#6 <AntibodyRegion name="framework"
framework_sele = AntibodyRegionSelector()
framework_sele.set_region(framework_region)

#7 <AntibodyRegion name="cdrs"
cdrs_sele = AntibodyRegionSelector()
cdrs_sele.set_region(cdr_region)


##########################
## Setup Simple Metrics ##
##########################
#1 <SasaMetric name="sasa" residue_selector="light_cdrs" />
#2 <SelectedResiduesPyMOLMetric name="cdr_selection" residue_selector="light_cdrs" />
#3 <SequenceMetric name="L1_seq" residue_selector="L1" custom_type="L1"/>
#4 <SequenceMetric name="L2_seq" residue_selector="L2" custom_type="L2"/>
#4 <SequenceMetric name="L3_seq" residue_selector="L3" custom_type="L3"/>
#5 <InteractionEnergyMetric name="cdr-int" residue_selector="light_cdrs" residue_selector2="antigen" />

#1 <SasaMetric name="sasa"
sasa = SasaMetric(light_sele)

#2 <SelectedResiduesPyMOLMetric name="cdr_selection"
cdr_selection = SelectedResiduesPyMOLMetric(light_sele)
L1_seq = SequenceMetric(l1_sele)
L1_seq.set_custom_type("L1")

#3 <SequenceMetric name="L1_seq"
L2_seq = SequenceMetric(l2_sele)
L2_seq.set_custom_type("L2")

#4 <SequenceMetric name="L2_seq"
L3_seq = SequenceMetric(l3_sele)
L3_seq.set_custom_type("L3")

#5 <InteractionEnergyMetric name="cdr-int"
cdr_int = InteractionEnergyMetric(light_sele, antigen_sele)
cdr_int.set_scorefunction(cart)


########################## 
## Setup TaskOperations ##
##########################
#1 <RestrictToCDRsAndNeighbors name="restrict_to_cdrs" cdrs="L1,L2,L3" design_framework="1" design_cdrs="1"/>
#2 <AddCDRProfilesOperation name="profiles" cdrs="L1,L2,L3" fallback_strategy="CONSERVATIVE" /> 

light_vec_bool = vector1_bool(6)
light_vec_bool[l1] = True
light_vec_bool[l2] = True
light_vec_bool[l3] = True

#1 <RestrictToCDRsAndNeighbors name="restrict_to_cdrs"
restrict_to_cdrs = RestrictToCDRsAndNeighbors()
restrict_to_cdrs.set_cdrs(light_vec_bool)
restrict_to_cdrs.set_allow_design_neighbor_framework(True)
restrict_to_cdrs.set_allow_design_cdr(True)

#2 <AddCDRProfilesOperation name="profiles"
profiles = AddCDRProfilesOperation()
profiles.set_cdrs(light_vec_bool)
profiles.set_fallback_strategy(design.seq_design_conservative)


##########################
## Setup MovemapFactory ##
##########################
#1 <MoveMapFactory name="movemap_cdrs" bb="0" chi="0">
#1  <Backbone residue_selector="light_cdrs" />
#1  <Chi residue_selector="light_cdrs" />
#1 </MoveMapFactory>

#1
mmf = MoveMapFactory()
mmf.all_bb(False) #First, disable everything, then enable from actions.
mmf.all_chi(False)
mmf.add_bb_action(mm_enable, light_sele)
mmf.add_chi_action(mm_enable, light_sele)
mmf.set_cartesian(True)


##################
## Setup Movers ##
##################
#1 <CDRDihedralConstraintMover name="dih_mover_L1" cdr="L1" use_cluster_csts="1" />
#2 <CDRDihedralConstraintMover name="dih_mover_L2" cdr="L2" use_cluster_csts="1" />
#3 <CDRDihedralConstraintMover name="dih_mover_L3" cdr="L3" use_cluster_csts="1" />
#4 <AntibodyCDRGrafter name="grafter" cdrs="L1,L2,L3" input_ab_scheme="AHO_Scheme" cdr_definition="North" donor_structure_from_pdb="malaria_cdrs.pdb" use_secondary_graft_mover="1" optimize_cdrs="0" optimize_cdr4_if_neighbor="1" />
#5 <PackRotamersMover name="packrot" task_operations="restrict_to_cdrs,profiles" />

# Minimize Light chain CDRs using the MoveMapFactory.  We do it in cartesian to prevent lever-arm effects
#6 <MinMover name="minmover" cartesian="1" movemap_factory="movemap_cdrs" tolerance=".1"/>

# Run our metrics
#7 <RunSimpleMetrics name="design_metrics" metrics="sasa,cdr_selection,L1_seq,L2_seq,L3_seq,cdr-int" prefix="design_" />
#8 <RunSimpleMetrics name="native_metrics" metrics="sasa,cdr_selection,L1_seq,L2_seq,L3_seq,cdr-int"  prefix="native_" />

#1 <CDRDihedralConstraintMover name="dih_mover_L1"
dih_csts_l1 = CDRDihedralConstraintMover()
dih_csts_l1.set_cdr(l1)

#2 <CDRDihedralConstraintMover name="dih_mover_L2"
dih_csts_l2 = CDRDihedralConstraintMover()
dih_csts_l2.set_cdr(l2)

#3 <CDRDihedralConstraintMover name="dih_mover_L3"
dih_csts_l3 = CDRDihedralConstraintMover()
dih_csts_l3.set_cdr(l3)

#4 <AntibodyCDRGrafter name="grafter"
new_pose = pose_from_pdb("malaria_cdrs.pdb")

grafter = AntibodyCDRGrafter(ab_info)
grafter.set_donor_structure(new_pose)
grafter.set_cdrs(light_vec_bool)
grafter.set_use_secondary_graft_mover_if_needed(True)
grafter.set_optimize_cdrs(False)
grafter.set_include_cdr4_in_optimization(True)

#5 <PackRotamersMover name="packrot"
tf = TaskFactory()
tf.push_back(restrict_to_cdrs)
tf.push_back(profiles)

packer = PackRotamersMover()
packer.task_factory(tf)
packer.score_function(cart)

#6 <MinMover name="minmover"
minmover = MinMover()
minmover.cartesian(True)
minmover.movemap_factory(mmf)
minmover.tolerance(.1)
minmover.score_function(cart)

#7 <RunSimpleMetrics name="design_metrics" metrics="sasa,cdr_selection,L1_seq,L2_seq,L3_seq,cdr-int"
metrics = [sasa, cdr_selection, L1_seq, L2_seq, L3_seq, cdr_int]

design_metrics = RunSimpleMetricsMover()
design_metrics.set_prefix("design_")
[design_metrics.add_simple_metric(metric) for metric in metrics] #Save a few lines as the vector1 for metrics is not yet bound to PyRosetta.

#8 <RunSimpleMetrics name="native_metrics" metrics="sasa,cdr_selection,L1_seq,L2_seq,L3_seq,cdr-int"
native_metrics = RunSimpleMetricsMover()
native_metrics.set_prefix("native_")
[native_metrics.add_simple_metric(metric) for metric in metrics]


######################
## Run The Protocol ##
######################
#You also use a MoverContainer from previous tutorials here if you wanted to.

#1 Run Metrics on the native antibody
# <Add mover_name="native_metrics" />

#2 Graft Light chain CDRs from our donor into our current antibody.
# <Add mover_name="grafter" />

#3 Add Constraints to the CDRs so that min doesn't change them too much
# <Add mover_name="dih_mover_L1" />
# <Add mover_name="dih_mover_L2" />
# <Add mover_name="dih_mover_L3" />

#4 Design framework around light chain to accomodate the grafted CDRs.  Pack/Min/Pack
# <Add mover_name="packrot" />
# <Add mover_name="minmover" />
# <Add mover_name="packrot" /> 
# <Add mover_name="minmover" />

#5 Run Simple Metrics on the result
# <Add mover_name="design_metrics" />

if not os.getenv("DEBUG"):
    #1 Run Metrics on the native antibody
    native_metrics.apply(pose)

    #2 Graft Light chain CDRs from our donor into our current antibody.
    grafter.apply(pose)

    #3 Add Constraints to the CDRs so that min doesn't change them too much
    dih_csts_l1.apply(pose)
    dih_csts_l2.apply(pose)
    dih_csts_l3.apply(pose)

    #4 Design framework around light chain to accomodate the grafted CDRs.  Pack/Min/Pack
    packer.apply(pose)
    minmover.apply(pose)
    packer.apply(pose)
    minmover.apply(pose)

    #5 Run Simple Metrics on the result
    design_metrics.apply(pose)


# Optimizing dG
if not os.getenv("DEBUG"):
    pose = original_pose.clone()
    rabd = XmlObjects.static_get_mover('<AntibodyDesignMover name="RAbD" graft_design_cdrs="L1,L3" seq_design_cdrs="L1,L3" mc_optimize_dG="1" light_chain="kappa"/>')
    rabd.apply(pose)


# Optimizing Interface Energy and Total Score (opt-dG and opt-E)
if not os.getenv("DEBUG"):
    pose = original_pose.clone()
    rabd = XmlObjects.static_get_mover('<AntibodyDesignMover name="RAbD" seq_design_cdrs="L1,L3" graft_design_cdrs="L1,L3" mc_optimize_dG="1" mc_total_weight=".01" mc_interface_weight=".99" light_chain="kappa"/>')
    rabd.apply(pose)


# CDR Instruction File
if not os.getenv("DEBUG"):
    pose = original_pose.clone()
    rabd = XmlObjects.static_get_mover('<AntibodyDesignMover name="RAbD" instruction_file="my_instruction_file.txt" seq_design_cdrs="L1,L3,H1,H2,H3" graft_design_cdrs="L1,H2,H1" random_start="1" light_chain="kappa"/>')
    rabd.apply(pose)


# Mintype
if not os.getenv("DEBUG"):
    pose = original_pose.clone()
    rabd = XmlObjects.static_get_mover('<AntibodyDesignMover name="RAbD" seq_design_cdrs="L1,L3,H1,H3" mintype="relax" light_chain="kappa"/>')
    rabd.apply(pose)

