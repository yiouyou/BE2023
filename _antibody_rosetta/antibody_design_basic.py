import pyrosettacolabsetup
import pyrosetta

from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.teaching import *

from pyrosetta.rosetta.protocols.rosetta_scripts import *
from pyrosetta.rosetta.protocols.antibody import *
from pyrosetta.rosetta.protocols.antibody.design import *
from pyrosetta.rosetta.utility import *

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
init('-no_fconfig @/home/songz/PyRosetta.notebooks/notebooks/inputs/rabd/common')

pose = pose_from_pdb("/home/songz/PyRosetta.notebooks/notebooks/inputs/rabd/my_ab.pdb")
original_pose = pose.clone()


# Sequence Design
rabd = XmlObjects.static_get_mover('<AntibodyDesignMover name="RAbD" seq_design_cdrs="L1,L3" light_chain="kappa"/>')
if not os.getenv("DEBUG"):
    rabd.apply(pose)

pose = original_pose.clone()
rabd2 = AntibodyDesignMover()
cdrs = vector1_protocols_antibody_CDRNameEnum()
cdrs.append(l1)
cdrs.append(l3)
rabd2.set_seq_design_cdrs(cdrs)
rabd2.set_light_chain("kappa")
if not os.getenv("DEBUG"):
    rabd2.apply(pose)

from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

if not os.getenv("DEBUG"):
    iam = InterfaceAnalyzerMover("LH_ABCDEFGIJKZ")
    iam.set_pack_separated(True)
    iam.apply(pose)
    iam.apply(original_pose)
    dg_term = "dG_separated"
    print("dG Diff:", pose.scores[dg_term] - original_pose[dg_term])

df = load_json_scorefile("/home/songz/PyRosetta.notebooks/notebooks/expected_outputs/rabd/tutA1_score.sc")
df = drop_cluster_columns(df, keep_cdrs=["L1", "L3"])
print(df)


# Graft Design
df_a2 = load_json_scorefile("/home/songz/PyRosetta.notebooks/notebooks/expected_outputs/rabd/tutA2_score.sc")
df_a2 = drop_cluster_columns(df_a2, keep_cdrs=["L1", "L3"])
print(df_a2)

df_tut_a12 = pandas.concat([df, df_a2], ignore_index=True).sort_values("dG_separated", ascending=True)
print(df_tut_a12)

if not os.getenv("DEBUG"):
    print("L1", original_pose.scores["cdr_cluster_ID_L1"])
    print("L3", original_pose.scores["cdr_cluster_ID_L3"])


# Basic De-novo run
if not os.getenv("DEBUG"):
    pose = original_pose.clone()
    rabd = XmlObjects.static_get_mover('<AntibodyDesignMover name="RAbD" graft_design_cdrs="L1,L3" seq_design_cdrs="L1,L3" random_start="1"/> light_chain="kappa"')
    rabd.apply(pose)
    # OR (REUSE code from above)
    pose = original_pose.clone()
    rabd2.set_seq_design_cdrs(cdrs)
    rabd2.set_graft_design_cdrs(cdrs)
    rabd2.set_random_start(True)
    rabd2.set_light_chain("kappa")
    rabd2.apply(pose)


# RAbD Framework Components
if not os.getenv("DEBUG"):
    pose = original_pose.clone()
    parser = RosettaScriptsParser()
    protocol = parser.generate_mover_and_apply_to_pose(pose, "/home/songz/PyRosetta.notebooks/notebooks/inputs/rabd/ab_design_components.xml")
    protocol.apply(pose)

