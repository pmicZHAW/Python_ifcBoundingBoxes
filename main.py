# 21.01.2022 from Adne:
# - To install trimesh in python, run "python -m pip install trimesh[easy]"
# - Attached is the json file we discussed at the meeting with Dag. The structure is a python dict where the keys is the
#   guids in the IFC-file and each value is a dict with a transformation matrix ("T") and size of the box ("S").
#   I also attached the ifc-file used to generate the oriented bounding boxes.

# --- imports ----------------------------------------------------------------------------------------------------------

import trimesh
import json
import trimesh.creation
import numpy as np
import math
from scipy.spatial.transform import Rotation
import pandas as pd

# --- define some parameters -------------------------------------------------------------------------------------------

# file_name = "20220120_zhaw_lab/objects.json"
# file_name = "20220207_zhaw_lab/objects.json"
file_name = "20220215_robotic_lab/objects.json"

rpy0 = np.array([0, 0, 90]) * math.pi / 180 * 0
p0 = np.array([[1], [2], [3]]) * 0


# --- custom functionality (Michael is lazy) ---------------------------------------------------------------------------

def getTransform(R, p):
    zero = np.zeros((1, 3))
    one = np.eye(1)
    T = np.concatenate((np.concatenate((R, p), axis=1), np.concatenate((zero, one), axis=1)), axis=0)
    return np.squeeze(np.asarray(T))


def getTransformFromEuler(rpy, p):
    R = Rotation.from_euler('xyz', rpy, degrees=False).as_matrix()
    return getTransform(R, p)


def getRotationMatrixAndTranslationVector(T):
    R = T[0:3, 0:3]
    p = np.matrix(T[0:3, 3]).transpose()
    return R, p


def getInverseTransform(T):
    R, p = getRotationMatrixAndTranslationVector(T)
    zero = np.zeros((1, 3))
    one = np.eye(1)
    T = np.concatenate(
        (np.concatenate((R.transpose(), -R.transpose() @ p), axis=1), np.concatenate((zero, one), axis=1)), axis=0)
    return np.squeeze(np.asarray(T))


# --- doing stuff ------------------------------------------------------------------------------------------------------

T0 = getTransformFromEuler(rpy0, p0)

with open(file_name, "r") as file:
    bounds = json.load(file)

# solid boxes
solids = []
for k, v in bounds["solids"]:  # bounds.items():
    solids.append(trimesh.creation.box(v["S"], v["T"]))

solid_cntr = 0
solid_attributes = pd.DataFrame(columns=("xMin", "xMax", "yMin", "yMax", "zMin", "zMax", 'volume'))
for box in solids:
    box.apply_transform(T0)
    solid_attributes.loc[solid_cntr] = [np.amin(np.array(box.vertices[:, 0])), np.amax(np.array(box.vertices[:, 0])),
                                        np.amin(np.array(box.vertices[:, 1])), np.amax(np.array(box.vertices[:, 1])),
                                        np.amin(np.array(box.vertices[:, 2])), np.amax(np.array(box.vertices[:, 2])),
                                        box.volume]
    solid_cntr += 1

# boxes to open the solid boxes
openings = []
for k, v in bounds["openings"]:
    openings.append(trimesh.creation.box(v["S"], v["T"]))

opening_cntr = 0
for box in openings:
    box.apply_transform(T0)
    opening_cntr += 1

# remove ceiling, we assume it is the box with the largest volume out of the N_ceiling_candidate solids with the highest vertices
N_ceiling_candidates = 10
zMax_list = np.flipud(np.argsort(solid_attributes["zMax"].to_numpy()))
ceiling_idx = zMax_list[np.argmax(solid_attributes["volume"].to_numpy()[zMax_list[0:N_ceiling_candidates]])]
solids.pop(ceiling_idx)

# remove floor, we assume it is the box with the largest volume out of the N_floor_candidate solids with the lowest vertices
N_floor_candidates = 10
zMin_list = np.argsort(solid_attributes["zMax"].to_numpy())
floor_idx = zMin_list[np.argmax(solid_attributes["volume"].to_numpy()[zMin_list[0:N_floor_candidates]])]
floor_zVals = np.array(solids[floor_idx].vertices[:, 2])
zMax_list = np.flipud(np.argsort(floor_zVals))
floor_zVal = floor_zVals[zMax_list[0]]
solids.pop(floor_idx)

# # calculate the overall bounding box
# all_vertices = np.empty((0, 3))  # all_vertices are needed to calculate the overall bounding box
# for box in solids:
#     all_vertices = np.append(all_vertices, np.array(box.vertices), axis=0)
# T, S = trimesh.bounds.oriented_bounds(all_vertices, angle_digits=1, ordered=True, normal=None)
# T = getInverseTransform(T)
# bounding_box = trimesh.creation.box(S, T)

# remove anything that is below slice_zMin and above slice_zMax relative to floor_zVal
slice_zMin = 0.01
slice_zMax = 1.2
slice_points = [[0, 0, floor_zVal + slice_zMin], [0, 0, floor_zVal + slice_zMax]]
slice_normals = [[0, 0, 1], [0, 0, -1]]
solids_sliced_z = []
for box in solids:
    solids_sliced_z.append(box.slice_plane(slice_points, slice_normals))

# # remove the openings
# solids_sliced = []
# for solid in solids_sliced_z:
#     for opening in openings:
#         # solids_sliced.append(solid.slice_plane(opening.facets_origin, opening.facets_normal))
#         facets_origin = opening.facets_origin
#         facets_normal = opening.facets_normal
#         for i in np.arange(len(facets_origin)):
#             solids_sliced.append(solid.slice_plane(facets_origin[i], facets_normal[i]))
solids_sliced = solids_sliced_z

scene = trimesh.Scene()
# scene.add_geometry(solids_sliced)
scene.add_geometry(solids)
scene.add_geometry(trimesh.creation.axis())
scene.show()

# # --- example for slicing starts here ----------------------------------------------------------------------------------
# m = trimesh.creation.icosphere()
# # box = trimesh.creation.box(extents=[1.5, 1.5, 1.5])
# box = trimesh.creation.box(extents=[3, 3, 3])
# result = m.slice_plane(box.facets_origin, -box.facets_normal)
# scene = trimesh.Scene()
# scene.add_geometry(result)
# scene.add_geometry(trimesh.creation.axis())
# scene.show()

# # --- example for box creation and rotation starts here ----------------------------------------------------------------
#
# # extents: extents S = [sx, sy, sz]
# #   from the center of body frame (see T) the box is extended si/2 towards
# #   the positive and the negative axis i for i = x, y, z
# S = [2, 1, 0.5]
#
# # transform: transformation matrix T = [[R, p], [0, 1]]
# #   R descibes the rotation and p the translation of the body frame (B) w.r.t. the global frame (G)
# #   e.g. G_corner = C_GB * B_corner + G_p_GB where C_GB = R and G_p_GB points from (G) to (B) w.r.t. (G)
# ex = [1, 0, 0]
# ey = [0, 1, 0]
# ez = [0, 0, 1]
# I = np.eye(3)
#
# # R <- rotation matrix from euler angles rpy = [-30, 45, 60]*pi/180
# R = np.array([[0.35355339, -0.92677670, -0.12682648],
#               [0.61237244, 0.12682648, 0.78033009],
#               [-0.70710678, -0.35355339, 0.61237244]])
# # p <- an arbitary shift
# p = np.array([[1], [2], [3]])
#
# zero = np.zeros((1, 3))
# one = np.eye(1)
# T = np.concatenate((np.concatenate((R, p), axis=1), np.concatenate((zero, one), axis=1)), axis=0)
#
# # build the box and the transformed (translated and rotated) box
# box = trimesh.creation.box(S, np.eye(4))  # T = I
# t_box = trimesh.creation.box(S, T)
#
# print("box --------------------------------")
# print(f"box.extends {box.extents}")
# print(f"box.centroid {box.centroid}")  # center of the box
# print(f"box.bounds {box.bounds}")
# print(f"box.vertices {box.vertices}")  # vertices are the edges of the box
# print(f"box.area {box.area}")
# print(f"box.area_faces {box.area_faces}")
# print(f"box.euler_number {box.euler_number}")
# print(f"box.volume {box.volume}")  # volume
#
# print(f"box projected yz {box.projected(I[:, 0]).bounds}")  # ex
# print(f"box projected xz {box.projected(I[:, 1]).bounds}")  # ey
# print(f"box projected xy {box.projected(I[:, 2]).bounds}")  # ez
#
# sphere_ = box.bounding_sphere  # ?
#
# print("box transformed --------------------------------")
# print(f"t_box.extends {t_box.extents}")
# print(f"t_box.centroid {t_box.centroid}")  # center of the box
# print(f"t_box.bounds {t_box.bounds}")
# print(f"t_box.vertices {t_box.vertices}")  # vertices are the edges of the box
# print(f"t_box.area {t_box.area}")
# print(f"t_box.area_faces {t_box.area_faces}")
# print(f"t_box.euler_number {t_box.euler_number}")
# print(f"t_box.volume {t_box.volume}")  # volume
#
# print(f"t_box projected yz {t_box.projected(I[:, 0]).bounds}")  # ex
# print(f"t_box projected xz {t_box.projected(I[:, 1]).bounds}")  # ey
# print(f"t_box projected xy {t_box.projected(I[:, 2]).bounds}")  # ez
#
# boxes = []
# boxes.append(box)
# boxes.append(t_box)
#
# trimesh.Scene(boxes).show()
