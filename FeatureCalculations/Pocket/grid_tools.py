#    MAHOMES II app
#    Copyright (C) 2021 University of Kansas
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published
#    by the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

# @file  
# @brief 
# @author Ryan Feehan <RFeehan93@gmail.com>

import os
import sys
import subprocess
import numpy as np
import scipy.spatial 

#added 11/16 MWF from CustomModules files
def rotation_matrix(vector1, origin):
    #vector1 should be a vector between two points; origin should be the vector to which vector1 is being moved
    vector1 = vector1/np.linalg.norm(vector1)
    origin = origin/np.linalg.norm(origin)
    v = np.cross(vector1, origin)
    c = np.dot(vector1, origin)
    s = np.linalg.norm(v)
    I = np.identity(3)
    k = np.matrix([(0, -v[2], v[1]), (v[2], 0, -v[0]), (-v[1], v[0], 0)])#.reshape(3,3)
    #print(k)
    r = I + k + np.matmul(k,k) * ((1 -c)/(s**2))
    return(r)

#added 11/16 MWF from CustomModules files
def translation(new_point, origin):
    return (new_point - origin)

#added 11/16 MWF from CustomModules files
def rotate_trans_coords(a, b, a_orig, b_orig, res_coord_array):
    #a, b are current points; a_orig, b_orig are where they should move to; res_coord_array contains coordinates to be moved
    #print(a_orig - b_orig)
    rot_mat = rotation_matrix(a-b, a_orig-b_orig)
    #print(rot_mat)
    trans_vect = translation(np.asarray(np.dot(rot_mat, a))[0], a_orig )
    #print(trans_vect)
    matrix_coords = np.transpose(res_coord_array)
    matrix_rot = np.matmul(rot_mat, matrix_coords)
    #print(matrix_rot)
    rotated_coords = np.transpose(matrix_rot)
    #print(rotated_coords)
    rotated_coords = rotated_coords - trans_vect
    return(rotated_coords) #returns a matrix of n by 3 of the translated/rotated vectors

#following functions blatantly copied from http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html
def fitEllipse(x,y):
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    #print("D:", D)
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2
    C[1,1] = -1
    #print("S:", S)
    #print("C:", C)
    #print( np.linalg.pinv(S) )
    E, V =  np.linalg.eig(np.dot(np.linalg.pinv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]
    return a

def ellipse_center(a):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])

def ellipse_angle_of_rotation( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    if b == 0:
        if a > c:
            return 0
        else:
            return np.pi/2
    else: 
        if a > c:
            return np.arctan(2*b/(a-c))/2
        else:
            return np.pi/2 + np.arctan(2*b/(a-c))/2

def ellipse_axis_length( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    #print(up, down1, down2)
    res1=np.sqrt(np.absolute(up/down1)) #added absolute value 1-10-18 MWF
    res2=np.sqrt(np.absolute(up/down2)) #added absolute value 1-10-18 MWF
    #print(res1, res2)
    return np.array([res1, res2])

#added 11/16 MWF

#given an array of points that have been pre-rotated to x,y,z and a pair of z-boundaries in that array
#return features describing that slice of points
#   farthest_point is the max of the distance between all the points in the slice
#   pocket_area is the area of the convex hull; in otherwords, the approx cross-sectional area of the slice
#   pocket_offset is how far the center of the ellipse is from the origin where the points have been rotated to
#   long/short_axis is the ellipse dimensions
def calc_slice_features(pocket_array, z_depth, interval = 1.4):
    upper_slice_boundary = z_depth + interval
    lower_slice_boundary = z_depth - interval
    #the z-coordinate is not considered for any of these calculations; the slice of points is flattened into just x,y space
    pocket_slice = pocket_array[ (pocket_array[:,2] <= upper_slice_boundary) & (pocket_array[:,2] >= lower_slice_boundary)]
    #print(pocket_slice)
    if len(pocket_slice) > 0:
        all_dist = scipy.spatial.distance.cdist(pocket_slice[:,0:2], pocket_slice[:,0:2])
        farthest_point = np.max(all_dist) 
        #print("farthest:", farthest_point)
        if len(pocket_slice) > 2:
            #calculate the area of the convex hull that encloses the points
            hull = scipy.spatial.ConvexHull(pocket_slice[:,0:2], qhull_options="QJ")
            #print(hull)
            #fit an ellipse to the hull points
            pocket_area = hull.area
            #print("pocket area", pocket_area)
            a = fitEllipse(hull.points[hull.vertices][:,0], hull.points[hull.vertices][:,1])
            center = ellipse_center(a)
            pocket_offset = scipy.spatial.distance.cdist([[0,0]], [center])
            phi = ellipse_angle_of_rotation(a)
            axes = ellipse_axis_length(a)
            long_axis = max(axes)
            short_axis = min(axes)
            #print(long_axis, short_axis)
            return [float(farthest_point), float(pocket_area), float(pocket_offset[0][0]), float(long_axis), float(short_axis)]
        else:
            return[farthest_point, 0, 0, 0, 0]
    else:
        return [0,0,0,0,0]

##Takes in a pocket pdb file and parses down to the relevant grid points    
def process_pocket(pocket_file):
    tp = []
    tpb = []
    tpe = []
    tps = []
    stp = []
    sts = []
    stb = []
    ste = []
    ts = []
    whole_pocket = []

    ##Read in the pocket PDB file
    ##Place the points into the respective lists. 
    ##TP should be key target. but, you can often have issues in pits so stp is useful too. 
    with open(pocket_file, "r") as pocketPDB:
        for entry in pocketPDB:
            if (entry[0:6] == "ATOM  "):
                m_id = int(entry[8:11])
                x = float(entry[31:38])
                y = float(entry[39:46]) 
                z = float(entry[47:54])
                resn = entry[13:16].strip() ##???
                t = np.array([x, y, z])

                if (resn == "PR"): #protein
                    pass
                elif (resn == "TP"): #target pocket
                    tp.append(t)
                    whole_pocket.append(t)
                elif (resn == "TPB"): #target pocket buried
                    tpb.append(t)
                    whole_pocket.append(t)
                elif (resn == "TPE"): #target pocket edge
                    tpe.append(t)
                    #whole_pocket.append(t)
                elif (resn == "TPS"): #target pocket surface
                    tps.append(t)
                    whole_pocket.append(t)
                elif (resn == "STP"): #small pocket
                    stp.append(t)
                    whole_pocket.append(t)
                elif (resn == "STS"): #small pocket surface
                    sts.append(t)
                    whole_pocket.append(t)
                elif (resn == "STB"): #small pocket buried
                    stb.append(t)
                    whole_pocket.append(t)
                elif (resn == "TS"): #target surface? 
                    ts.append(t)
                    whole_pocket.append(t)
                elif (resn == "STE"): #small pocket edge
                    ste.append(t)
                    #whole_pocket.append(t)
    return [tp, tpb, tpe, tps, stp, sts, stb, ts, whole_pocket] #each of these is a list of np arrays; these can be converted to an n by 3 array with np.asarray(X)

def get_vectors(pocket_defs, origin, closest_res):
    #pocket_defs should be a n by 3 numpy array
    #origin should be a [x, y, z] array
    #closest res is also in the form of an [x, y, z] array

    #Find average of subset of points selected to be the pocket opening, i.e. TP/STP residues in grid file
    ##Average of all points in pocket_defs 
    z_target = np.average(pocket_defs, axis = 0)
    #print(z_target)
    #Get Z
    z = z_target - origin
    z /= np.linalg.norm(z)
    #Get X
    x = closest_res - origin
    #print(x)
    #get the orthogonal vector from x
    x -= x.dot(z) * z

    ##Now the cross product to get the Y vector
    y = np.cross(x, z)

    #returns unit vectors of x, y, z
    return [x, y, z, z_target - origin]

def vectors_from_pdb(pocket_grid_name, metal_ion_center, m_closest_residue):
    pocket_defs = process_pocket(pocket_grid_name)
    endpoints = pocket_defs[0] + pocket_defs[4]
    vecs = get_vectors(np.asarray(endpoints), metal_ion_center, m_closest_residue)
 
    return vecs

#added 11-20-18 MWF
def calc_pocket_features(pdb_id, site_id, pocket, origin, nearest_res_coords, center_of_mass, directory):
    
    feature_labels = ["farPt", "PocketArea", "Offset", "LongAx", "ShortAx"]
    labels1 = [x+"Low" for x in feature_labels]
    labels2 = [x+"Mid" for x in feature_labels]
    labels3 = [x+"High" for x in feature_labels]
    all_labels = ["SEPocket", "Depth", "Vol", "LongPath"] + labels1 + labels2 + labels3

    target_pocket = pocket[0] + pocket[4] #TP and STP residue codes in pocket pdb - target pocket
    #print(target_pocket)
    if len(target_pocket) == 0:
        if len(pocket[3] + pocket[5]) > 0: #use TPS/STS - pocket surface
            target_pocket = pocket[3] + pocket[5]
            pocket_type = 2
            #print("USING TPS residues")
        #print(len(pocket[1]))
        elif len(pocket[1]) > 0: #if we have no TP or STP residues, use the TPB ones instead
            pocket_type = 1
            target_pocket = farthest_residue(center_of_mass, pocket[1]) #rather than the average of all of the TPB, 
            #print("USING TPB residues from C.O.M.", target_pocket)
        else:
            print("NO SURFACE POCKET", pdb_id, site_id)
            if len(pocket[7]) > 0 :
                return(all_labels, [0] * 19)
            else:
                return(all_labels, [4] + [0]*18 )
    else:
        pocket_type = 3
        
    vectors = get_vectors(np.asarray(target_pocket), origin, nearest_res_coords)
    #print(vectors)
    depth_in_pocket = np.linalg.norm(vectors[-1])
    volume = subprocess.check_output(["sed", "1q;d", "%s/%s_Vol.txt"%(directory, site_id)] )
    #print(volume)
    volume = volume.decode("utf-8").strip().split("\n")[0].split()[-1]
    volume = float(volume) #volume of pocket from pocket_grid
    longest_path_out = scipy.spatial.distance.cdist([vectors[-1]], [origin], "cityblock" )[0][0] #if you traverse a grid from origin to the target pocket, how far is it?

    #vectors[2] is the unit z-axis as returned from get_vectors, [0,0,0] is the origin for this z-axis; the pocket array must be pre-translated
    rotated_pocket = np.asarray(rotate_trans_coords(np.array([0,0,0]), vectors[2], np.array([0, 0, 0]), np.array([0,0,1]), np.asarray(pocket[-1])-origin))

    #the pocket grid is currently done on 0.5A intervals; 
    #an interval of 1.4 is approximately the distance between 3 stacked grid points on the diagonal; 
    #a slice containing +- an interval should capture a resonable collection of points
    middle_features = calc_slice_features(rotated_pocket, depth_in_pocket/2)
    upper_features = calc_slice_features(rotated_pocket, depth_in_pocket)
    lower_features = calc_slice_features(rotated_pocket, 0)
    
    this_SITE_features = [pocket_type, depth_in_pocket, volume, longest_path_out]
    this_SITE_features.extend(lower_features)
    this_SITE_features.extend(middle_features)
    this_SITE_features.extend(upper_features)
    
    #solvent-exposed pocket, depth, vol, long_path, [farthest_point, pocket_area, pocket_offset[0][0], long_axis, short_axis] for low/mid/upper
    return(all_labels, this_SITE_features)
    
def nearest_neigh(site_center, metal_resIDs, myProtein):
    distances, near_res = myProtein.KDTree.query(site_center, k = 40, distance_upper_bound = 10)
    distances = [x for x in distances if np.isinf(x) == False]
    if len(distances) == 0:
        return("NA")
        
    near_res = near_res[0:len(distances)]
    new_res_ids = [myProtein.Coords_index[x] for x in near_res]
    new_res_ids = sorted(set(new_res_ids), key = new_res_ids.index)
    new_res_ids = [x for x in new_res_ids if x not in metal_resIDs]
    #print(new_res_ids)
    if len(new_res_ids) == 0:
        return("NotBound")
    else:
        for x in range(0, len(new_res_ids)):
            this_res = myProtein.residues[ myProtein.res_nums.index(new_res_ids[x]) ]
            #print(new_res_ids[x], this_res.resnum, this_res.type)
            if this_res.type == "protein":
                return(new_res_ids[x])
        return("NotBound")

def closest_residue(origin, residue_list):
    distances = residue_list - origin
    nearest_neighbor = residue_list[np.argmin( np.sum(distances**2, axis=1))] #return the row that corresponds to the minimum distance
    return nearest_neighbor
    
def farthest_residue(origin, residue_list):
    distances = residue_list - origin
    nearest_neighbor = residue_list[np.argmax( np.sum(distances**2, axis=1))] #return the row that corresponds to the minimum distance
    return nearest_neighbor



#adjusted version for using ghecom pocket
def calc_pocket_features_ghecom(pocket, origin, metal_res, this_protein, res_nums):
    site_neigh = nearest_neigh(origin, metal_res, this_protein)
    nearest_neigh_coords = closest_residue(origin, this_protein.residues[ res_nums.index(int(site_neigh[:-1])) ].Coords)

    feature_labels = ["farPt", "PocketArea", "Offset", "LongAx", "ShortAx"]
    labels1 = [x+"Low" for x in feature_labels]
    labels2 = [x+"Mid" for x in feature_labels]
    labels3 = [x+"High" for x in feature_labels]
    all_labels = ["Depth", "Vol", "LongPath"] + labels1 + labels2 + labels3

    target_pocket = pocket
        
    vectors = get_vectors(np.asarray(target_pocket), origin, nearest_neigh_coords)
    #print(vectors)
    depth_in_pocket = np.linalg.norm(vectors[-1])
    
    volume = len(pocket)*0.512
    longest_path_out = scipy.spatial.distance.cdist([vectors[-1]], [origin], "cityblock" )[0][0] #if you traverse a grid from origin to the target pocket, how far is it?

    #vectors[2] is the unit z-axis as returned from get_vectors, [0,0,0] is the origin for this z-axis; the pocket array must be pre-translated
    rotated_pocket = np.asarray(rotate_trans_coords(np.array([0,0,0]), vectors[2], np.array([0, 0, 0]), np.array([0,0,1]), np.asarray(pocket)-origin))

    #the pocket grid is currently done on 0.5A intervals; 
    #an interval of 1.4 is approximately the distance between 3 stacked grid points on the diagonal; 
    #a slice containing +- an interval should capture a resonable collection of points
    middle_features = calc_slice_features(rotated_pocket, depth_in_pocket/2)
    upper_features = calc_slice_features(rotated_pocket, depth_in_pocket)
    lower_features = calc_slice_features(rotated_pocket, 0)
    
    this_SITE_features = [depth_in_pocket, volume, longest_path_out]
    this_SITE_features.extend(lower_features)
    this_SITE_features.extend(middle_features)
    this_SITE_features.extend(upper_features)
    
    #solvent-exposed pocket, depth, vol, long_path, [farthest_point, pocket_area, pocket_offset[0][0], long_axis, short_axis] for low/mid/upper
    return(all_labels, this_SITE_features)
