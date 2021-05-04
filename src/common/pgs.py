'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: Pilgrim
Version     : 2021.3
License     : MIT/x11

Copyright (c) 2021, David Ferro Costas (david.ferro@usc.es) and
Antonio Fernandez Ramos (qf.ramos@usc.es)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software
is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
---------------------------

*----------------------------------*
| Module     :  common             |
| Sub-module :  pgs                |
| Last Update:  2020/02/03 (Y/M/D) |
| Main Author:  A. Fernandez-Ramos |
*----------------------------------*

This module contains some functions related
the calculation of the point group of symmetry
'''

#=============================================#
import numpy as np
#=============================================#

# /////  Matrix representation of symmetry operators \\\\\
def inversion():
    i_center = ([[-1.,0.,0.],
                 [0.,-1.,0.],
                 [0.,0.,-1.]])
    return i_center

def reflex(uvw,main_plane):
  # reflection about a plane perpendicular to uvw
    plane = np.zeros( (3,3) )
    if main_plane == 'xy': uvw = [0.,0.,1.]
    if main_plane == 'xz': uvw = [0.,1.,0.]
    if main_plane == 'yz': uvw = [1.,0.,0.]
    x = uvw[0]
    y = uvw[1]
    z = uvw[2]
    plane[0,0] = 1.-2.*x**2 
    plane[1,1] = 1.-2.*y**2 
    plane[2,2] = 1.-2.*z**2 
    plane[0,1] = -2.*x*y
    plane[0,2] = -2.*x*z
    plane[1,2] = -2.*y*z
    plane[1,0] = plane[0,1] 
    plane[2,0] = plane[0,2] 
    plane[2,1] = plane[1,2] 
    return plane
  
def Cngen(n,u):
  # rotacion about the u vector
    Cngen = np.zeros( (3,3) )
    t=2.*np.pi/float(n)
    ux = u[0]
    uy = u[1]
    uz = u[2]
    cosu = np.cos(t)
    sinu = np.sin(t)
    Cngen[0,0] = cosu + ux**2*(1.-cosu)
    Cngen[1,1] = cosu + uy**2*(1.-cosu)
    Cngen[2,2] = cosu + uz**2*(1.-cosu)
    Cngen[0,1] = ux*uy*(1.-cosu)-uz*sinu
    Cngen[0,2] = ux*uz*(1.-cosu)+uy*sinu
    Cngen[1,2] = uy*uz*(1.-cosu)-ux*sinu
    Cngen[1,0] = ux*uy*(1.-cosu)+uz*sinu
    Cngen[2,0] = ux*uz*(1.-cosu)-uy*sinu
    Cngen[2,1] = uy*uz*(1.-cosu)+ux*sinu
    return Cngen

def cdmass(mass,xyz,natom):
  # center of masses
    cdm = [0.,0.,0.]
    tot_mass = sum(mass)
    for i in range(0,natom):
      cdm[0] += mass[i]*xyz[i,0]/tot_mass
      cdm[1] += mass[i]*xyz[i,1]/tot_mass
      cdm[2] += mass[i]*xyz[i,2]/tot_mass
    for i in range(0,natom):
      xyz[i,0] += - cdm[0]
      xyz[i,1] += - cdm[1]
      xyz[i,2] += - cdm[2]
    return xyz
  
def I_diff(fx,fy,tol):
    abs_rel = np.fabs(fx-fy)/fy*100.
    if abs_rel <= tol: 
      return True
    else:
      return False
  
def calc_tensor(mass,xyz,natom):
  # tensor of inertia 
    tensor_iner = np.zeros( (3,3) )
    for i in range(0,natom):
      tensor_iner[0,0] += mass[i]*(xyz[i,1]**2+xyz[i,2]**2)
      tensor_iner[1,1] += mass[i]*(xyz[i,0]**2+xyz[i,2]**2)
      tensor_iner[2,2] += mass[i]*(xyz[i,0]**2+xyz[i,1]**2)
      tensor_iner[0,1] += -mass[i]*xyz[i,0]*xyz[i,1]
      tensor_iner[0,2] += -mass[i]*xyz[i,0]*xyz[i,2]
      tensor_iner[1,2] += -mass[i]*xyz[i,1]*xyz[i,2]
    tensor_iner[1,0] = tensor_iner[0,1]
    tensor_iner[2,0] = tensor_iner[0,2]
    tensor_iner[2,1] = tensor_iner[1,2]
    #eigvals, eigvecs = np.linalg.eigvalsh(tensor_iner)
    Im,Ivec =  np.linalg.eigh(tensor_iner)
    Ivec = Ivec.transpose()
    return Im,Ivec


def dist_mat(xyz,natom):
   # Matrix of distances
   dmx = []
   for i in range(0,natom):
     vtmp = xyz[i]
     for j in range (0,natom):
       vtmp2 = xyz[j]
       vx=vtmp-vtmp2
   # Norm + storage
       dmx.append(np.linalg.norm(vx))
   dmx=np.reshape(dmx,(natom,natom))
   return dmx
  
def dist_c(v1,v2):
    L_novec = False
    vxn = []
    vx = (v1+v2)/2
    vxn=np.linalg.norm(vx)
    if vxn == 0.: 
      L_novec = True
      return vx,L_novec
    vx /= vxn
    return vx,L_novec 
     
def dist_v(v1,v2):
    L_nosig = False
    vxn = []
    vx = (v1-v2)
    vxn=np.linalg.norm(vx)
    if vxn == 0.: 
      L_nosig = True
      return vx,L_nosig
    vx /= vxn
    return vx,L_nosig
    
  
def get_sea(dict_sea,mat_geo,idx,atom_mass):
    # Mass and geometry of each SEA
    coor_sea = []
    mass_sea = []
    idum=len(dict_sea[idx])
    for i in range(idum):
          k = dict_sea[idx][i]
          coor_sea.append(mat_geo[k])
          mass_sea.append(atom_mass[k])
    coor_sea=np.reshape(coor_sea,(idum,3))
    return idum,mass_sea,coor_sea
  
def compare_geom(A,B):
    nrows, ncols = A.shape
    tot_error = 0.0
    indices_comparison = []

    for row_i in range(nrows):
        Arow = A[row_i,:]
        diff    = float("inf")
        idx     = None
        for row_j in range(nrows):
            Brow = B[row_j,:]
            error = np.linalg.norm(Arow - Brow)
            if error < diff:
               diff = error
               idx  = row_j
        tot_error += diff
        indices_comparison.append( (row_i,idx)  )
    tot_error = 1.0 * tot_error / nrows
    return tot_error 

def get_sym(dict_sea,mat_geo,mat_sym,num_set,strsym,tolsym,atom_mass):
    # returns True or False, i.e. if the molecule is invariant after a given
    # symmetry operation (strsym) 
    total_error = 0.
    mat_tmp = mat_geo.T 
    mat_tmp = mat_sym*mat_tmp
    mat_tmp = mat_tmp.T
    L_sym = False
    for i in range(num_set):
      idx = i 
      isea,mass_sea,coor_sea =  get_sea(dict_sea,mat_geo,idx,atom_mass)
      isea,mass_tmp,coor_tmp =  get_sea(dict_sea,mat_tmp,idx,atom_mass)
      coor_sea = np.array(coor_sea)
      error_sym_oper = compare_geom(coor_sea,coor_tmp)
      total_error += error_sym_oper
    total_error /= num_set
    if total_error < tolsym: L_sym = True

    return L_sym  
  
def get_c2(dict_sea,mat_geo,num_set,uvw,L_round,epsilon,atom_mass,dxm):
    #
    # uvw is for all type of molecules except spheric
    # L_round = True for spherical molecules 
    # search for C2 axis: 1) passing through the middle of a bond
    #                     2) passing through atoms
    coor_sea = []
    mass_sea = []
    new_c2 = 0
    new_c2_red = 0  
    list_tot = []
    list_red = []
    dxc2 = []

    for l in range(num_set):
      idx = l
      isea,mass_sea,coor_sea = get_sea(dict_sea,mat_geo,idx,atom_mass)
      if isea > 1:
        for i in range(isea):
          for j in range(i+1,isea):
            dx,L_novec = dist_c(coor_sea[i],coor_sea[j])
            if L_novec:
              continue
            elif L_round:
              dxm.append(dx)
              new_c2_red += 1
            else:
              vcros = np.cross(dx,uvw) 
              vcrosn=np.linalg.norm(vcros)
              if vcrosn > 1.-epsilon: 
                new_c2_red += 1
                dxm.append(dx)
          vnorm=np.linalg.norm(coor_sea[i])
          if vnorm > epsilon:
            dxm.append(coor_sea[i]/vnorm)
            new_c2_red += 1
  #
  # ---remove the redundancies from the C2 axes
  #
    for i in range(new_c2_red): list_tot.append(i)  
    set_tot = set(list_tot)
    if new_c2_red > 1:
      for i in range(new_c2_red):
        vdxi = dxm[i]
        for j in range(i+1,new_c2_red): 
          vdxj = dxm[j]
          vcros = np.cross(vdxi,vdxj) 
          vcrosn=np.linalg.norm(vcros)
          if j in list_red:
            continue
          else:
            if vcrosn < epsilon: 
              list_red.append(j)
    set_red = set(list_red)
    set_not_red = set_tot - set_red
    dxc2.append([dxm[i] for i in set_not_red])
    new_c2 = len(set_not_red)
    dxc2=np.reshape(dxc2,(new_c2,3))
    return new_c2,dxc2

def get_sigma_d(dict_sea,mat_geo,num_set,mat_c2_axis,numc2,tolsym,atom_mass):
    '''
    Diagonal planes: considers that are perpendicular to the
      vectors defined by the C2 axes and the coordinate of
      an atom of the SEA 
      Important for Td molecules
    '''
    coor_sea = []
    mass_sea = []
    mat_sigma_d = [] 
    num_sigma_d = 0
    vdxi = []
    vdxj = []
    for l in range(num_set):
      idx = l
      isea,mass_sea,coor_sea = get_sea(dict_sea,mat_geo,idx,atom_mass)
      if isea > 1:
        for i in range(isea):
          vdxi = coor_sea[i]
          for j in range(numc2):
            vdxj = mat_c2_axis[j]
            vnormal = np.cross(vdxi,vdxj)
            vnorm=np.linalg.norm(vnormal)
            if vnorm > tolsym:
              vnormal /= vnorm 
              mat_sigma_d = reflex(vnormal,' ')
              L_sigma_d = get_sym(dict_sea,mat_geo,mat_sigma_d,num_set,'sigma_d',tolsym,atom_mass)
              if L_sigma_d:
                num_sigma_d += 1
            else:
              continue
    return num_sigma_d 
# ///// END \\\\

def get_pgs(atom_num,atom_mass,geom_xyz,toldist=0.05,tolsym=3e-2,epsilon=3e-2):
  '''
  This module finds the symmetry point
  group symmetry of a molecule 
  
  *---------------------------------------*
  | Main Author:  Antonio Fernandez-Ramos |
  | Last Update:  May 10th 2017 (by DFC)  |
  *---------------------------------------*
  '''
  
  #  -------Parameters----------
  # atomic numbers
  # atomic masses
  # Geometry as [[x1,y1,z1],[x2,y2,z2],[x3,y3,z3] ....
  # toldist : tolerance of difference of distances
  # tolsym : tolerance in the symmetry operation
  # epsilon : tolerance in the cross products
  #           usually used to check if two C2 axis are parallel

  
  # in case geom_xyz is a 3N list 
  if len(geom_xyz) == 3*len(atom_num):
     xyz = []
     for idx in range(0,len(geom_xyz),3):
         x,y,z = geom_xyz[idx+0:idx+3]
         xyz.append( [x,y,z] )
     geom_xyz = xyz

  natom = len(atom_num)

  if natom == 1: return 'K',1 # point group of a Sphere
  
    
  
  #####initial variables##################
  
  L_cn = False
  L_c2_p = False
  L_c2_sp = False
  L_sigma_v = False
  L_sigma_h = False
  L_sigma_x = False
  L_i = False
  L_s2n = False
  L_cs = False
  cnrot = []
  cn_from_I = [1]
  sigma_h = []
  sigma_v = []
  sigma_x = []
  sn_mat = []
  i_center = inversion()
  udum = []
  uvw = [0.,0.,0.]
  axis = ['x','y','z']
  plane = ['yz','xz','xy']
  dxm = [] 
  
  
  mat_geo = np.matrix(geom_xyz)
  # Puts the molecule in the center of masses
  mat_geo = cdmass(atom_mass,mat_geo,natom)
  # Evaluates the tensor of inertia
  Ipmol,Ipvec = calc_tensor(atom_mass,mat_geo,natom)
  # Rotates the molecule so Ia, Ib and Ic coincide with i,j,k
  mat_geo = Ipvec * mat_geo.T
  mat_geo = mat_geo.T
  
  
  dist = dist_mat(mat_geo,natom)
  distv = []
  # sum the distances of each row
  for i in range(natom):
     rkk = sum(dist[i,:])/float(natom)
     distv.append(rkk)
  
  # now it checks for possible SEA
  # index of SEAs stored in the dictionary
  dict_sea = { }
  idx = -1
  sea_at_idx = []
  for i in range(natom):
     next = False
     if i not in sea_at_idx: 
       idx += 1
     for j in range(i,natom):
       
       if atom_mass[i] == atom_mass[j] and atom_num[i] == atom_num[j] and np.fabs(distv[i] - distv[j])< toldist and j not in sea_at_idx:
         sea_at_idx.append(j)
         dict_sea[idx]  = dict_sea.get(idx,[]) + [j]
  
  # number of SEAs
  # maximum number of atoms in SEA
  num_set = idx+1
  atoms_sea = []
  for i in range(num_set):
    atoms_sea.append(len(dict_sea[i])) 
  atom_sea_max = max(atoms_sea)
  
  
  # Is it linear?
  if np.fabs(Ipmol[0]) <= epsilon:
  # Linear molecule
     linear = True
     L_i=get_sym(dict_sea,mat_geo,i_center,num_set,'inversion',tolsym,atom_mass)
     if L_i == True:
       return 'Dinfv',2
     else:
       return 'Cinfv',1
  
  # Non linear molecule
      # spheric symmetry (Ia = Ib = Ic)
      # achatado : oblate (Ia = Ib < Ic)
      # alargado : prolate (Ia < Ib = Ic)
      # asymmetric (Ia != Ib != Ic)
  L_sphera = False
  L_asym = False
  L_oblate = I_diff(Ipmol[0],Ipmol[1],epsilon*10)
  L_prolate = I_diff(Ipmol[1],Ipmol[2],epsilon*10)
  if L_oblate and L_prolate:
    L_sphera = True
    L_asym = False
  if L_sphera:
    L_prolate = False
    L_oblate = False
  if not L_oblate and not L_prolate and not L_sphera: 
    L_asym = True

  #
  # --- SPHERICAL MOLECULES
  #
  if L_sphera: 
    # ---C2 axes
    uvw_asym = [0.,0.,0]
    new_c2,dxc2 = get_c2(dict_sea,mat_geo,num_set,uvw_asym,L_sphera,epsilon,atom_mass,dxm)
    c2sp = []
    dxc2_save = []
    number_c2 =0
    for i in range(new_c2):
      c2sp = Cngen(2,dxc2[i])
      L_c2_sp=get_sym(dict_sea,mat_geo,c2sp,num_set,'C2',tolsym,atom_mass)
      if L_c2_sp:
        number_c2 += 1
        dxc2_save.append(dxc2[i]) 
    # ---diagonal planes
    if number_c2 <= 3:
      num_sigma_d = get_sigma_d(dict_sea,mat_geo,num_set,dxc2_save,number_c2,tolsym,atom_mass)
    # ---inversion
    L_i=get_sym(dict_sea,mat_geo,i_center,num_set,'inversion',tolsym,atom_mass)
    if number_c2 == 3:
      if L_i:
        return 'Th',12
      elif num_sigma_d >= 6: 
        return 'Td',12 
      else:
        return 'T',12 
    if number_c2 == 9:
      if L_i:
        return 'Oh',24
      else:
        return 'O',24 
    if number_c2 > 9:
      if L_i:
        return 'Ih',60
      else:
        return 'I',60 
  
  #
  # --- SYMMETRIC MOLECULES (PROLATE or OBLATE)
  #
  #      Cn axis search
  #
  if L_oblate or L_prolate:
      L_cn = False
      L_cn_tmp = False
      L_cn2_tmp = False
      if L_oblate: 
         main_axis = axis[2]
         main_plane = plane[2]
         uvw = [0.,0.,1.]
      if L_prolate: 
         main_axis = axis[0]
         main_plane = plane[0]
         uvw = [1.,0.,0.]
      for nrot in range(6,1,-1):
        cn_axis = nrot
        cad_axis = 'c'+str(cn_axis)+main_axis
        cnrot =  Cngen(cn_axis,uvw)
        L_cn=get_sym(dict_sea,mat_geo,cnrot,num_set,cad_axis,tolsym,atom_mass)
        if L_cn: 
          max_cn_axis = cn_axis
          break
      if not L_cn:
        L_asym = True
  
  #
  # --- ASYMMETRIC MOLECULES
  #
  #      C2 axis search
  #
  if L_asym:
    max_cn_axis = 1
    for iaxis in range(3):
      init_axis = axis[iaxis]
      init_plane = plane[iaxis]
      if iaxis == 0: uvw = [1.,0.,0.]
      if iaxis == 1: uvw = [0.,1.,0.]
      if iaxis == 2: uvw = [0.,0.,1.]
      cn_axis = 2
      cad_axis = 'c'+str(cn_axis)+init_axis
      cnrot = []
      cnrot =  Cngen(cn_axis,uvw)
      L_cn_tmp=get_sym(dict_sea,mat_geo,cnrot,num_set,cad_axis,tolsym,atom_mass)
      if L_cn_tmp:
        L_cn = L_cn_tmp
        max_cn_axis = 2
        main_axis = init_axis
        main_plane = init_plane
        cad_plane = 'sigma_'+main_plane
        dxm.append(uvw) 
        uvw_asym = uvw
  
  # --- C2 axis perpendicular to the main axis with redundancies
  # ... Check for vertical planes too
  # --- If not C2 then looks for a Cs plane

  #
  # ---finds C2 axis 
  #
  if L_cn:
    if L_prolate or L_oblate: uvw_asym = uvw 
    new_c2,dxc2 = get_c2(dict_sea,mat_geo,num_set,uvw_asym,L_sphera,epsilon,atom_mass,dxm)
    new_c2_p = 0
    if new_c2 > 0:
      for i in range(new_c2):
        cnrot = Cngen(2,dxc2[i])
        L_kk=get_sym(dict_sea,mat_geo,cnrot,num_set,'C2p',tolsym,atom_mass)
        if L_kk: new_c2_p += 1
    if new_c2_p >= max_cn_axis : L_c2_p = True
  #
  # ---finds vertical planes sigma_v
  #
    coor_sea = []
    mass_sea = []
    num_sigma_v = 0
    for l in range(num_set):
     idx = l
     isea,mass_sea,coor_sea =  get_sea(dict_sea,mat_geo,idx,atom_mass)
     if isea > 1:
      for i in range(isea):
        for j in range(i+1,isea):
          dv,L_nosig = dist_v(coor_sea[i],coor_sea[j])
          if L_nosig:
            continue
          else: 
            sigma_v = reflex(dv,' ')
            L_sigma_tmp=get_sym(dict_sea,mat_geo,sigma_v,num_set,' ',tolsym,atom_mass)
            if(L_sigma_tmp): num_sigma_v += 1
    
  
    if num_sigma_v > 0: L_sigma_v = True
  #
  #--- look for horizontal reflection
  #
    sigma_h = np.matrix(reflex(udum,main_plane))
    cad_plane = 'sigma_h'
    L_sigma_h=get_sym(dict_sea,mat_geo,sigma_h,num_set,cad_plane,tolsym,atom_mass)
  #--- look for inversion
    L_i=get_sym(dict_sea,mat_geo,i_center,num_set,'inversion',tolsym,atom_mass)
  # ---look for S2n symmetry operation
    sn_axis = 2*max_cn_axis
    cnrot =  Cngen(sn_axis,uvw)
    snrot = sigma_h * cnrot
    cad_saxis = 's'+str(sn_axis)
    L_s2n=get_sym(dict_sea,mat_geo,snrot,num_set,cad_saxis,tolsym,atom_mass)
  else:
   #-- looking for a reflection plane
    for iaxis in range(3):
      init_axis = axis[iaxis]
      init_plane = plane[iaxis]
      cad_cs = 'Cs in axis ',axis[iaxis]
      csplane = reflex(udum,init_plane)
      L_cs=get_sym(dict_sea,mat_geo,csplane,num_set,cad_cs,tolsym,atom_mass)
      if L_cs:
        break
    uvw = [[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]
    num_sigma_cs = 0
    num_sigma_cs = get_sigma_d(dict_sea,mat_geo,num_set,uvw,3,tolsym,atom_mass)
    if num_sigma_cs > 0: L_cs = True
   #-- looking for inversion
    L_i=get_sym(dict_sea,mat_geo,i_center,num_set,'inversion',tolsym,atom_mass)

# Any Cn axis (n >=2)?
#    YES
#       Any C2 axis perpendicular to the Cn? 
#         YES  
#            Horizontal plane?
#               YES -> Dnh
#               NO
#                  Improper rotation?
#                     YES -> Dnd
#                     NO  -> Dn
#         NO 
#            Horizontal plane?
#               YES -> Cnh
#               NO  
#                  Vertical planes?
#                     YES -> Cnv
#                     NO
#                        Improper rotation?
#                           YES -> S2n
#                           NO  -> Cn
#    NO  
#      Any plane?
#         YES -> Cs
#            Inversion?
#               YES -> Ci
#               NO  -> C1

  if L_cn:
    if L_c2_p:
      if L_sigma_h:
        return 'D'+str(max_cn_axis)+'h',2*max_cn_axis
      else:
        if L_s2n:
          return 'D'+str(max_cn_axis)+'d',2*max_cn_axis
        else:
          return 'D'+str(max_cn_axis),2*max_cn_axis
    else:
      if L_sigma_h:
        return 'C'+str(max_cn_axis)+'h',max_cn_axis
      else: 
        if L_sigma_v:
          return 'C'+str(max_cn_axis)+'v',max_cn_axis
        else:
          if L_s2n:
            return 'S'+str(sn_axis),max_cn_axis
          else:
            return 'C'+str(max_cn_axis),max_cn_axis
  else:
    if L_cs: 
      return 'Cs',1
    else:
      if L_i: 
        return 'Ci',1
      else:
        return 'C1',1

