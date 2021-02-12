'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: Pilgrim
Version     : 2021.1
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
| Sub-module :  dicts              |
| Last Update:  2020/02/03 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains different dictionaries.
They are, mainly, related to periodic table
''' 

#=============================================#
from numpy           import deg2rad
#---------------------------------------------#
from common.physcons import AMU, ANGSTROM
#=============================================#

#=============================================#
# PERIODIC TABLE: symbols, atnums, masses and #
#                 covalent radii              #
# dpt = dict for periodic table               #
#=============================================#

# Z ---> symbol
dpt_z2s = {  0:'XX' ,                                             \
             1:'H ' ,   2:'He' ,   3:'Li' ,   4:'Be' ,   5:'B ' , \
             6:'C ' ,   7:'N ' ,   8:'O ' ,   9:'F ' ,  10:'Ne' , \
            11:'Na' ,  12:'Mg' ,  13:'Al' ,  14:'Si' ,  15:'P ' , \
            16:'S ' ,  17:'Cl' ,  18:'Ar' ,  19:'K ' ,  20:'Ca' , \
            21:'Sc' ,  22:'Ti' ,  23:'V ' ,  24:'Cr' ,  25:'Mn' , \
            26:'Fe' ,  27:'Co' ,  28:'Ni' ,  29:'Cu' ,  30:'Zn' , \
            31:'Ga' ,  32:'Ge' ,  33:'As' ,  34:'Se' ,  35:'Br' , \
            36:'Kr' ,  37:'Rb' ,  38:'Sr' ,  39:'Y ' ,  40:'Zr' , \
            41:'Nb' ,  42:'Mo' ,  43:'Tc' ,  44:'Ru' ,  45:'Rh' , \
            46:'Pd' ,  47:'Ag' ,  48:'Cd' ,  49:'In' ,  50:'Sn' , \
            51:'Sb' ,  52:'Te' ,  53:'I ' ,  54:'Xe' ,  55:'Cs' , \
            56:'Ba' ,  57:'La' ,  58:'Ce' ,  59:'Pr' ,  60:'Nd' , \
            61:'Pm' ,  62:'Sm' ,  63:'Eu' ,  64:'Gd' ,  65:'Tb' , \
            66:'Dy' ,  67:'Ho' ,  68:'Er' ,  69:'Tm' ,  70:'Yb' , \
            71:'Lu' ,  72:'Hf' ,  73:'Ta' ,  74:'W ' ,  75:'Re' , \
            76:'Os' ,  77:'Ir' ,  78:'Pt' ,  79:'Au' ,  80:'Hg' , \
            81:'Tl' ,  82:'Pb' ,  83:'Bi' ,  84:'Po' ,  85:'At' , \
            86:'Rn' ,  87:'Fr' ,  88:'Ra' ,  89:'Ac' ,  90:'Th' , \
            91:'Pa' ,  92:'U ' ,  93:'Np' ,  94:'Pu' ,  95:'Am' , \
            96:'Cm' ,  97:'Bk' ,  98:'Cf' ,  99:'Es' , 100:'Fm' , \
           101:'Md' , 102:'No' , 103:'Lr'                         \
          }

# Z --> mass
dpt_z2m = {  0:  0.000000 ,                                                    \
             1:  1.007825 ,   2:  4.002600 ,   3:  7.016000 ,   4:  9.012180 , \
             5: 11.009310 ,   6: 12.000000 ,   7: 14.003070 ,   8: 15.994910 , \
             9: 18.998400 ,  10: 19.992440 ,  11: 22.989800 ,  12: 23.985040 , \
            13: 26.981530 ,  14: 27.976930 ,  15: 30.973760 ,  16: 31.972070 , \
            17: 34.968850 ,  18: 39.948000 ,  19: 38.963710 ,  20: 39.962590 , \
            21: 44.955920 ,  22: 47.900000 ,  23: 50.944000 ,  24: 51.940500 , \
            25: 54.938100 ,  26: 55.934900 ,  27: 58.933200 ,  28: 57.935300 , \
            29: 62.929800 ,  30: 63.929100 ,  31: 68.925700 ,  32: 73.921900 , \
            33: 74.921600 ,  34: 79.916500 ,  35: 78.918300 ,  36: 83.911500 , \
            37: 84.911700 ,  38: 87.905600 ,  39: 89.905400 ,  40: 89.904300 , \
            41: 92.906000 ,  42: 97.905500 ,  43: 97.000000 ,  44:101.903700 , \
            45:102.904800 ,  46:105.903200 ,  47:106.904100 ,  48:113.903600 , \
            49:114.904100 ,  50:119.902200 ,  51:120.903800 ,  52:129.906700 , \
            53:126.904400 ,  54:131.904200 ,  55:132.905400 ,  56:137.905200 , \
            57:138.906300 ,  58:139.905400 ,  59:140.907600 ,  60:141.907700 , \
            61:144.912700 ,  62:151.919700 ,  63:152.921200 ,  64:157.924100 , \
            65:158.925300 ,  66:163.929200 ,  67:164.930300 ,  68:165.930300 , \
            69:168.934200 ,  70:173.938900 ,  71:174.940800 ,  72:179.946500 , \
            73:180.948000 ,  74:183.950900 ,  75:186.955700 ,  76:191.961500 , \
            77:192.962900 ,  78:194.964800 ,  79:196.966500 ,  80:201.970600 , \
            81:204.974400 ,  82:207.976600 ,  83:208.980400 ,  84:208.982400 , \
            85:209.987100 ,  86:222.017600 ,  87:223.019700 ,  88:226.025400 , \
            89:227.027800 ,  90:232.038100 ,  91:231.035900 ,  92:238.050800 , \
            93:237.048200 ,  94:244.064200 ,  95:243.061400 ,  96:247.070300 , \
            97:247.070300 ,  98:251.079600 ,  99:252.082900 , 100:257.075100 , \
           101:258.098600 , 102:259.100900 , 103:260.105300 }

# Z --> covalent radius
dpt_z2cr = {  0:-10.00 ,                                             \
              1:0.31 ,   2:0.28 ,   3:1.28 ,   4:0.96 ,   5:0.84 , \
              6:0.73 ,   7:0.71 ,   8:0.66 ,   9:0.57 ,  10:0.58 , \
             11:1.66 ,  12:1.41 ,  13:1.21 ,  14:1.11 ,  15:1.07 , \
             16:1.05 ,  17:1.02 ,  18:1.06 ,  19:2.03 ,  20:1.76 , \
             21:1.70 ,  22:1.60 ,  23:1.53 ,  24:1.39 ,  25:1.50 , \
             26:1.42 ,  27:1.38 ,  28:1.24 ,  29:1.32 ,  30:1.22 , \
             31:1.22 ,  32:1.20 ,  33:1.19 ,  34:1.20 ,  35:1.20 , \
             36:1.16 ,  37:2.20 ,  38:1.95 ,  39:1.90 ,  40:1.75 , \
             41:1.64 ,  42:1.54 ,  43:1.47 ,  44:1.46 ,  45:1.42 , \
             46:1.39 ,  47:1.45 ,  48:1.44 ,  49:1.42 ,  50:1.39 , \
             51:1.39 ,  52:1.38 ,  53:1.39 ,  54:1.40 ,  55:2.44 , \
             56:2.15 ,  57:2.07 ,  58:2.04 ,  59:2.03 ,  60:2.01 , \
             61:1.99 ,  62:1.98 ,  63:1.98 ,  64:1.96 ,  65:1.94 , \
             66:1.92 ,  67:1.92 ,  68:1.89 ,  69:1.90 ,  70:1.87 , \
             71:1.87 ,  72:1.75 ,  73:1.70 ,  74:1.62 ,  75:1.51 , \
             76:1.44 ,  77:1.41 ,  78:1.36 ,  79:1.36 ,  80:1.32 , \
             81:1.45 ,  82:1.46 ,  83:1.48 ,  84:1.40 ,  85:1.50 , \
             86:1.50 ,  87:2.60 ,  88:2.21 ,  89:2.15 ,  90:2.06 , \
             91:2.00 ,  92:1.96 ,  93:1.90 ,  94:1.87 ,  95:1.80 , \
             96:1.69 }

# dict for isotopic elements
dpt_im = {'C13 ': 13.00335500 , 'C14 ': 14.00324200 , 'Cl37': 36.96590300 , \
          'D   ':  2.01410178 , 'T   ':  3.01604928 , 'N15 ': 15.00010900 , \
          'O17 ': 16.99913200 , 'O18 ': 17.99916000 }

# symbol --> Z
dpt_s2z  = dict( (s.strip(), z        ) for z,s in dpt_z2s.items())

# symbol --> mass
dpt_s2m  = dict( (dpt_z2s[z].strip(),m) for z,m in dpt_z2m.items() )

# symbol --> covalent radius
dpt_s2cr = dict( (dpt_z2s[z].strip(),r) for z,r in dpt_z2cr.items() )

# add also lower case
symbols = list(dpt_s2z.keys())
for key in symbols:
    if key in dpt_s2z.keys() : dpt_s2z[key.lower()]  = dpt_s2z[key]
    if key in dpt_s2m.keys() : dpt_s2m[key.lower()]  = dpt_s2m[key]
    if key in dpt_s2cr.keys(): dpt_s2cr[key.lower()] = dpt_s2cr[key]

#=====================================#
# Dict of points for 5,6-member rings #
#=====================================#

# 5RING
dicCONF5 = {}

# 6RING; each point = (theta,phi); theta in [0,pi], phi in [0,2pi]
dicCONF6_v2 = {'C1' :(  0.0,  0.0), \
               'C2' :(180.0,  0.0), \

               'B1 ':( 90.0,  0.0), \
               'TB1':( 90.0, 90.0), \
               'B2 ':( 90.0,180.0), \
               'TB2':( 90.0,270.0), \

               'HB1':( 45.0,  0.0), \
               'HC1':( 45.0, 90.0), \
               'HB2':( 45.0,180.0), \
               'HC2':( 45.0,270.0), \

               'HB3':(135.0,  0.0), \
               'HC3':(135.0, 90.0), \
               'HB4':(135.0,180.0), \
               'HC4':(135.0,270.0), \
                }


# 6RING; each point = (theta,phi); theta in [0,pi], phi in [0,2pi]
dicCONF6 = {'B25':( 90.0, 60.0) , 'B36':( 90.0,300.0) , 'B41':( 90.0,180.0) , \
            'E1 ':(125.3,180.0) , 'E2 ':( 54.7, 60.0) , 'E3 ':(125.3,300.0) , \
            'E4 ':( 54.7,180.0) , 'E5 ':(125.3, 60.0) , 'E6 ':( 54.7,300.0) , \
            '1C4':(  0.0,  0.0) , '1E ':( 54.7,  0.0) , '1H2':( 50.8, 30.0) , \
            '1H6':( 50.8,330.0) , '1S2':( 67.5, 30.0) , '1S6':( 67.5,330.0) , \
            '1T3':( 90.0,330.0) , '2E ':(125.3,240.0) , '2H1':(129.2,210.0) , \
            '2H3':(129.2,270.0) , '2S1':(112.5,210.0) , '2S3':(112.5,270.0) , \
            '2T4':( 90.0,210.0) , '2T6':( 90.0,270.0) , '3E ':( 54.7,120.0) , \
            '3H2':( 50.8, 90.0) , '3H4':( 50.8,150.0) , '3S2':( 67.5, 90.0) , \
            '3S4':( 67.5,150.0) , '3T1':( 90.0,150.0) , '4C1':(180.0,  0.0) , \
            '4E ':(125.3,  0.0) , '4H3':(129.2,330.0) , '4H5':(129.2, 30.0) , \
            '4S3':(112.5,330.0) , '4S5':(112.5, 30.0) , '4T2':( 90.0, 30.0) , \
            '5E ':( 54.7,240.0) , '5H4':( 50.8,210.0) , '5H6':( 50.8,270.0) , \
            '5S4':( 67.5,210.0) , '5S6':( 67.5,270.0) , '6E ':(125.3,120.0) , \
            '6H1':(129.2,150.0) , '6H5':(129.2, 90.0) , '6S1':(112.5,150.0) , \
            '6S5':(112.5, 90.0) , '6T2':( 90.0, 90.0) , '14B':( 90.0,  0.0) , \
            '25B':( 90.0,240.0) , '36B':( 90.0,120.0) }

#=============================#
# Preparation of dictionaries #
#=============================#
dpt_im   = {key.strip():value for key,value in dpt_im.items()}
dicCONF6 = {key.strip():value for key,value in dicCONF6.items() }
dicCONF6_v2 = {key.strip():value for key,value in dicCONF6_v2.items() }

#====================#
# print dictionaries #
#====================#
def string_dict(thedict,name="thedict",f1="'%2s'",f2="%3i",each=10):
    '''
    a function just to print a dictionary
    '''
    cc = 0
    string = "%s = {"%name
    for key in sorted(thedict.keys()):
        string += f1%key + ":" + f2%thedict[key] + " , "
        cc += 1
        if cc % each == 0:
           string += "\\\n" + " "*len(name) + "    "
    string =  string[:-2]+"}\n"
    return string

def print_dicts():
    print(string_dict(dpt_z2s ,name="dpt_z2s" ,f1="%3i"   ,f2="'%2s'"        ,each=5))
    print(string_dict(dpt_z2m ,name="dpt_z2m" ,f1="%3i"   ,f2="%10.6f"       ,each=4))
    print(string_dict(dpt_z2cr,name="dpt_z2cr",f1="%3i"   ,f2="%4.2f"        ,each=5))
    print(string_dict(dpt_im  ,name="dpt_im"  ,f1="'%-4s'",f2="%12.8f"       ,each=3))
    print(string_dict(dicCONF6,name="dicCONF6",f1="'%-7s'",f2="(%5.1f,%5.1f)",each=3))
    print(string_dict(dpt_s2z ,name="dpt_s2z" ,f1="'%2s'" ,f2="%3i"          ,each=5))
    print(string_dict(dpt_s2m ,name="dpt_s2m" ,f1="'%2s'" ,f2="%10.6f"       ,each=4))
    print(string_dict(dpt_s2cr,name="dpt_s2cr",f1="'%2s'" ,f2="%4.2f"        ,each=5))

if __name__ == '__main__': print_dicts()

#=======================================#
# Unit conversion (to au, radians, etc) #
#=======================================#
dpt_z2m  = {k:v/AMU           for k,v       in dpt_z2m.items() }
dpt_s2m  = {k:v/AMU           for k,v       in dpt_s2m.items() }
dpt_im   = {k:v/AMU           for k,v       in dpt_im.items()  }
dpt_z2cr = {k:v/ANGSTROM      for k,v       in dpt_z2cr.items()}
dpt_s2cr = {k:v/ANGSTROM      for k,v       in dpt_s2cr.items()}
dicCONF6 = {k:(deg2rad(v1),deg2rad(v2)) for k,(v1,v2) in dicCONF6.items()   }

