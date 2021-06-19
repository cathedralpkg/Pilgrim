'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: Pilgrim
Version     : 2021.4
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
| Module     :  modpilgrim         |
| Sub-module :  optSUMMARY         |
| Last Update:  2021/04/20 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

Module for the --summary option
of Pilgrim
'''

#==================================================#
import random
import os
import numpy as np
#--------------------------------------------------#
from   common.physcons    import KCALMOL,KB
from   common.fncs        import exp128
#--------------------------------------------------#
from   modpilgrim.diverse import status_check
from   modpilgrim.diverse import ffchecking
from   modpilgrim.names   import get_pof
#==================================================#


#==================================================#
#           RELATED TO PFN OUTPUT FILES            #
#==================================================#
def readout_pfn(pof):
    with open(pof,'r') as asdf: lines = asdf.readlines()
    data = {}
    data["all"] = {}
    for idx,line in enumerate(lines):
        # weight for each itc
        if "| weight |" in line:
            idx2  = idx+2
            while True:
                if "------" in lines[idx2]: break
                itc,relV0,relV1,zpe,mass,weight = lines[idx2].replace("|"," ").split()[0:6]
                data[itc] = {}
                data[itc]["weight"] = weight
                data[itc]["relV0" ] = relV0
                data[itc]["relV1" ] = relV1
                idx2 += 1
        # Conformer
        if "Conformation:" in line:
            itc = line.split()[-1]
            #data[itc] = {}
        # V0, V1 and ZPE
        if "Electronic energy " in line:
            data[itc]["V0" ] = line.split()[-2]
            data[itc]["V1" ] = line.split()[-2]
        if "V0 + zero-point en" in line: data[itc]["V1" ] = line.split()[-2]
        if "zero-point energy:" in line: data[itc]["zpe"] = line.split()[-3]
        # part fncs
        if "| Partition functions (pfns):" in line:
            data[itc]["T"   ] = []
            data[itc]["Qtr" ] = []
            data[itc]["Qrot"] = []
            data[itc]["Qvib(V1)"] = []
            data[itc]["Qel" ] = []
            data[itc]["Qtot(V1)"] = []
            idx2  = idx+4
            while True:
                if "------" in lines[idx2]: break
                T,Qtr,Qrot,Qvib,Qel,Qtot = lines[idx2].replace("|"," ").split()
                data[itc]["T"  ].append(T)
                data[itc]["Qtr"  ].append(Qtr )
                data[itc]["Qrot" ].append(Qrot)
                data[itc]["Qvib(V1)" ].append(Qvib)
                data[itc]["Qel"  ].append(Qel )
                data[itc]["Qtot(V1)" ].append(Qtot)
                idx2 += 1
        # Gibbs
        if "| Gibbs free energy (hartree):" in line:
            data[itc]["T"   ] = []
            data[itc]["G(V)"] = []
            data[itc]["G(p)"] = []
            idx2  = idx+4
            while True:
                if "------" in lines[idx2]: break
                T,GV0,Gp0 = lines[idx2].replace("|"," ").split()
                GV0 = "%.3f"%((float(GV0)-float(data[itc]["V0"]))*KCALMOL)
                Gp0 = "%.3f"%((float(Gp0)-float(data[itc]["V0"]))*KCALMOL)
                data[itc]["T"   ].append(T)
                data[itc]["G(V)"].append(GV0)
                data[itc]["G(p)"].append(Gp0)
                idx2 += 1
        # min(V0), min(V1)
        if "min(V0) =" in line: data["all"]["V0" ] = line.split()[-2]
        if "min(V1) =" in line: data["all"]["V1" ] = line.split()[-2]
        # total part fncs
        if "Total multi-structural HO" in line:
            data["all"]["T"   ] = []
            data["all"]["Qtot(V1)" ] = []
            data["all"]["G(V)"] = []
            data["all"]["G(p)"] = []
            idx2  = idx+4
            while True:
                if "------" in lines[idx2]: break
                T,QMSHO,GV0,Gp0 = lines[idx2].replace("|"," ").split()
                GV0 = "%.3f"%((float(GV0)-float(data["all"]["V0"]))*KCALMOL)
                Gp0 = "%.3f"%((float(Gp0)-float(data["all"]["V0"]))*KCALMOL)
                data["all"]["T"   ].append(T)
                data["all"]["Qtot(V1)" ].append(QMSHO)
                data["all"]["G(V)"].append(GV0)
                data["all"]["G(p)"].append(Gp0)
                idx2 += 1
        # Anharmonicity
        if " ZPE_MSHO:" in line: data["all"]["ZPE_MSHO"] = line.split()[-2]
        if " ZPE_ANH :" in line: data["all"]["ZPE_ANH" ] = line.split()[-2]
        if "Calculating anh. ratio" in line:
            data["all"]["ANH. RATIO"   ] = []
            idx2  = idx+5
            while True:
                if "------" in lines[idx2]: break
                T,anh = lines[idx2].replace("|"," ").split()
                data["all"]["ANH. RATIO"   ].append(anh)
                idx2 += 1
    # Rovibrational itc
    for itc in data.keys():
        if itc == "all": continue
        Qrot = np.array([float(v) for v in data[itc]["Qrot"]])
        Qvib = np.array([float(v) for v in data[itc]["Qvib(V1)"]])
        Qrv  = ["%.3E"%v for v in Qrot*Qvib]
        data[itc]["Qrv(V1)" ] = Qrv
    # Rovibrational all
    try:
        Qtot = np.array([float(v) for v in data["all"]["Qtot(V1)"]])
        Qtr  = np.array([float(v) for v in data[itc]["Qtr"]])
        Qel  = np.array([float(v) for v in data[itc]["Qel"]])
        Qrv  = ["%.3E"%v for v in Qtot/(Qtr*Qel)]
        data["all"]["Qrv(V1)"] = Qrv
        data["all"]["Qtr"] = list(data[itc]["Qtr"])
        data["all"]["Qel"] = list(data[itc]["Qel"])
    except: pass
    # vib with regards to V0
    for itc in data.keys():
        try:
           # correction factor from V1 to V0
           V0, V1 = float(data[itc]["V0"]), float(data[itc]["V1"])
           V1toV0 = np.array([exp128((V0-V1)/KB/float(T)) \
                              for T in data[itc]["T"]])
           # correct Qvib
           if itc != "all":
              Qvib = np.array([float(v) for v in data[itc]["Qvib(V1)"]])
              data[itc]["Qvib"] = ["%.3E"%v for v in Qvib * V1toV0]
           # correct Qrv
           Qrv  = np.array([float(v) for v in data[itc]["Qrv(V1)"]])
           data[itc]["Qrv" ] = ["%.3E"%v for v in Qrv * V1toV0]
           # correct Qtot
           Qtot = np.array([float(v) for v in data[itc]["Qtot(V1)"]])
           data[itc]["Qtot" ] = ["%.3E"%v for v in Qtot * V1toV0]
        except: pass
    # ANH partition function
    try:
       Qar  = np.array([float(v) for v in data["all"]["ANH. RATIO"]])
       Qtot = np.array([float(v) for v in data["all"]["Qtot"]])
       T    = np.array([float(v) for v in data["all"]["T"]])
       Gcor = - (np.log(Qar)*T*KB)*KCALMOL
       GV   = np.array([float(v) for v in data["all"]["G(V)"]])
       Gp   = np.array([float(v) for v in data["all"]["G(p)"]])
       data["all"]["Qanh"   ] = ["%.3E"%v for v in Qar*Qtot]
       data["all"]["Ganh(V)"] = ["%.3f"%v for v in GV+Gcor]
       data["all"]["Ganh(p)"] = ["%.3f"%v for v in Gp+Gcor]
    except: pass
    # return all
    return data
#--------------------------------------------------#
def genpfntable1(data,itcs):
    props = "conf,weight,V0,V1,relV0,relV1"
    frow  = " %7s | %6s | %14s | %14s | %6s | %6s "
    head  = frow%tuple(props.split(","))

    string  = "="*len(head)+"\n"
    string += "Columns:\n"
    string += "     - conf     : conformer\n"
    string += "     - weight   : conformer weight\n"
    string += "     - V0       : electronic energy [hartree]\n"
    string += "     - V1       : V0+zero-point energy [hartree]\n"
    string += "     - relV0    : relative V0 [kcal/mol]\n"
    string += "     - relV1    : relative V1 [kcal/mol]\n"
    string += "\n"

    string += "-"*len(head)+"\n" + head +"\n" + "-"*len(head)+"\n"
    for itc in itcs:
        data_i = data[itc]
        vals   = [itc]
        for prop in props.split(",")[1:]:
            vals.append(data_i.get(prop," - "))
        string += frow%tuple(vals)+"\n"
    string += "-"*len(head)+"\n\n"
    string += "="*len(head)+"\n\n"
    # print
    for line in string.split("\n"): print("     %s"%line)
#--------------------------------------------------#
def genpfntable2(data,itcs,stemps):
    if len(itcs) == 1: case, props, frow = 1, "T (K),"," %7s |"
    else             : case, props, frow = 2, "conf," ," %7s |"
    if case == 2: itcs = itcs + ["all"]
    props += "Qtr,Qel,Qrot,Qvib,Qrv,Qtot,Qanh"
    frow  += " %10s | %10s | %10s | %10s | %10s | %10s | %10s "
    head = frow%tuple(props.split(","))
    sinfo  = "="*len(head)+"\n"
    sinfo += "Columns (pf stands for partition function):\n"
    sinfo += "     - Qtr  : traslational     pf (per unit volume)\n"
    sinfo += "     - Qel  : electronic       pf\n"
    sinfo += "     - Qrot : rotational       pf\n"
    sinfo += "     - Qvib : vibrational      pf\n"
    sinfo += "     - Qrv  : rovibrational    pf\n"
    sinfo += "     - Qtot : total MS-HO      pf (per unit volume)\n"
    sinfo += "     - Qanh : total anharmonic pf (per unit volume)\n"
    sinfo += "\n"
    sinfo += "     * Qvib, Qrv, Qtot and Qanh --> relative to V0\n"
    sinfo += "\n"
    for line in sinfo.split("\n"): print("     %s"%line)

    if case == 1:
       string = "-"*len(head)+"\n" + head +"\n" + "-"*len(head)+"\n"
    else: string = ""
    for sT in stemps:
        if case == 2:
           string += "TEMPERATURE: %s K\n"%sT
           string += "-"*len(head)+"\n" + head +"\n" + "-"*len(head)+"\n"
        else:
           vals = [sT]
        for itc in itcs:
            data_i = data[itc]
            if case == 2: vals = [itc]
            if sT not in data_i.get("T",[]): continue
            idxT = data_i["T"].index(sT)
            for prop in props.split(",")[1:]:
                try   : val = data_i[prop]; bool1 =True
                except: val = " - "
                # temp dependency
                if type(val) in[list,tuple]: val = val[idxT]
                vals.append(val)
            string += frow%tuple(vals)+"\n"
        if case == 2: string += "-"*len(head)+"\n\n\n"
    if case == 1: string += "-"*len(head)+"\n"
    string += "="*len(head)+"\n\n"
    # print
    for line in string.split("\n"): print("     %s"%line)
    print("")
#--------------------------------------------------#
def genpfntable3(data,itcs,stemps):
    if len(itcs) == 1: case, props, frow = 1, "T (K),"," %7s |"
    else             : case, props, frow = 2, "conf," ," %7s |"
    if case == 2: itcs = itcs + ["all"]
    props += "G(V),G(p),Ganh(V),Ganh(p)"
    frow  += " %10s | %10s | %10s | %10s "
    head = frow%tuple(props.split(","))
    sinfo  = "="*len(head)+"\n"
    sinfo += "Columns:\n"
    sinfo += "     - G    : MS-HO Gibbs free energy [kcal/mol]\n"
    sinfo += "     - Ganh : anharmonic Gibbs free energy [kcal/mol]\n"
    sinfo += "\n"
    sinfo += "     * Values are relative to V0\n"
    sinfo += "     * (V)  --> for a volume per molecule of V = 1cm^3\n"
    sinfo += "     * (p)  --> for a volume per molecule of V = kB*T/p0\n"
    sinfo += "                with p0=1 bar\n"
    sinfo += "\n"
    for line in sinfo.split("\n"): print("     %s"%line)

    if case == 1:
       string = "-"*len(head)+"\n" + head +"\n" + "-"*len(head)+"\n"
    else: string = ""
    for sT in stemps:
        if case == 2:
           string += "TEMPERATURE: %s K\n"%sT
           string += "-"*len(head)+"\n" + head +"\n" + "-"*len(head)+"\n"
        else:
           vals = [sT]
        for itc in itcs:
            data_i = data[itc]
            if case == 2: vals = [itc]
            if sT not in data_i.get("T",[]): continue
            idxT = data_i["T"].index(sT)
            for prop in props.split(",")[1:]:
                try   : val = data_i[prop]; bool1 =True
                except: val = " - "
                # temp dependency
                if type(val) in[list,tuple]: val = val[idxT]
                vals.append(val)
            string += frow%tuple(vals)+"\n"
        if case == 2: string += "-"*len(head)+"\n\n\n"
    if case == 1: string += "-"*len(head)+"\n"
    # print
    string += "="*len(head)+"\n\n"
    for line in string.split("\n"): print("     %s"%line)
    print("")
#--------------------------------------------------#
def summary_pfn(targets,ltemp,dctc,dlevel):
    # Temperature(s) & CTCs
    ctcs   = []
    stemps = []
    for target in targets:
        try   : stemps.append(float(target))
        except: ctcs.append(target)
    # selected system
    if len(ctcs  ) == 0:
       ctcs = list(dctc.keys())
       if True:
          print("    --> System NOT selected  :")
          ml = max([len(ctc) for ctc in ctcs])
          for idx in range(0,len(ctcs),4):
              print("        "+"  ".join("%%-%is"%ml%ctc for ctc in ctcs[idx:idx+4]))
          return
    # Only 1 system!
    try   : ctc = random.choice(ctcs)
    except: return
    # Selected temperature(s)
    if len(stemps) == 0: stemps = list(ltemp)
    # float --> str
    stemps = ["%.2f"%T for T in sorted(stemps)]
    print("    --> Molecule             : %s"%ctc)
    if len(stemps) == 1:
       print("    --> Selected temperature : %s K"%stemps[0])
    print("")

    #-------------------#
    # Read output files #
    #-------------------#
    # whole ctc or individual itc?
    if "." in ctc: ctc,the_itc = ctc.split(".")
    else: the_itc = None
    # read ctc file
    pof = get_pof(dlevel,"pfn",ctc)
    print("    Reading file: %s"%pof)
    if not os.path.exists(pof):
       print("    --> NOT FOUND")
       return
    data = readout_pfn(pof)
    print("")
    # list of conformers
    if the_itc is None: itcs = list([itc for itc in data.keys() if itc != "all"])
    else              : itcs = [the_itc]

    #-----------------#
    # Generate tables #
    #-----------------#
    genpfntable1(data,itcs)
    genpfntable2(data,itcs,stemps)
    genpfntable3(data,itcs,stemps)
#==================================================#

#==================================================#
#           RELATED TO PATH OUTPUT FILES           #
#==================================================#
def readout_path(pof,mep=False):
    with open(pof,'r') as asdf: lines = asdf.readlines()
    for idx,line in enumerate(lines):
        # initialize variables
        if "Variables for first step" in line:
           data = {}
           continue
        # Imaginary frequency
        if "Vibrational frequencies [" in line:
            data["ifreq"] = lines[idx+1].split()[0].replace("-","")+"i"
        # Imag freqs along MEP
        if "Fine! There are no ima" in line: data["iMEP"] = False
        if "WARNING! There are ima" in line: data["iMEP"] = True
        # CVT
        if "s_CVT" in line and "Gamma^CVT" in line:
            data["T"] = []
            data["s_CVT"] = []
            data["Gamma^CVT"] = []
            idx2  = idx+2
            while True:
                if "------" in lines[idx2]: break
                T, s, gamma = lines[idx2].replace("|"," ").split()[0:3]
                data["T"].append(T)
                data["s_CVT"].append(s)
                data["Gamma^CVT"].append(gamma)
                idx2 += 1
        # Max of VaG
        if "Maximum of VaG (VAG)" in line:
            data["sAG"] = line.split()[ 0]
            data["VAG"] = line.split()[-3]
        # E0 and (again) VAG
        if "E0  =" in line: data["E0" ] = line.split()[-2]
        if "VAG =" in line: data["VAG"] = line.split()[-2]
        # ZCT
        if "ZCT transmission coefficient:" in line:
            data["T"  ] = []
            data["kappa^ZCT"  ] = []
            data["RTE(ZCT)"] = []
            idx2  = idx+4
            while True:
                if "------" in lines[idx2]: break
                T, i1, i2, kappa, rte = lines[idx2].replace("|"," ").replace("*","").replace("++","").split()
                data["T"  ].append(T)
                data["kappa^ZCT"  ].append(kappa)
                data["RTE(ZCT)"].append(rte)
                idx2 += 1
        # SCT
        if "SCT transmission coefficient:" in line:
            data["T"  ] = []
            data["kappa^SCT"  ] = []
            data["RTE(SCT)"] = []
            idx2  = idx+4
            while True:
                if "------" in lines[idx2]: break
                T, i1, i2, kappa, rte = lines[idx2].replace("|"," ").replace("*","").replace("++","").split()
                data["T"  ].append(T)
                data["kappa^SCT"  ].append(kappa)
                data["RTE(SCT)"].append(rte)
                idx2 += 1
        # CAG
        if "Calculating CAG coefficient..." in line:
            data["T"  ] = []
            data["TST/CAG"] = []
            data["CVT/CAG"] = []
            idx2  = idx+5
            while True:
                if "------" in lines[idx2]: break
                T, d1, cagtst, d2, cagcvt = lines[idx2].replace("|"," ").split()
                data["T"  ].append(T)
                data["TST/CAG"].append(cagtst)
                data["CVT/CAG"].append(cagcvt)
                idx2 += 1
        # P(E0)
        if "P^ZCT(E)" in line:
            tmp   = lines[idx+3].replace("|"," ")
            tmp   = tmp.replace("B","").replace("L","").replace("R","")
            s1s2  = tmp.split()[3].replace("[","").replace("]","")
            data["E0"]        = tmp.split()[0]
            data["P^ZCT(E0)"] = tmp.split()[1]
            data["P^SCT(E0)"] = tmp.split()[2]
            data["s<"]        = s1s2.split(",")[0]
            data["s>"]        = s1s2.split(",")[1]
        # Get data along MEP
        if mep and "Eref + VaG (au)" in line:
            data["s"      ] = []
            data["V_MEP"  ] = []
            data["VaG(cc)"] = []
            data["VaG(ic)"] = []
            idx2  = idx+2
            while True:
                if "------" in lines[idx2]: break
                s,Vmep,zpecc,zpeic = lines[idx2].replace("|"," ").split()[0:4]
                data["s"      ].append(s)
                data["V_MEP"  ].append(Vmep)
                Vadicc = "%.3f"%(float(Vmep)+float(zpecc))
                try   : Vadiic = "%.3f"%(float(Vmep)+float(zpeic))
                except: Vadiic = " - "
                data["VaG(cc)"].append(Vadicc)
                data["VaG(ic)"].append(Vadiic)
                idx2 += 1
        if mep and "- Effective mass (mueff) in a.u." in line:
            data["mueff/mu"] = []
            idx2  = idx+4
            while True:
                if "------" in lines[idx2]: break
                mueff = lines[idx2].replace("|"," ").replace("*","").split()[-1]
                data["mueff/mu"].append(mueff)
                idx2 += 1
    # gamma^CVT/SCT
    try:
        p1 = np.array([float(v) for v in data["Gamma^CVT"]])
        p2 = np.array([float(v) for v in data["kappa^SCT"]])
        p3 = np.array([float(v) for v in data["CVT/CAG"]])
        data["gamma^CVT/SCT"] = ["%.3E"%v for v in p1*p2*p3]
    except:
        pass
    # return all
    return data
#--------------------------------------------------#
def genpathtable1(data):
    props = "conf,ifreq,E0,VAG,sAG,s<,s>,P^SCT(E0)"
    frow = " %5s | %9s | %10s | %10s | %8s | %8s | %8s | %10s "
    head = frow%tuple(props.split(","))
    string  = "="*len(head)+"\n"
    string += "Columns:\n"
    string += "     - conf     : conformer\n"
    string += "     - ifreq    : imaginary frequency for transition state [1/cm]\n"
    string += "     - E0       : lower energy limit for tunnelling [kcal/mol]\n"
    string += "     - VAG      : upper energy limit for tunnelling [kcal/mol]\n"
    string += "     - sAG      : MEP position corresponding to VAG [Bohr]\n"
    string += "     - s<       : MEP limit towards reactant(s) [Bohr]\n"
    string += "     - s>       : MEP limit towards product(s) [Bohr]\n"
    string += "     - P^SCT(E0): SCT tunnel probability at E0\n"
    string += "\n"
    string += "-"*len(head)+"\n"
    string +=         head +"\n"
    string += "-"*len(head)+"\n"
    any_imep = False
    for itc,data_i in data.items():
        vals = [itc]
        imep = data_i.get("iMEP",False)
        for prop in props.split(",")[1:]: vals.append( data_i.get(prop," - ") )
        string += frow%tuple(vals)
        if data_i.get("iMEP",False):
           any_imep = True
           string += "[*]"
        string += "\n"
    string += "-"*len(head)+"\n"
    if any_imep: string += "  [*] There are imaginary frequencies along this MEP\n"
    string += "\n"
    string += "="*len(head)+"\n"
    # print table
    for line in string.split("\n"): print("     %s"%line)
    print("")
#--------------------------------------------------#
def genpathtable2(data,stemps):
    if len(data.keys()) == 1: case, props = 1, "T (K),"
    else                    : case, props = 2, "conf,"
    props += "s_CVT,Gamma^CVT,kappa^SCT,gamma^CVT/SCT,RTE(SCT)"
    frow = " %7s | %9s | %11s | %11s | %13s | %8s "
    head = frow%tuple(props.split(","))
    sinfo  = "="*len(head)+"\n"
    sinfo += "Columns:\n"
    sinfo += "     - s_CVT        : MEP position of CVT generalized transition state [Borh]\n"
    sinfo += "     - Gamma^CVT    : CVT recrossing transmission coefficient\n"
    sinfo += "     - kappa^SCT    : SCT tunnelling transmission coefficient\n"
    sinfo += "     - gamma^CVT/SCT: CVT/SCT transmission coefficient\n"
    sinfo += "                      = Gamma^CVT * gamma^CVT/SCT * kappa^{CAG/CVT}\n"
    sinfo += "     - RTE(SCT)     : SCT representative tunnelling energy\n"
    sinfo += "\n"
    for line in sinfo.split("\n"): print("     %s"%line)

    itcs = list(data.keys())

    if case == 1:
       string = "-"*len(head)+"\n" + head +"\n" + "-"*len(head)+"\n"
    else: string = ""
    for sT in stemps:
        if case == 2:
           string += "TEMPERATURE: %s K\n"%sT
           string += "-"*len(head)+"\n" + head +"\n" + "-"*len(head)+"\n"
        else:
           vals = [sT]
        for itc in itcs:
            data_i = data[itc]
            if case == 2: vals = [itc]
            if sT not in data_i.get("T",[]): continue
            idxT = data_i["T"].index(sT)
            for prop in props.split(",")[1:]:
                try   : val = data_i[prop][idxT]; bool1 = True
                except: val = " - "
                vals.append(val)
            string += frow%tuple(vals)+"\n"
        if case == 2: string += "-"*len(head)+"\n\n\n"
    if case == 1: string += "-"*len(head)+"\n"
    # print
    for line in string.split("\n"): print("     %s"%line)
    print("")
    print("     "+"="*len(head))
    print("")
#--------------------------------------------------#
def genpathtable3(data):
    props  = "s,V_MEP,VaG(cc),VaG(ic),mueff/mu"
    frow   = " %7s | %11s | %11s | %11s | %8s "
    head = frow%tuple(props.split(","))
    sinfo  = "="*len(head)+"\n"
    sinfo += "Columns:\n"
    sinfo += "     - s        : MEP coordinate [Bohr]\n"
    sinfo += "     - V_MEP    : MEP total energy [kcal/mol]\n"
    sinfo += "     - VaG(cc)  : Adiabatic potential in Cartesian coordinates [kcal/mol] \n"
    sinfo += "     - VaG(ic)  : Adiabatic potential in internal coordinates [kcal/mol] \n"
    sinfo += "     - mueff/mu : ratio between effective mass and scalar mass\n"
    sinfo += "\n"
    for line in sinfo.split("\n"): print("     %s"%line)

    string = "-"*len(head)+"\n" + head +"\n" + "-"*len(head)+"\n"
    for idx,sval in enumerate(data["s"]):
        vals = [sval]
        for prop in props.split(",")[1:]:
            try   : val = data[prop][idx]
            except: val = " - "
            vals.append(val)
        string += frow%tuple(vals)+"\n"
    string += "-"*len(head)+"\n"
    string += "\n"
    string += "="*len(head)+"\n"
    # print
    for line in string.split("\n"): print("     %s"%line)
    print("")
#--------------------------------------------------#
def summary_path(targets,ltemp,dpath,dctc,dlevel):
    # Temperature(s) & transition states
    tss    = []
    stemps = []
    for target in targets:
        try   : stemps.append(float(target))
        except: tss.append(target)
    # Selected temperature(s)
    if len(stemps) == 0: stemps = list(ltemp)
    if len(stemps) == 0: print("    --> Temperature(s) NOT found!\n"); return
    stemps = ["%.2f"%T for T in sorted(stemps)]
    # Selected TS
    if   len(tss) == 0:
       tss = list(dpath.keys())
       if len(tss) == 0:
          print("    --> path files NOT found!")
          return
       if len(tss) != 1:
          print("    --> System NOT selected  :")
          ml = max([len(ts) for ts in tss])
          for idx in range(0,len(tss),4):
              print("        "+"  ".join("%%-%is"%ml%ts for ts in tss[idx:idx+4]))
          return
    if   len(tss) == 0:
       print("    --> Transition state(s) NOT found!\n"); return
    elif len(tss) != 1:
       ts = random.choice(tss)
       print("    --> Transition state     : %s [randomly selected]"%ts)
    else:
       ts = tss[0]
       print("    --> Transition state     : %s "%ts)
    # selected temperature
    if len(stemps) == 1:
       print("    --> Selected temperature : %s K"%stemps[0])
    else:
       print("    --> Number temperatures  : %i"%len(stemps))
    print("")

    #-------------------#
    # Read output files #
    #-------------------#
    ref_pof = get_pof(dlevel,"path",ts+".[xxx]")
    if ts in dctc.keys():
       clusterconf = dctc[ts]
       itcs = [itc for itc,weight in clusterconf._itcs]
       pofs = [ref_pof.replace("[xxx]",itc) for itc in itcs]
    else:
       folder  = "/".join(ref_pof.split("/")[0:-1])
       part1,part2 = ref_pof.split("/")[-1].split("[xxx]")
       pofs = [folder+"/"+pof for pof in os.listdir(folder) \
                              if part1 in pof and part2 in pof]
    # only available outputs
    pofs = [pof for pof in pofs if os.path.exists(pof)]
    pofs.sort()
    if len(pofs) == 0:
       print("    --> Unable to find any output file that")
       print("        matches with the selected arguments...\n")
       return

    data = {}
    print("    Reading output files (num files = %i):"%len(pofs))
    for pof in pofs:
        itc = pof.split(".")[-3]
        print("       --> %s"%pof)
        try:
           if len(pofs) == 1: data_itc = readout_path(pof,mep=True )
           else             : data_itc = readout_path(pof,mep=False)
        except:
           print("           problems reading file...")
           continue
        data[itc] = data_itc
    print("")

    #-----------------#
    # Generate tables #
    #-----------------# 
    genpathtable1(data)
    if len(pofs) == 1: genpathtable3(data[itc])
    genpathtable2(data,stemps)
#==================================================#


#==================================================#
#          RELATED TO RCONS OUTPUT FILES           #
#==================================================#
def kdir_from_lines(lines,rctype,which="k"):
    record = False
    btotal = False
    data   = {}
    for line in lines:
        line = line.strip()
        if line == ""       : continue
        if "|" not in line  : continue
        if "-------" in line: continue
        if "Current" in line: continue
        if "Elapsed" in line: continue
        # beginning of table
        if "(K)" in line:
            record = False
            # read MS-TST from MP table
            if rctype == "mstst" and \
              ("MS-CVT" in line or "MS-TST/ZCT" in line): continue
            # prepare
            cols   = line.replace("-","").replace("/","").lower()
            cols   = [col.strip() for col in cols.split("|")]
            btotal = "ts"   in cols
            record = rctype in cols
            # idx in table
            if record: col = cols.index(rctype)
            continue

        # get data
        if record:
           if btotal and "total" not in line: continue
           cols = line.split("|")
           temp = cols[0].strip()
           val  = cols[col].strip()
           if which.lower() == "k": data[temp] = val
           if which.lower() == "g": data[temp] = data.get(temp,[])+[val]
    return data
#--------------------------------------------------#
def readout_rcons(ofile,rctype):
    # Data to return
    V0dir, V1dir = None, None
    V0inv, V1inv = None, None
    kdir  = {}
    Gdir  = {}
    kinv  = {}
    Ginv  = {}
    # Read file
    with open(ofile,'r') as asdf: lines = list(asdf.readlines())
    # localize parts of importance
    idx1,idx2,idx3,idx4,idx5,idx6 = -1,-1,-1,-1,-1,-1
    idx7 = None
    nR   = 0
    for idx,line in enumerate(lines):
        if   "reactant(s)      ==>"    in line: nR   = len(line.split("==>")[-1].split("+"))
        if   "min{V1(i)} of reactants" in line: idx1 = idx
        elif "FORWARD  RATE CONSTANTS" in line: idx2 = idx
        elif "FORWARD  GIBBS FREE"     in line: idx3 = idx
        elif "BACKWARD RATE CONSTANTS" in line: idx4 = idx
        elif "BACKWARD GIBBS FREE"     in line: idx5 = idx
        elif "Updating plot file:"     in line: idx6 = idx
        elif "Reactants = Products"    in line: idx7 = idx
    # Read barrier
    if idx1 != None:
       count = 0
       for line in lines[idx1:]:
           if "SP: stationary point" in line: break
           if "-----" in line: count += 1; continue
           if count == 3:
              if V0dir is None: V0dir = float("inf")
              if V1dir is None: V1dir = float("inf")
              V0dir_i, V1dir_i = line.split("|")[1:3]
              V0dir = min(V0dir,float(V0dir_i))
              V1dir = min(V1dir,float(V1dir_i))
           if count == 4:
              if V0inv is None: V0inv = float("inf")
              if V1inv is None: V1inv = float("inf")
              V0inv_i, V1inv_i = line.split("|")[1:3]
              V0inv = min(V0inv,float(V0inv_i))
              V1inv = min(V1inv,float(V1inv_i))
    # Read rate constant
    if idx2 != -1:
       kdir = kdir_from_lines(lines[idx2:idx3],rctype,"k")
       Gdir = kdir_from_lines(lines[idx3:idx4],rctype,"g")
    if idx4 != -1:
       kinv = kdir_from_lines(lines[idx4:idx5],rctype,"k")
       Ginv = kdir_from_lines(lines[idx5:idx6],rctype,"g")
    # REACTANTS = PRODUCTS (duplicate rate constant and correct Gibbs with -RTln(2))
    if idx7 is not None:
       kdir = {T:"%.3E"%(2*float(k)) for T,k in kdir.items()}
       Gdir = {T:["%.3f"%(float(G1)-KB*float(T)*np.log(2)*KCALMOL), \
                  "%.3f"%(float(G2)-KB*float(T)*np.log(2)*KCALMOL)] \
                  for T,[G1,G2] in Gdir.items()}
       kinv = None
       Ginv = None
    # correct V0inv, V1inv
    if V0inv is not None: V0inv = V0dir - V0inv
    if V1inv is not None: V1inv = V1dir - V1inv
    return (V0dir,V1dir,kdir,Gdir),(V0inv,V1inv,kinv,Ginv)
#--------------------------------------------------#
def get_rconstr(rcname,stemp,dV0dir,dV1dir,dkdir,dGdir,dV0inv,dV1inv,dkinv,dGinv):
    try:
       V0dir = "%.2f"%dV0dir[rcname]
       V1dir = "%.2f"%dV1dir[rcname]
       kdir  = dkdir[rcname]
       Gdir  = dGdir[rcname]
    except:
       print("       --> Data NOT found for: %s"%rcname)
       raise Exception
    try   : kdir = kdir[stemp]
    except:
       print("       --> Temperature NOT found: %s K"%stemp)
       raise Exception
    try   : G1dir,G2dir = Gdir[stemp]
    except: G1dir,G2dir = " - ", " - "
    # inverse
    try:
       V0inv = "%.2f"%dV0inv[rcname]
       V1inv = "%.2f"%dV1inv[rcname]
       kinv        = dkinv[rcname][stemp]
       G1inv,G2inv = dGinv[rcname][stemp]
       bool_inv = True
    except:
       V0inv = " - "
       V1inv = " - "
       kinv = " - "
       G1inv,G2inv = " - ", " - "
       bool_inv = False
    # return data
    return (V0dir,V1dir,kdir,G1dir,G2dir), (V0inv,V1inv,kinv,G1inv,G2inv), bool_inv
#--------------------------------------------------#
def summary_rcons(targets,ltemp,dchem,dlevel):

    #-------------#
    # Preparation #
    #-------------#
    valid_rcons  = "mstst,mststzct,mststsct,mscvt,mscvtzct,mscvtsct,"
    valid_rcons +=       "mptstzct,mptstsct,mpcvt,mpcvtzct,mpcvtsct"
    valid_rcons  = valid_rcons.split(",")
    # Get selected rctype
    if len(targets) == 0:
       print("    --> Rate const REQUIRED  :")
       print("         %s  %s  %s  %s  %s  %s"%tuple(valid_rcons[0:6]))
       print("                %s  %s  %s  %s  %s"%tuple(valid_rcons[6: ]))
       return
    rctype = targets.pop(0).lower().strip()
    if rctype == "tst": rctype = "mstst"
    if rctype not in valid_rcons:
       print("    --> Invalid constant type: %s\n"%rctype)
       print("    --> VALID rate constants :")
       print("         %s  %s  %s  %s  %s  %s"%tuple(valid_rcons[0:6]))
       print("                %s  %s  %s  %s  %s"%tuple(valid_rcons[6: ]))
       return

    print("    --> Rate constant type   : %s"%rctype)

    # Temperature(s) & Reaction(s)
    srcnames = []
    stemps   = []
    for target in targets:
        try   : stemps.append(float(target))
        except: srcnames.append(target)
    if len(stemps  ) == 0: stemps   = list(ltemp)
    if len(srcnames) == 0: srcnames = list(dchem.keys())
    # float --> str
    stemps = ["%.2f"%T for T in sorted(stemps)]

    if len(stemps)   == 1:
       print("    --> Selected temperature : %s K"%stemps[0])
    if len(srcnames) == 1:
       print("    --> Selected reaction    : %s"%srcnames[0])
    print("")

    #-------------------#
    # Read output files #
    #-------------------#
    print("    Reading files:")
    dV0dir,dV1dir,dkdir,dGdir = {},{},{},{}
    dV0inv,dV1inv,dkinv,dGinv = {},{},{},{}
    count = 0
    for rcname in srcnames:
        pof = get_pof(dlevel,"rcons",rcname)
        if not os.path.exists(pof):
           print("       --> %s [NOT FOUND]"%pof)
           continue
        print("       --> %s"%pof)
        data_dir,data_inv = readout_rcons(pof,rctype)
        V0dir_i,V1dir_i,kdir_i,Gdir_i = data_dir
        V0inv_i,V1inv_i,kinv_i,Ginv_i = data_inv
        if kdir_i == {}:
           print("           %s rate constant NOT FOUND!!\n"%rctype)
           continue
        dV0dir[rcname] = V0dir_i
        dV1dir[rcname] = V1dir_i
        dkdir[rcname]  = kdir_i
        dGdir[rcname]  = Gdir_i

        dV0inv[rcname] = V0inv_i
        dV1inv[rcname] = V1inv_i
        dkinv[rcname]  = kinv_i
        dGinv[rcname]  = Ginv_i
        count += 1
    print("")
    if count == 0: return

    #-----------------------#
    # Print table with data #
    #-----------------------#
    if   len(stemps)   == 1:
         props,ml = "reaction,", max([len(n) for n in srcnames]+[8])
    elif len(srcnames) == 1:
         props,ml = "T (K),"   , 8
    else:
         props,ml = "reaction,", 8

    frow   = " %%%is | %%7s | %%7s | %%9s | %%9s | %%14s "%ml
    props1 = props + "V0_dir,V1_dir,G_dir(V),G_dir(p),k_dir"
    props2 = props + "V0_inv,V1_inv,G_inv(V),G_inv(p),k_inv"
    string  = "Columns:\n"
    string += "     - V0       : difference in electronic energy between the transition\n"
    string += "                  state and the reactant(s) [kcal/mol]\n"
    string += "     - V1       : V0 corrected with vibrational zero-point energy [kcal/mol]\n"
    string += "     - G(V)     : Gibbs free energy of activation [kcal/mol] calculated\n"
    string += "                  for a volume per molecule of V = 1cm^3\n"
    string += "     - G(p)     : Gibbs free energy of activation [kcal/mol] calculated\n"
    string += "                  for a volume per molecule of V = kB*T/p0 with p0=1 bar\n"
    string += "     - k        : reaction rate constant\n"
    string += "                  * unimolecular reaction (u) --> [1/s]\n"
    string += "                  * bimolecular  reaction (b) --> [cm^3/molecule/s]\n"
    string += "\n"
    string += "     - _dir     : direct /forward  reaction\n"
    string += "     - _inv     : inverse/backward reaction\n"
    string += "\n"
    for line in string.split("\n"): print("    %s"%line)

    # A) Several reactions, one temperature
    if   len(stemps  ) == 1: case = 1
    elif len(srcnames) == 1: case = 2
    else                   : case = 1
    srcnames = sorted(dV0dir.keys())

    string = ""
    head1  = frow%tuple(props1.split(","))
    head2  = frow%tuple(props2.split(","))
    if case == 2:
       string += "REACTION: %s\n\n"%srcnames[0]
       table1  = "-"*len(head1)+"\n"+head1 +"\n"+"-"*len(head1)+"\n"
       table2  = "-"*len(head2)+"\n"+head2 +"\n"+"-"*len(head2)+"\n"
       boolt2  = False
    for stemp in stemps:
        if case == 1:
           string += "TEMPERATURE: %s K\n\n"%stemp
           table1  = "-"*len(head1)+"\n"+head1 +"\n"+"-"*len(head1)+"\n"
           table2  = "-"*len(head2)+"\n"+head2 +"\n"+"-"*len(head2)+"\n"
           boolt2  = False
        for rcname in srcnames:
            # Get data
            idata = (rcname,stemp,dV0dir,dV1dir,dkdir,dGdir,\
                                  dV0inv,dV1inv,dkinv,dGinv)
            try   : odata = get_rconstr(*idata)
            except: return
            V0dir,V1dir,kdir,G1dir,G2dir = odata[0]
            V0inv,V1inv,kinv,G1inv,G2inv = odata[1]
            if odata[2]: boolt2 = True
            # unimolecular or bimolecular label
            if len(dchem[rcname][0]) == 1: kdir = kdir+" (u)"
            else                         : kdir = kdir+" (b)"
            if " - " not in kinv:
               if len(dchem[rcname][2]) == 1: kinv = kinv+" (u)"
               else                         : kinv = kinv+" (b)"
            # add to table
            if case == 1: idata = [rcname]
            else        : idata = [stemp ]
            idata1 = idata+[V0dir,V1dir,G1dir,G2dir,kdir]
            idata2 = idata+[V0inv,V1inv,G1inv,G2inv,kinv]
            table1 += frow%tuple(idata1)+"\n"
            table2 += frow%tuple(idata2)+"\n"
        if case == 1:
           table1 += "-"*len(head1)+"\n\n"
           table2 += "-"*len(head2)+"\n\n"
           string += table1
           if boolt2: string += table2
           string += "\n"
    if case == 2:
       table1 += "-"*len(head1)+"\n\n"
       table2 += "-"*len(head2)+"\n\n"
       string += table1
       if boolt2: string += table2
       string += "\n"

    for line in string.split("\n"): print("    %s"%line)
#==================================================#



#==================================================#
def main(idata,status,case,targets):

    stat2check = []
    mustexist  = []
    tocreate   = []
    #-------------------------------------------------------#
    # Read Pilgrim input files, check file/folder status    #
    # and expand tuple 'case'                               #
    #-------------------------------------------------------#
    # expand data
    (dctc,dimasses), ltemp, dpath, (dtesLL,dtesHL), dchem, dkmc, ddlevel = idata
    # status ok?
    fstatus = status_check(status,stat2check)
    if fstatus == -1: exit()
    # existency of folders
    fstatus = ffchecking(mustexist,tocreate)
    if fstatus == -1: exit()
    # expand case
    (dof,hlf,plotfile),dlevel,software = case
    #-------------------------------------------------------#

    #---------------------------------#
    # files with data and output file #
    #---------------------------------#
    if len(targets) < 1 :
       print("    --> Type of output REQUIRED! [pfn/path/rcons]")
       return

    option = targets.pop(0).lower().strip()
    print("    --> Option               : %s"%option)
    if   option == "pfn"  : summary_pfn(targets,ltemp,dctc,dlevel)
    elif option == "path" : summary_path(targets,ltemp,dpath,dctc,dlevel)
    elif option == "rcons": summary_rcons(targets,ltemp,dchem,dlevel)
    else                  : print("        Unknown option!")
    print("")



    return
#==================================================#



#---------------------------------------------------------------#


