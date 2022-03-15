import pandas as pd
import numpy as np
import scipy.stats
import scipy.optimize
import os
def fsigmoid(x, a, b):
    return 1.0 / (1.0 + np.exp(-a*(x-b)))
def flat_line(x, a):
    return(a)

def bluues(directory, pdb_id):
    output_filename = "%s/%s_ElectFeatures.txt"%(directory, pdb_id)
    #if os.path.isfile(output_filename) == True:
    #    return

    with open("%s/%s_bluues.titration"%(directory, pdb_id), "r") as inData:
        inData = inData.readlines()

    ph_values = []
    residues = []
    for line in inData:
        if "ATOM" in line:
            prev_resn = "1"
            prev_atom = "N"
            these_values = []
        else:
            line = line.strip().split()
            #print(line)
            #print(prev_resn, line)
            if (line[2] == prev_resn and line[0] == prev_atom):
                these_values.append(float(line[-1]))
                prev_resn = line[2]
                prev_atom = line[0]
                this_res = line[1] + "_" + line[2]
            elif len(these_values) == 41:
                residues.append(this_res)
                ph_values.append(these_values)
                #print(len(these_values), this_res)
                these_values = []
                these_values.append(float(line[-1]))
                prev_resn = line[2]
                prev_atom = line[0]
                this_res = line[1] + "_" + line[2]
            elif len(these_values) == 82:
                residues.append(this_res)
                ph_values.append(these_values[0:41])
                residues.append(this_res)
                ph_values.append(these_values[41:])
                these_values = []
                these_values.append(float(line[-1]))
                prev_resn = line[2]
                prev_atom = line[0]
                this_res = line[1] + "_" + line[2]
            else:
                these_values = []
                these_values.append(float(line[-1]))
                prev_resn = line[2]
                prev_atom = line[0]
                this_res = line[1] + "_" + line[2]
            #print(this_res, prev_atom, prev_resn)

    residues.append(this_res)
    ph_values.append(these_values)
    #lengths = [len(x) for x in ph_values]
    #print(lengths)
    ph_values = np.asarray(ph_values)

    collated_ph = pd.DataFrame(index = np.arange(-2,18.5, 0.5), columns = sorted(set(residues)))
    for entry in sorted(set(residues)):
        indices = [i for i, x in enumerate(residues) if x == entry]
        #print(entry)
        these_ph_vals = ph_values[indices, :]
        these_ph_vals = np.sum(these_ph_vals, axis = 0)
        #print(these_ph_vals)
        collated_ph[[entry]] = these_ph_vals.reshape(41,1)

    #collated_ph.to_csv("bluues/%s_AllTitr.txt"%pdb_id, sep = "\t", index = True, header = True)

    #calculate the destabilizing rank (higher number is more destablizing) and the stabilizing rank (inverse order of destabilizing)
    with open("%s/%s_bluues.solv_nrg"%(directory, pdb_id), "r") as inData:
        inData = inData.readlines()
    solv_energy = []; prev_resn = 0
    for i in range(0, len(inData)):
        line = inData[i]
        if "ATOM" in line:
            prev_resn = inData[i+1].strip().split()[5]#"1"
            these_values = []
        elif "SOLV NRG" in line:
            line = line.strip().split()
            #print(line)
            #print(prev_resn, line)
            if (line[5] == prev_resn):
                these_values.append(float(line[-1]))

            else:
                solv_energy.append([this_res, np.sum(these_values)])
                #print(len(these_values), this_res)
                these_values = []
                these_values.append(float(line[-1]))
            prev_resn = line[5]
            this_res = line[4] + "_" + line[5]
    solv_energy.append([this_res, np.sum(these_values)])
    new_solv_energy = pd.DataFrame.from_records(solv_energy, columns=["Res", "SolvEnergy"])
    new_solv_energy = new_solv_energy.sort_values(by=["SolvEnergy"], ascending=False).reset_index(drop=True)
    new_solv_energy.reset_index(inplace=True)
    new_solv_energy["StabRank"] = np.digitize( new_solv_energy["index"], bins = np.linspace( 0, np.max(new_solv_energy["index"])+1, 6 ) )/5 #instead of the actual bin (1-5), it converts bin 1 to a value of 1 and the last bin to a value of 0.2; higher StabRank numbers are more stabilizing
    new_solv_energy["DestabRank"] = np.digitize( new_solv_energy["index"], bins = np.linspace( np.max(new_solv_energy["index"])+1,0, 6 ) )/5
    new_solv_energy = new_solv_energy[["Res", "DestabRank", "StabRank", "SolvEnergy"]]
    new_solv_energy.set_index("Res", inplace = True)

    moments = pd.DataFrame(np.zeros((len(sorted(set(residues))), 10)), index = sorted(set(residues)), columns = ["mu2", "mu3", "mu4", "ResSigDev", "pKa_shift", "dpKa_desolv", "dpKa_bg", "dpKa_titr", "GBR6", "SolvExp" ])
    #calcualted from bluues *.pka output file
    #pKa: the calculated pKa of the titratable atom, 
    #pKa_shift: the model compound pKa - calculated pKa
    #dpKa_desolv: the shift in pKa due to desolvation, 
    #dpKa_bg: the shift in pKa due to the interaction with other charges in the molecule with all titratable sites in their neutral state, 
    #dpKa_titr: the shift in pKa due to the interaction between titratable sites, 
    #GBR6: the GB radius according to the GBR6 model (a measure of the depth of the titratable site in the molecular structure)
    #SolvExp: TBD
    collated_ph = collated_ph.T
    for index, entry in collated_ph.iterrows():
        #print(index, entry)
        # abs diff/0.5 are the derivatives at that point; the abs is necessary because this isn't net charge but ionization data
        deriv = np.abs(np.diff(entry))/0.5
        moments.at[index, "mu2"] = scipy.stats.variation(deriv)
        moments.at[index, "mu3"] = scipy.stats.skew(deriv)
        moments.at[index, "mu4"] = scipy.stats.kurtosis(deriv)
        ## below was used in 2021 pub. 
        ## A decent explanation of diff between calc. is here:
        ## https://www.mathworks.com/matlabcentral/answers/332353-difference-between-third-moment-skewness-and-e-x-3
        #moments.at[index, "mu2"] = scipy.stats.moment(deriv, moment = 2)
        #moments.at[index, "mu3"] = scipy.stats.moment(deriv, moment = 3)
        #moments.at[index, "mu4"] = scipy.stats.moment(deriv, moment = 4)

        #calculate deviation from a fitted sigmoidal curve
        popt, pcov = scipy.optimize.curve_fit(fsigmoid, list(collated_ph), entry, method='dogbox', bounds=([-3, -4],[3, 15])) #, bounds=([0., 600.],[0.01, 1200.])
        residuals = np.sum( (entry - fsigmoid(list(collated_ph), *popt))**2 )
        #check the flat line just in case for res like Tyr
        popt2, pcov2 = scipy.optimize.curve_fit(flat_line, list(collated_ph), entry, method='dogbox', bounds=([0],[1]))
        residuals2 = np.sum( (entry - flat_line(list(collated_ph), *popt2))**2 )
        moments.at[index, "ResSigDev"] = np.min([residuals, residuals2])


    with open("%s/%s_bluues.pka"%(directory, pdb_id), "r") as inData:
        inData = inData.readlines()
    for line in inData[1:]:
        #print(line)
        line = line.strip().split()
        moments.at[line[1] + "_" + line[2], "pKa_shift"] += float(line[4]) - float(line[3]) #pKa_0-pKa
        moments.at[line[1] + "_" + line[2], "dpKa_desolv"] += float(line[5]) #dpKa^self
        moments.at[line[1] + "_" + line[2], "dpKa_bg"] += float(line[6])#dpKa^bg
        moments.at[line[1] + "_" + line[2], "dpKa_titr"] += float(line[7])#dpKa^ii
        moments.at[line[1] + "_" + line[2], "GBR6"] += float(line[8]) #uses a correction from GBR6 in .gbr file
        moments.at[line[1] + "_" + line[2], "SolvExp"] += float(line[9]) #Solv. exp.
    #print(moments)
    moments = moments.merge(new_solv_energy, left_index = True, right_index = True)

    moments.to_csv(output_filename, sep = "\t", header = True, index = True)




import subprocess
def Rosetta(out_dir, pdb_id):
    output_filename = "%s/%s_ByResValues.txt"%(out_dir, pdb_id)
    #if os.path.isfile(output_filename) == True:
    #   return
    
    #print("\n\nBSA work:")
    filename1 = "%s/StdOutputScore.txt"%(out_dir)
    #buried surface area/dssp/rotamer probabilites and energy reductions are for each residue - including metals and can be added to the ByResScore file
    bsa = subprocess.check_output(["sed", "-n", "/============Begin report for BSA==================/,/TOTAL BURIED AREA/{/^pose/d;/TOTAL BURIED AREA/d;p;}", filename1])
    bsa = bsa.decode("utf-8").strip().split("\n")
    bsa = [x.split("[0m")[-1].split("\t") for x in bsa[2:-1]]
    bsa = [x[1] for x in bsa]
    #print(bsa)

    #print("\n\nSASA work:")
    filename2 = "%s/%s.pdb"%(out_dir, pdb_id)
    sasa = subprocess.check_output(["grep", "res_sasa_", filename2])
    sasa = sasa.decode("utf-8").strip().split("\n")
    sasa = [x.split()for x in sasa]
    sasa_nums = [int(x[0][9:]) for x in sasa]
    sasa = [x[1] for x in sasa]
    #print(sasa[0:10], sasa_nums[0:10])
    sasa = [x for _, x in sorted(zip(sasa_nums, sasa))]
    #print(sasa[0:10])

    #print("\n\nENERGY work:")
    by_res_terms = subprocess.check_output(["sed", "-n", "/^pose/,/END_POSE_ENERGIES/{/^pose/d;/END_POSE_ENERGIES/d;p;}", filename2])
    by_res_terms = by_res_terms.decode("utf-8").strip().split("\n")
    by_res_terms = ["\t".join(x.replace(":CtermProteinFull", "").replace(":NtermProteinFull", "").split()) for x in by_res_terms]
    #print(by_res_terms[0:4])
    
    #print("\n\nCONCAT work:")
    all_data = zip(by_res_terms, bsa, sasa)
    all_data = ["\t".join(x) for x in all_data]
    #for entry in all_data:
    #    print(entry)
    labels = subprocess.check_output(["grep", "label", filename2])
    labels = labels.decode("utf-8").strip().split("\n")[0].split()
    labels = "\t".join(labels) + "\tBSA\tSASA\n"
    with open(output_filename, "w+") as outData:    
        outData.write(labels)
        outData.write("\n".join(all_data))
