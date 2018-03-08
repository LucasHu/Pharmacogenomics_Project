from cmapPy.pandasGEXpress import parse
import pandas as pd
import sys
from rdkit import Chem
import numpy as np


import Feature_generation as FeatureGeneration

# Read CMap drug information file, and map...
def readDrugInformation(fname = "GSE92742_Broad_LINCS_pert_info.txt"):

    #broad_ID_list = []
    pert_ID_list = []
    InChiKey_list = []
    PubChemID_list = []
    smile_list = []


    #broad_ID_InChiKey_dict = {}
    pert_ID_InChiKey_dict = {}

    with open(fname) as fp:
        for line in fp:
            # if line.startswith("BRD"):
            #     items = line.split("\t")
            #     broad_ID_comb = (items[0].rstrip()).split("-")
            #     broad_ID = broad_ID_comb[0] + "-" + broad_ID_comb[1]
            #     broad_ID_list.append(broad_ID)
            #     smile_list.append(items[8].rstrip())
            #     InChiKey_list.append(items[9].rstrip())
            #     PubChemID_list.append(items[10].rstrip())
            #     broad_ID_InChiKey_dict[broad_ID] = items[9].rstrip()
            items = line.split("\t")
            pert_ID = items[0]
            pert_ID_list.append(pert_ID)
            inchi_key = items[5]
            smile = items[6]
            InChiKey_list.append(inchi_key)
            pert_ID_InChiKey_dict[pert_ID] = inchi_key




    #return broad_ID_list, InChiKey_list, PubChemID_list, smile_list, broad_ID_InChiKey_dict
    return pert_ID_InChiKey_dict

# a helper function to read CMap signature file and generate the list of signature ID based on the name of small molecules from DrugBank...
def read_CMap_generate_sig_IDs(fname = "GSE92742_Broad_LINCS_sig_info.txt"):

    with open(fname) as fp:
        next(fp)
        for line in fp:
            line = line.rstrip()
            items = line.split("\t")
            sig_ID = items[0]
            drug_name = items[2]
            print drug_name

# A helper function to generate a dictionary that the key is pertID and the value is drugbank ID...
def generate_broadpert_DrugBank_dict():

    _, drug_name_set = FeatureGeneration.read_DrugBank()

    # read pertID_InCHiKey_Dict mapping file...
    pert_ID_InChiKey_dict = readDrugInformation() #pertID and InChiKey overlap is 51384

    # read DrugBankID and InChiKey mapping file...
    _, _, _, INCHI_KEY_ID_drugbank_ID_dict = FeatureGeneration.read_sdf_file() #DrugBankID and InChiKey overlap is 8738

    # map pertID to DrugBankID...
    pert_ID_DrugBank_ID_dict = {}

    # the drugs here are the one in drug-protein interaction list
    for drug in pert_ID_InChiKey_dict:
        if pert_ID_InChiKey_dict[drug] in INCHI_KEY_ID_drugbank_ID_dict:
            if INCHI_KEY_ID_drugbank_ID_dict[pert_ID_InChiKey_dict[drug]] in drug_name_set:
                pert_ID_DrugBank_ID_dict[drug] = INCHI_KEY_ID_drugbank_ID_dict[pert_ID_InChiKey_dict[drug]] # pertID and DrugBankID overlap is 3786
        else:
            pass

    return pert_ID_DrugBank_ID_dict


if __name__ == '__main__':
    def main():
        #read drug protein interaction data form drugbank.
        protein_drug_list, drug_name_list = FeatureGeneration.read_DrugBank()

        pertID_drugbankID_dict = generate_broadpert_DrugBank_dict()

        #put all pertID with drugBank drug-protein interaction data into a dictionary...
        pertID_with_drugbank_interaction_dict = {}

        selected_drug_protein_list_with_CMap_data = {}

        for keys in pertID_drugbankID_dict:
            if pertID_drugbankID_dict[keys] in drug_name_list:
                pertID_with_drugbank_interaction_dict[keys] = pertID_drugbankID_dict[keys]

            #find the drug-protein pairs with drugs found in CMap...
            for drug_protein_pair in protein_drug_list:
                if pertID_drugbankID_dict[keys] == drug_protein_pair[1]:
                    selected_drug_protein_list_with_CMap_data[keys] = drug_protein_pair

        print "selected_drug_protein_list_with_CMap_data_overlap: " + str(len(selected_drug_protein_list_with_CMap_data))
        sys.exit()

        # a way to generate mol objecy from SMILE string directly...
        m2 = Chem.MolFromSmiles('C1CCC1')
        Mogen2_matrix = FeatureGeneration.generate_fingerprint("Morgan2", [m2])

        # play with GEO dataset..
        sig_info = pd.read_csv("GSE92742_Broad_LINCS_sig_info.txt", sep="\t")

        selected_sig_id_list = []
        test = []
        # get the ids for signature IDs for those perturbation drugs in both drug-target interaction pairs and CMap... ~ 2700
        for key in pertID_with_drugbank_interaction_dict:

            selected_sig_id_list.append(sig_info["sig_id"][sig_info["pert_id"] == key])

        gene_info = pd.read_csv("GSE92742_Broad_LINCS_gene_info.txt", sep="\t", dtype=str)

        landmark_gene_row_ids = gene_info["pr_gene_id"][gene_info["pr_is_lm"] == "1"]

        my_col_metadata = parse("GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx", col_meta_only=True)

        print my_col_metadata
        print type(my_col_metadata)
        print np.shape(my_col_metadata)

        #vorinostat_only_gctoo = parse("GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx", cid=vorinostat_ids)



        #test = vorinostat_only_gctoo.data_df.as_matrix()
        #print type(test)
        #print np.shape(test)
        #print test


        sys.exit()



if __name__ == "__main__":
    main()

# read the drug name from DrugBank database...
protein_drug_from_DrugBank_list, drug_name_list = FeatureGeneration.read_DrugBank()






vorinostat_ids = sig_info["sig_id"][sig_info["pert_iname"] == "vorinostat"]
print vorinostat_ids

gene_info = pd.read_csv("GSE70138_Broad_LINCS_gene_info_2017-03-06.txt", sep="\t", dtype=str)

landmark_gene_row_ids = gene_info["pr_gene_id"][gene_info["pr_is_lm"] == "1"]

my_col_metadata = parse("GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx", col_meta_only=True)
vorinostat_only_gctoo = parse("GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx", cid=vorinostat_ids)
landmark_only_gctoo = parse("GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx", rid = landmark_gene_row_ids)

vorinostat_sig_id_info =  sig_info[sig_info["pert_iname"] == "vorinostat"]

vorinostat_sig_id_info.set_index("sig_id", inplace=True)
