from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors
import numpy as np
#from sklearn.svm import SVC
#from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
#from sklearn.model_selection import train_test_split, StratifiedKFold, GridSearchCV
#from sklearn.preprocessing import StandardScaler
#from sklearn.metrics import accuracy_score, cohen_kappa_score, matthews_corrcoef
#from sklearn.externals import joblib

def read_sdf_file(directory):
    mols = []
    for mol in Chem.SDMolSupplier(directory):
        if mol is not None:
            mols.append(mol)

    return mols

# convert fingerprint BitVec to np.array
def rdkit_numpy_convert(fp):
    output = []
    for f in fp:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(f, arr)
        output.append(arr)
    return np.asarray(output)

# generate fingerprint and return as a np.array
def generate_fingerprint(fingerprint_type, mols):
    # generate binary Morgan fingerprint with radius 2
    if fingerprint_type == "Morgan2":
        fp = [AllChem.GetMorganFingerprintAsBitVect(m, 2) for m in mols]
    else:
        print "the fingerprint is not supported here!"

    print fp
    x = rdkit_numpy_convert(fp)
    return x

#get extra chemical descriptors
def get_descriptors(mols):
    descr = []
    for m in mols:
        descr.append([Descriptors.MolLogP(m),
                    Descriptors.TPSA(m),
                    Descriptors.NHOHCount(m),
                    Descriptors.NOCount(m),
                    Descriptors.NumHAcceptors(m),
                    Descriptors.NumHDonors(m),
                    Descriptors.NumRotatableBonds(m),
                    Descriptors.NumHeteroatoms(m),
                    Descriptors.FractionCSP3(m)])
    descr = np.asarray(descr)
    return descr


def main():
    # the path of sdf file
    fname = "/Users/lucasminghu/Desktop/Pharmacogenomics/structures.sdf"

    molecules = read_sdf_file(fname)
    


if __name__ == "__main__":
    main()
