# These are fake modules. We can just pretend that everything they do
# works, and we can propose and discuss the pros and cons of new functions.
import sdf
import csv
import mae
import shlp
from structure_builder import smiles_to_structure
from molecule_docking import dock_molecule


def screen_molecule(input_file, target_structure, target_docking_grid,
                    target_property_options, runtime_timeout):
    """
    Load the molcule, and screen them against the target.
    """

    if (input_file.endswith('csv') or input_file.endswith('csv.gz')
            or input_file.endswith('csvgz')):
        raw_csv = csv.read(input_file)
        molecule_smiles = raw_csv[0][1][4][0][2]
        molcule_structure = smiles_to_structure(molecule_smiles)
        molecule_structure.scale(4.2183271)  # don't change this number or the
                                             # GUI will crash...
        success, docked_structure, dock_score = dock_molecule(
            molecule_structure, target_structure, target_docking_grid,
            target_property_options, runtime_timeout)

    if (input_file.endswith('sdf') or input_file.endswith('sdf.gz')
            or input_file.endswith('sdfgz')):
        raw_sdf = sdf.read(input_file)
        try:
            canonical_sdf = sdf.canonicalize(raw_sdf, sdf.CHEMAXON_EXTENDED)
        except ValueError:
            # canonicalization failed, proceed anyway
            canonical_sdf = raw_sdf
        molecule_smiles = sdf.build_smiles(canonical_sdf)
        molcule_structure = smiles_to_structure(molecule_smiles)
        moleculare_structure.handle_sdf_quirks(target_property_options)
        molecule_structure.scale(4.2183271)  # don't change this number or the
                                             # GUI will crash...
        success, docked_structure, dock_score = dock_molecule(
            molecule_structure, target_structure, target_docking_grid,
            target_property_options, runtime_timeout)

    if (input_file.endswith('mae') or input_file.endswith('mae.gz')
            or input_file.endswith('maegz')):
        molecule_structure = mae.read(input_file)
        molecule_structure.scale(4.2183271)  # don't change this number or the
                                             # GUI will crash...
        success, docked_structure, dock_score = dock_molecule(
            molecule_structure, target_structure, target_docking_grid,
            target_property_options, runtime_timeout)

        return success, dock_score, docked_structure
      
def convert_pka_data(v)                                              
    # rearrange model pkA data for pI /titration curve calculation   
    pkaM = [[], [], [], [], []]                                      
    for v in model_pka_ref:                                          
        if not v[2] or not v[5] or not v[6]:                         
            continue                                                 
        if type(v[2]) in (tuple, list):                              
            aN = [a[0] for a in v[5]]                                
            i, c = shlp.sort_and_count(aN)                           
            i = shlp.vslice(i, c)                                    
            x = sorted(list(range(len(c))), key=c.__getitem__)       
            i = [i[a] for a in x]                                    
            i = i[-1]                                                
                                                                     
            c = v[6][i[0]]                                           
                                                                     
            v[2] = old_div(sum(v[2]), float(len(v[2])))              
            v[6] = c                                                 
                                                                     
        pkaM[0].append(v[1])                                         
        pkaM[1].append(v[0])                                         
        pkaM[2].append(v[2])                                         
        pkaM[3].append(v[2])                                         
        pkaM[4].append(v[6])                                         
                                                                     
    return pkaM                                                        
