import os
from django.shortcuts import render
from .forms import QueryForm
from django.http import HttpResponse
from django.db import connection
import pandas as pd
import math 
from Bio.PDB import PDBParser


def select_pdb_file(uniprot_id, pdb_directory):
    # Get a list of all files in the specified directory
    files = os.listdir(pdb_directory)
    matching_files = 0
    # Iterate through each file
    for file_name in files:
        # Check if the file name contains the Uniprot ID
        if uniprot_id in file_name:
            # Check if the file is a PDB file
            if file_name.endswith(".pdb"):
                # Return the full path to the PDB file
                
                A = f'/static/dbQuery/alphafold_pdb/{file_name}'
                B = f'{pdb_directory}{file_name}'
                
                matching_files +=1
    if matching_files > 1:
        A = None
        B = None    
    return A, B
    
    # If no matching PDB file is found, return None

def extract_protein_sequence(pdb_file_path):
    sequence = ''
    amino_acids = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
        'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    with open(pdb_file_path, 'r') as file:
        record_lines = []
        for line in file:
            if line.startswith('SEQRES'):
                record_lines.append(line)

        for line in record_lines:
            parts = line.split()[4:]  # Skip SEQRES record number, chain ID, and residue count
            for part in parts:
                sequence += amino_acids.get(part, '?')  # Use '?' for unknown or non-standard residues

    return sequence


def calculate_distance(coord1, coord2):
    """Calculate the Euclidean distance between two 3D points."""
    return math.sqrt((coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2 + (coord1[2] - coord2[2])**2)

def extract_atom_info_and_min_distance(pdb_file_path, template_residue_number,template_chain, difference_in_position):
    """
    Extracts each atom's residue number and 3D coordinates, then calculates the
    shortest distance between atoms in all residues and the atoms in a specified variant residue.
    """
    atoms_info = {}
    variant_atom_coordinates = []

    with open(pdb_file_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM"):
                atom_type = line[12:16].strip()
                chain_id = line[21].strip()
                residue_number = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coordinate = (x, y, z)

                # Collect all atom coordinates for the variant residue on the specified chain
                if residue_number == template_residue_number and chain_id == template_chain:
                    variant_atom_coordinates.append(coordinate)

                # Use (residue_number, chain_id) as a unique identifier
                residue_key = (residue_number, chain_id)
                if residue_key not in atoms_info:
                    atoms_info[residue_key] = {
                        "residue_number": residue_number + difference_in_position,
                        "chain_id": chain_id,
                        "atom_coordinates": [],
                        "min_distance_to_variant": None  # To be calculated later
                    }
                atoms_info[residue_key]["atom_coordinates"].append(coordinate)

    # Calculate the shortest distance for each residue to the variant residue
    if variant_atom_coordinates:
        for residue_key, atom_info in atoms_info.items():
            min_distance = float('inf')  # Initialize with a large number
            for coord in atom_info["atom_coordinates"]:
                for variant_coord in variant_atom_coordinates:
                    distance = calculate_distance(coord, variant_coord)
                    if distance < min_distance:
                        min_distance = distance
            atom_info["min_distance_to_variant"] = min_distance

    # Convert to a DataFrame
    atoms_info_DF = pd.DataFrame.from_dict(atoms_info, orient="index")
    atoms_info_DF = atoms_info_DF.rename(columns={"residue_number": "Position"})
    atoms_info_DF["min_distance_to_variant"] = pd.to_numeric(
        atoms_info_DF["min_distance_to_variant"], errors="coerce"
    )
    atoms_info_DF = atoms_info_DF.round({"min_distance_to_variant": 2})
    
    # Condense rows: keep only the shortest distance for each residue number
    condensed_atoms_info_DF = atoms_info_DF.loc[
        atoms_info_DF.groupby("Position")["min_distance_to_variant"].idxmin()
    ].reset_index(drop=True)

    return atoms_info, condensed_atoms_info_DF

def categorize_atoms_by_distance_adjusted(atoms_info, difference_in_position):
    """
    Categorizes atom residues into groups based on their distance to the target CA atom.
    Ensures a residue number is included in all matched groups.
    The groups are: <=5, <=10, <=15, <=20, <=25 angstroms.
    """
    # Initialize groups with sets to avoid duplicate residue numbers
    distance_groups = {
        'group_1': set(),
        'group_2': set(),
        'group_3': set(),
        'group_4': set(),
        'group_5': set(),
        'group_6': set(),
    }

    for atom_info in atoms_info.values():
        distance = atom_info.get('min_distance_to_variant', None)
        residue_number = atom_info.get('residue_number') - difference_in_position

        # Use if-else structure to include residue numbers in all matching groups
        if distance is not None:
            if distance >0:
                if distance <= 5:
                    distance_groups['group_1'].add(residue_number)
                if distance <= 10:
                    distance_groups['group_2'].add(residue_number)
                if distance <= 15:
                    distance_groups['group_3'].add(residue_number)
                if distance <= 20:
                    distance_groups['group_4'].add(residue_number)
                if distance <= 25:
                    distance_groups['group_5'].add(residue_number)
                if distance <= 30:
                    distance_groups['group_6'].add(residue_number)

    # Convert sets back to lists for easier processing/use later
    for group in distance_groups:
        distance_groups[group] = list(distance_groups[group])
        
    group_1 = distance_groups['group_1']
    group_2 = distance_groups['group_2']
    group_3 = distance_groups['group_3']
    group_4 = distance_groups['group_4']
    group_5 = distance_groups['group_5']
    group_6 = distance_groups['group_6']

    return group_1,group_2,group_3,group_4,group_5,group_6



def query_topology(params):
    with connection.cursor() as cursor:
        # Query to get topology and pdb file path
        cursor.execute("""
        SELECT
            CASE
                WHEN gene_names.uniprot_id = %(user_protein)s THEN %(user_protein)s
                WHEN gene_names.gene_name = %(user_protein)s THEN gene_names.uniprot_id
            END AS uniprot_id,
            gene_names.gene_name,
            topology.region_type,
            topology.region_number,
            topology.start,
            topology.end,
            structure.structure AS structure_path
        FROM
            gene_names
        LEFT JOIN
            topology ON gene_names.uniprot_id = topology.uniprot_id
        LEFT JOIN
            structure ON gene_names.uniprot_id = structure.uniprot_id
        WHERE
            gene_names.uniprot_id = %(user_protein)s OR gene_names.gene_name = %(user_protein)s;

        """, params)

        # get results
        rows = cursor.fetchall()

        topology_df = pd.DataFrame(rows, columns = [col[0] for col in cursor.description])
        topology_df.rename(columns={"region_type": "Type", "start": "Start", "end": "End"}, inplace=True)

        if topology_df["Type"].isna().any():
            return topology_df
        
        topology_df.sort_values(by="Start", inplace=True)
        topology_df["Start"] = topology_df["Start"].astype(int)
        topology_df["End"] = topology_df["End"].astype(int)

    return topology_df


def query_ptm(params):
    with connection.cursor() as cursor:
        # Query to get all ptms
        cursor.execute("""
        SELECT
            CASE
                WHEN gene_names.uniprot_id = %(user_protein)s THEN %(user_protein)s
                WHEN gene_names.gene_name = %(user_protein)s THEN gene_names.uniprot_id
            END AS uniprot_id,
            gene_names.gene_name,
            alternative_gene_names.alternative_gene_name,
            protein_names.protein_name,
            ptm.ptm_type,
            ptm.start AS ptm_start,
            ptm.end AS ptm_end,
            ptm.residue AS wt_residue,
            CONCAT(ptm_source.source, ": ", ptm_source.id) AS ptm_source,
            ptm_source.evidence_code, 
            CASE
                WHEN ptm.residue = %(user_wt)s THEN 'True'
                ELSE 'False'
            END AS wt_check
        FROM
            gene_names
        LEFT JOIN 
            alternative_gene_names ON gene_names.uniprot_id = alternative_gene_names.uniprot_id
        LEFT JOIN 
            protein_names ON gene_names.uniprot_id = protein_names.uniprot_id
        LEFT JOIN
            ptm ON protein_names.uniprot_id = ptm.uniprot_id
        LEFT JOIN 
            ptm_source ON ptm.uniprot_id = ptm_source.uniprot_id
            AND ptm.ptm_type = ptm_source.ptm_type
            AND ptm.start = ptm_source.start
            AND ptm.end = ptm_source.end
            AND ptm.residue = ptm_source.residue
        WHERE
            gene_names.uniprot_id = %(user_protein)s OR gene_names.gene_name = %(user_protein)s;
        """, params)

        # get results from query and return as a DataFrame
        rows = cursor.fetchall()
        # print(rows)
        # checks if valid Uniprot ID or gene name exits in database
        if len(rows) == 0:
            return(None)

        ptm_df = pd.DataFrame(rows, columns = [col[0] for col in cursor.description])
        # ptm_df.sort_values(by = "ptm_start", inplace=True)

        # Condense the DataFrame such that:
        #   alternative gene names are in a list 
        #   ptm_info is a list of dictionaries with {ptm_type, ptm_start, ptm_end, wt_residue, ptm_source, wt_check, evidence_code}
        grouped = ptm_df.groupby(['uniprot_id', 'gene_name', 'protein_name'])

        # Initialize an empty dictionary to store the condensed data
        condensed_dict = {}

        # Iterate through each group
        for group_name, group_data in grouped:
            # extract unique identifiers
            uniprot_id, gene_name, protein_name = group_name
            
            ptm_info_set = set() # create a set to store unique ptm_info dictionaries
            alternative_gene_names = set()  # to store unique alternative gene names
            
            for _, row in group_data.iterrows():
                # Create a hashable version of ptm_info dictionary
                ptm_info_hashable = tuple(row[['ptm_type', 'ptm_start', 'ptm_end', 'wt_residue', 'ptm_source', 'wt_check', 'evidence_code']])
                
                # Add ptm_info to set if it's not already present
                if ptm_info_hashable not in ptm_info_set:
                    ptm_info_set.add(ptm_info_hashable)
                
                # Collect unique alternative gene names
                alternative_gene_names.add(row['alternative_gene_name'])
            
            # Convert set of ptm_info back to list of dictionaries
            ptm_info_list = [dict(zip(['ptm_type', 'ptm_start', 'ptm_end', 'wt_residue', 'ptm_source', 'wt_check', 'evidence_code'], item)) for item in ptm_info_set]
            
            # Create the final dictionary entry
            entry = {
                'uniprot_id': uniprot_id,
                'gene_name': gene_name,
                'alternative_gene_names': list(alternative_gene_names),
                'protein_name': protein_name,
                'ptm_info': ptm_info_list
            }
            
            # Use uniprot_id as the dictionary key
            condensed_dict[uniprot_id] = entry

    return entry

def get_variantDBLink(row):
    dbLink_list = []
    for db in row["db"]:
        link = db
        if db == "ClinVar":
            link = f'https://www.ncbi.nlm.nih.gov/clinvar/{row["snp_id"]}/'
        elif db == "UniProt":
            link = f'https://web.expasy.org/variant_pages/{row["snp_id"]}.html'
        elif db == "dbSNP":
            link = f'https://www.ncbi.nlm.nih.gov/snp/{row["snp_id"]}'
        elif db == "gnomAD":
            link = f'https://gnomad.broadinstitute.org/variant/{row["snp_id"]}'
        elif db == "ExAC":
            link = f'https://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?type=rs&rs={row["snp_id"]}'
        elif db == "TOPMed":
            link = f'https://www.ncbi.nlm.nih.gov/snp/{row["snp_id"]}#frequency_tab'
        elif db == "ClinGen":
            link = f'https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/alleles?name={row["snp_id"]}&skip=0&limit=100'
        elif db == "NCI-TCGA":
            link = f'https://www.ncbi.nlm.nih.gov/snp/{row["snp_id"]}'
            #dbLink_list.append(f'https://portal.gdc.cancer.gov/exploration?filters=%7B%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22variants.dbsnp_rs%22%2C%22value%22%3A%5B%22{row["snp_id"]}%22%5D%7D%7D%5D%2C%22op%22%3A%22and%22%7D')
        elif db == "1000Genomes":
            link = f'https://www.ensembl.org/homo_sapiens/Variation/Explore?v={row["snp_id"]}'
        elif db == "ESP":
            link = f'https://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?type=rs&rs={row["snp_id"]}'
        elif db == "NCI-TCGA Cosmic":
            link = f'https://cancer.sanger.ac.uk/cosmic'
        elif db == "cosmic curated":
            link = f'https://cancer.sanger.ac.uk/cosmic'
        if link != db:
            dbLink_list.append(
                    f'<a href="{link}" target="_blank" title = "Click here to go to the database for {row["snp_id"]}." data-toggle="tooltip">'
                    f'{db} '
                    f'<i class="bi-box-arrow-up-right"></i>'
                    f'</a>'
                )
        else: 
            dbLink_list.append(db)
    return(dbLink_list)


def get_diseaseOmimIDLink(row):
    """
    Given a DataFrame row with any of:
      - row['disease']                (old single‑column style)
      - row['disease_rcv'], row['disease_nonrcv']
    plus row['omim_id'] (a list of IDs),
    returns a flat list of strings of the form
      "DiseaseName (OMIM link…)"
    or just "DiseaseName" if there's no valid OMIM ID.
    """
    omimIDLink_list = []

    # 1) Figure out which 'disease' columns we actually have
    disease_cols = [c for c in row.index if c == 'disease' or c.startswith('disease_')]

    # 2) Gather them all into one list
    diseases = []
    for col in disease_cols:
        vals = row[col]
        if vals is None:
            continue
        # if it's a scalar, make it a list
        if not isinstance(vals, (list, tuple, pd.Series)):
            vals = [vals]
        diseases.extend(vals)
    # drop None, "-" and duplicates while preserving order
    seen = set()
    clean_diseases = []
    for d in diseases:
        if d and d != '-' and d not in seen:
            seen.add(d)
            clean_diseases.append(d)

    # 3) Same for OMIM IDs
    omims = row.get('omim_id', [])
    if omims is None:
        omims = []
    if not isinstance(omims, (list, tuple, pd.Series)):
        omims = [omims]

    # 4) Build links
    for disease in clean_diseases:
        # if there are any valid OMIM IDs, append one entry per ID
        any_linked = False
        for om in omims:
            if om and om != '-':
                any_linked = True
                url = f'https://www.omim.org/entry/{om}'
                omimIDLink_list.append(
                    f'{disease} ('
                    f'<a href="{url}" target="_blank" title="Click to view OMIM entry">'
                    f'OMIM ID: {om}'
                    f' <i class="bi-box-arrow-up-right"></i>'
                    f'</a>)'
                )
        # if we never added an OMIM link for this disease, just show the name
        if not any_linked:
            omimIDLink_list.append(disease)

    return omimIDLink_list


def query_variantAtPos(params):
    with connection.cursor() as cursor:
        cursor.execute("""
        SELECT
            CASE
                WHEN gn.uniprot_id = %(user_protein)s THEN %(user_protein)s
                WHEN gn.gene_name  = %(user_protein)s THEN gn.uniprot_id
            END AS uniprot_id,
            v.snp_id,
            v.position,
            v.wt,
            v.mutant,
            v.ensembl_id,
            vd.database_name AS db,
            vsr.GermlineClassification      AS clinical_significance_rcv,
            vsr.Diseases                    AS disease_rcv,
            vsr.phenotype                    AS phenotype_rcv,
            vsr.GermlineReviewStatus        As review_status_rcv,
            versus.clinical_significance    AS clinical_significance_nonrcv,
            versus.Diseases                AS disease_nonrcv,
            versus.phenotype                AS phenotype_nonrcv,
            versus.GermlineReviewStatus     As review_status_nonrcv,
            p.omim_id,
            ps.pubmed_id                    AS pheno_source,
            t.region_type,
            t.region_number,
            t.start,
            t.end,
            CASE WHEN v.wt = %(user_wt)s THEN 'True' ELSE 'False' END AS wt_check
        FROM gene_names AS gn
        LEFT JOIN variants AS v
            ON gn.uniprot_id = v.uniprot_id
           AND %(user_position)s = v.position
        LEFT JOIN variant_databases AS vd
            ON v.snp_id = vd.snp_id
        LEFT JOIN variant_significance_rcv AS vsr
            ON v.snp_id = vsr.accession_id
        LEFT JOIN variant_significance_nonrcv AS versus
            ON v.snp_id = versus.accession_id
        LEFT JOIN phenotype AS p
            ON v.snp_id = p.snp_id
        LEFT JOIN phenotype_sources AS ps
            ON v.snp_id = ps.snp_id
        LEFT JOIN topology AS t
            ON v.uniprot_id = t.uniprot_id
           AND v.position BETWEEN t.start AND t.end
        WHERE
            gn.uniprot_id = %(user_protein)s
         OR gn.gene_name = %(user_protein)s;
        """, params)

        results = cursor.fetchall()
        cols    = [col[0] for col in cursor.description]
        df      = pd.DataFrame(results, columns=cols)

    # no variants found?
    if df['snp_id'].isna().all():
        return None

    # Convert wt_check to boolean
    df['wt_check'] = df['wt_check'] == 'True'

    # Aggregate one-to-many fields
    grouped = (
        df
        .groupby(['uniprot_id','snp_id','wt','mutant','position','wt_check'], as_index=False)
        .agg({
            'ensembl_id':                    lambda x: list(set(filter(pd.notna, x))),
            'db':                             lambda x: list(set(filter(pd.notna, x))),
            'clinical_significance_rcv':      lambda x: list(set(filter(pd.notna, x))),
            'disease_rcv':                    lambda x: list(set(filter(pd.notna, x))),
            'phenotype_rcv':                    lambda x: list(set(filter(pd.notna, x))),
            'review_status_rcv':             lambda x: list(set(filter(pd.notna, x))),
            'clinical_significance_nonrcv':   lambda x: list(set(filter(pd.notna, x))),
            'disease_nonrcv':                 lambda x: list(set(filter(pd.notna, x))),
            'phenotype_nonrcv':                 lambda x: list(set(filter(pd.notna, x))),
            'review_status_nonrcv':          lambda x: list(set(filter(pd.notna, x))),
            'omim_id':                        lambda x: list(set(filter(pd.notna, x))),
            'pheno_source':                   lambda x: list(set(filter(pd.notna, x))),
            'region_type':                    lambda x: list(set(filter(pd.notna, x))),
            'region_number':                  lambda x: list(set(filter(pd.notna, x))),
            'start':                          lambda x: list(set(filter(pd.notna, x))),
            'end':                            lambda x: list(set(filter(pd.notna, x)))
        })
    )

    condensed_list = []

    for _, row in grouped.iterrows():
        # Build all the variantInfo entries for this row
        vi_entries = []

        if row['clinical_significance_rcv']:
            vi_entries.append({
                'snp_id':               row['snp_id'],
                'db':                   row['db'],
                'db_source':            get_variantDBLink(row),
                'clinical_significance': row['clinical_significance_rcv'],
                'conditions':         row['phenotype_rcv'],
                'disease':              row['disease_rcv'],
                'disease_OMIMID':       get_diseaseOmimIDLink(row),
                'omim_id':              row['omim_id'],
                'pheno_source':         row['pheno_source'],
                'source':               'RCV',
                'review_status':        row['review_status_rcv']
            })

        else:
            vi_entries.append({
                'snp_id':               row['snp_id'],
                'db':                   row['db'],
                'db_source':            get_variantDBLink(row),
                'clinical_significance': row['clinical_significance_nonrcv'],
                'conditions':         row['phenotype_nonrcv'],
                'disease':              row['disease_nonrcv'],
                'disease_OMIMID':       get_diseaseOmimIDLink(row),
                'omim_id':              row['omim_id'],
                'pheno_source':         row['pheno_source'],
                'source':               'non‑RCV',
                'review_status':        row['review_status_nonrcv']
            })

        if not vi_entries:
            continue

        # Create the top‑level dict for this position/mutation
        variant_dict = {
            'position':    row['position'],
            'wt':          row['wt'],
            'mutant':      row['mutant'],
            'wt_check':    row['wt_check'],
            'variantInfo': vi_entries
        }

        # If we've already seen this position/wt/mutant, merge the lists
        existing = next(
            (item for item in condensed_list
             if item['position']==row['position']
             and item['wt']==row['wt']
             and item['mutant']==row['mutant']),
            None
        )

        if existing:
            # extend—not append—a list of entries
            existing['variantInfo'].extend(vi_entries)
        else:
            condensed_list.append(variant_dict)

        # Clean up the condensed list
        for var in condensed_list:
            for info in var['variantInfo']:
                # normalize any [''] → []
                if info.get('conditions') == ['']:
                    info['conditions'] = []
                if info.get('disease') == ['']:
                    info['disease'] = []
                if info.get('clinical_significance') == ['']:
                    info['clinical_significance'] = []
                if info.get('review_status') == ['']:
                    info['review_status'] = []
    return condensed_list


def query_additionalPTMTextInfo(params):
    with connection.cursor() as cursor:
        # Gets all additional ptms
        cursor.execute("""
        SELECT
            CASE
                WHEN gene_names.uniprot_id = %(user_protein)s THEN %(user_protein)s
                WHEN gene_names.gene_name = %(user_protein)s THEN gene_names.uniprot_id
            END AS uniprot_id,
            gene_names.gene_name,
            protein_names.protein_name,
            additional_ptm.description AS addtl_ptm,
            CONCAT (additional_ptm_sources.source, ': ', additional_ptm_sources.id) AS addtl_ptm_source
        FROM
            gene_names
        LEFT JOIN 
            protein_names ON gene_names.uniprot_id = protein_names.uniprot_id
        LEFT JOIN
            additional_ptm on protein_names.uniprot_id = additional_ptm.uniprot_id
        LEFT JOIN
            additional_ptm_sources on additional_ptm.description = additional_ptm_sources.description
            AND additional_ptm.uniprot_id = additional_ptm_sources.uniprot_id
        WHERE
            gene_names.uniprot_id = %(user_protein)s OR gene_names.gene_name = %(user_protein)s;
        """, params)

        # Get results
        results = cursor.fetchall()
        # Get column names
        col_names = [column[0] for column in cursor.description]

        # Post processing data
        resultsdf = pd.DataFrame(results, columns=col_names)
        resultsdf["addtl_ptm_source_htmlLink"] = None
        if resultsdf["addtl_ptm_source"].notnull().any():
            for index, row in resultsdf.iterrows():
                split = row["addtl_ptm_source"].split(": ")
                if split[0] == "PubMed":
                    link = f'https://pubmed.ncbi.nlm.nih.gov/{split[1]}/'
                    row["addtl_ptm_source_htmlLink"] = (
                        f'<a href="{link}" target="_blank" title = "Click here to go to the source." data-toggle="tooltip">'
                        f'{row["addtl_ptm_source"]} '
                        f'<i class="bi-box-arrow-up-right"></i>'
                        f'</a>'
                        
                        )
                elif split[0] == "UniProtKB":
                    link = f'https://www.uniprot.org/uniprotkb/{split[1]}/'
                    row["addtl_ptm_source_htmlLink"] = (
                    f'By sequence similarity - ' 
                    f'<a href="{link}" target="_blank" title = "Click here to go to the source." data-toggle="tooltip">'
                    f'{row["addtl_ptm_source"]} '
                    f'<i class="bi-box-arrow-up-right"></i>'
                    f'</a>'
                    )



        # Group the DataFrame by the 'addtl_ptm' column and apply the custom function
        addtl_ptms = resultsdf.groupby(['addtl_ptm'])['addtl_ptm_source_htmlLink'].agg(list).reset_index()


        return(addtl_ptms)

def query_ptmVarCosite(params):

    with connection.cursor() as cursor:
        cursor.execute("""
        SELECT
            CASE
                WHEN gene_names.uniprot_id = %(user_protein)s THEN %(user_protein)s
                WHEN gene_names.gene_name = %(user_protein)s THEN gene_names.uniprot_id
            END AS uniprot_id,
            gene_names.gene_name
        FROM
            gene_names,
            ptm_variant_sites
        WHERE gene_names.uniprot_id = %(user_protein)s OR gene_names.gene_name = %(user_protein)s;


        """, params)

    # Get results
    results = cursor.fetchall()
    # Get column names
    col_names = [column[0] for column in cursor.description]

    # Post processing data
    resultsdf = pd.DataFrame(results, columns=col_names)
    # print("#### PTM Variant Cosite ####")
    # print(resultsdf)
    return(resultsdf)
    

def condense_ptmInfo(ptm_all):
    """
    Condenses the data from ptm_all such that ptm_source, Link, and htmlLink are in a list from each unique
    ptm_type and Position

    :param ptm_all, pandas DataFrame, post-processed dataframe outputted from query_ptm()
    """
    # Defining a function to remove duplicates and aggregate into a list
    #print(ptm_all)
    aggregate_list = lambda x: list(set(x))

    # Grouping by 'Type' and 'Position(s)' and aggregating into lists
    condensed_df = ptm_all.groupby(['ptm_type', 'wt_residue', 'wt_check', 'ptm_start', 'Position']).agg({
        'ptm_source': aggregate_list,
        'Link': aggregate_list,
        'htmlLink': aggregate_list,
        'evidence_code': aggregate_list, #TODO
    }).reset_index()
    #print(condensed_df)
    condensed_df.sort_values(by="ptm_start", inplace = True)

    return condensed_df

def calculate_position(row):
    """
    Handles positions of PTMs that are not just 1 residue and 
    converts those into a string of ptm_start - ptm end
    :param row, Pandas Series, row from ptm_all dataframe 
    """

    if row['ptm_start'] == row['ptm_end']:
        return int(row['ptm_start'])
    else:
        return f"{row['ptm_start']}-{row['ptm_end']}"
    
def get_sourceLink(row):
    """
    Creates a link for given ptm_source and ptm_id 
    :param row, Pandas Series, row from ptm_all dataframe 
    """
    if row["ptm_source"] is not None:
        parts = row["ptm_source"].split(": ")
        if parts[0] == "PubMed":
            return f'https://pubmed.ncbi.nlm.nih.gov/{parts[1]}/'
        elif parts[0] == "PDB":
            return f'https://www.rcsb.org/structure/{parts[1]}'
        elif parts[0] == "PRIDE":
            return f'https://www.ebi.ac.uk/pride/archive?keyword={parts[1]}'
        else:
            return row["ptm_source"]
    else:
        return None
    
def remove_duplicate_letters(s):
    seen = set()
    output = []
    for char in s:
        if char not in seen:
            seen.add(char)
            output.append(char)
    return ''.join(output)


# Function to check differences between uniprot and alphafold
def check_differences(table, residue_by_position):
    results = []
    for _, row in table.iterrows():
        position = row['Position']
        wt_residue = row['wt_residue']
        residue = residue_by_position.get(position, None)
        if residue is None:
            results.append({
                'Position': position,
                'wt_residue': wt_residue,
                'residue': 'None',
                'status': 'Position not found in AlphaFold Model.'
            })
        elif wt_residue != residue:
            results.append({
                'Position': position,
                'wt_residue': wt_residue,
                'residue': residue,
                'status': 'Difference found in AlphaFold Model.'
            })
    return pd.DataFrame(results)

def generate_jsmol_plddt_script(pdb_file_path):


    # Define AlphaFold pLDDT coloring scheme
    color_ranges = [
        (0, 50, "red"),
        (50, 70, "yellow"),
        (70, 90, "cyan"),
        (90, 100, "blue"),
    ]

    # Parse the PDB file
    residues = {}
    with open(pdb_file_path, "r") as file:
        for line in file:
            if line.startswith("ATOM"):
                residue_number = int(line[22:26].strip())  # Extract residue number
                bfactor = float(line[60:66].strip())  # Extract B-factor (pLDDT)
                residues[residue_number] = bfactor
        

    # Generate JSmol commands for each color range
    script_lines = []
    for min_val, max_val, color in color_ranges:
        residues_in_range = [
            str(res_num)
            for res_num, bfactor in residues.items()
            if min_val <= bfactor < max_val
        ]
        if residues_in_range:
            residue_list = ",".join(residues_in_range)
            script_lines.append(f"select {residue_list}; color {color};")

    # Combine the commands into a single script
    return " ".join(script_lines)




# Function to filter the variant information
# Function to filter and modify the variant information
def filter_and_add_dbsnp(variants, allowed_dbs):
    if not variants:
        return []  # or return None, depending on what your view expects

    for variant in variants:
        for info in variant["variantInfo"]:
            if info["snp_id"].startswith("rs"):
                # Add dbSNP if missing
                if "dbSNP" not in info["db"]:
                    info["db"].append("dbSNP")
                    dbsnp_url = f'<a href="https://www.ncbi.nlm.nih.gov/snp/{info["snp_id"]}" target="_blank" title = "Click here to go to the database for {info["snp_id"]}." data-toggle="tooltip">dbSNP <i class="bi-box-arrow-up-right"></i></a>'
                    info["db_source"].append(dbsnp_url)

                # Filter db and db_source based on allowed_dbs
                filtered_dbs = []
                filtered_sources = []
                for db, source in zip(info["db"], info["db_source"]):
                    if db in allowed_dbs:
                        filtered_dbs.append(db)
                        filtered_sources.append(source)
                info["db"] = filtered_dbs
                info["db_source"] = filtered_sources
    return variants





    
