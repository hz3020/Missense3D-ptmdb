from django.shortcuts import render
from .forms import QueryForm
from django.http import HttpResponse, FileResponse
from django.db import connection
import pandas as pd
import os
import json
from Bio.PDB import PDBParser
import tempfile
from django.db import connection

from .utils import select_pdb_file, extract_protein_sequence, query_topology, query_ptm, get_variantDBLink, get_diseaseOmimIDLink, query_variantAtPos, condense_ptmInfo, calculate_position, get_sourceLink, query_additionalPTMTextInfo,query_ptmVarCosite
def process_form(request):
    '''
    Gets user input form and processes it to be displayed in results.html
    '''
    if request.method == "POST":
        form = QueryForm(request.POST, request.FILES)
        if form.is_valid():
            # print("Form is valid")
            uniprot_id = form.cleaned_data["UNIPROT"]
            gene_name = form.cleaned_data["GENENAME"]
            position = form.cleaned_data['POSITION']
            residue = form.cleaned_data["RESIDUE"]


            ### this is for user to upload their own pdb file in the future, right now no upload option yet
            upload_choice = form.cleaned_data['upload_choice']
            if upload_choice == 'upload' and 'upload_file' in request.FILES:
                start_for_pdb_match = form.cleaned_data["start_position"]
                uploaded_file = request.FILES['upload_file']
                Target_Chain = form.cleaned_data["Target_Chain"]
                
                temp_path ="/project/home23/hz3020/variantPSite/dbQuery/static/dbQuery/temp/"
                for filenameeee in os.listdir(temp_path):
                    file_path = os.path.join(temp_path, filenameeee)
                    os.unlink(file_path)
                    
                public_path = "/project/home23/hz3020/variantPSite/dbQuery/static/dbQuery/temp/" + uploaded_file.name
                with open(public_path, 'wb') as public_file:
                    for chunk in uploaded_file.chunks():
                        public_file.write(chunk)
                pdb_file_path2 = "/static/dbQuery/temp/" + uploaded_file.name
            else:
                Target_Chain = 'A'
            
            if uniprot_id :
                input_data = {'user_protein':uniprot_id, 'user_position': position, 'user_wt': residue}
            elif gene_name:
                input_data = {'user_protein':gene_name, 'user_position': position, 'user_wt': residue}




            # Define results dict 
            results = {}
            
            # PTM Query -------
            ptmQueryRes = query_ptm(input_data)

            # checks if valid Uniprot ID or gene name exits in database
            if ptmQueryRes == None:
                return render(request, "dbQuery/error.html", {"uniprot_id": uniprot_id, "gene_name": gene_name, "position": position, "residue": residue})

            # process PTM query to pull out uniprot ID, gene name, and alternative gene names 
            uniprot_id = ptmQueryRes["uniprot_id"]
            gene_name = ptmQueryRes["gene_name"]
            alt_gene_names = ptmQueryRes["alternative_gene_names"]
            protein_name = ptmQueryRes["protein_name"]

            # post-process PTM query to pull out all PTM information
            ptm_all = pd.DataFrame(ptmQueryRes["ptm_info"])
            is_ptmAll_empty = ptm_all['ptm_type'].isna().all()
            results["ptm_all"] = ptm_all
            results["is_ptmAll_empty"] = is_ptmAll_empty
            if is_ptmAll_empty == False:
                ptm_all['Position'] = ptm_all.apply(calculate_position, axis=1) # get the range of position(s) (start-end)           
                ptm_all["Link"] = ptm_all.apply(get_sourceLink, axis = 1) # get the source link 
                ptm_all['htmlLink'] = ptm_all.apply(lambda row: f'<a href="{row["Link"]}" target="_blank" data-toggle="tooltip">{row["ptm_source"]}<i class="bi-box-arrow-up-right"></i></a>' if row['ptm_source'] is not None else "-", axis=1) # render links in html
                condensed_ptm = condense_ptmInfo(ptm_all)
                results["ptm_all"] = condensed_ptm
                evidence_codes = {
                    "ECO:0000269": "Experimental evidence manually curated with published experimental support.",
                    "ECO:0000303": "Non-traceable author statement evidence based on scientific opinion.",
                    "ECO:0000305": "Curator inference based on scientific knowledge or article content.",
                    "ECO:0000250": "Sequence similarity evidence propagated from a related protein.",
                    "ECO:0000255": "Sequence model evidence used in manual assertions.",
                    "ECO:0000256": "Sequence model evidence used in automatic assertions.",
                    "ECO:0000259": "Descendant of ECO:0000256 used in automatic assertions.",
                    "ECO:0000312": "Imported information evidence used in manual assertions.",
                    "ECO:0000313": "Imported information evidence used in automatic assertions.",
                    "ECO:0007744": "Combinatorial evidence inferred from experimental and computational evidence (manual).",
                    "ECO:0007829": "Combinatorial evidence inferred from experimental and computational evidence (automatic).",
                    "ECO:0008006": "Deep learning neural network evidence used in automatic assertions."
                }
                # Function to map evidence codes to descriptions
                def add_descriptions(codes):
                    return [f"{code} ({evidence_codes.get(code, 'Description not found')})" for code in codes]

                # Apply the function to the column
                condensed_ptm['evidence_code'] = condensed_ptm['evidence_code'].apply(add_descriptions)
                condensed_ptm_target = condensed_ptm[condensed_ptm["Position"] == int(input_data["user_position"])]

                # Initialize list to store motif positions around the target
                condensed_ptm_motif = []
                # Iterate through the rows of the DataFrame
                for _, row in condensed_ptm.iterrows():
                    if row["Position"] != int(input_data["user_position"]):
                        if int(row["Position"]) -5 <= int(input_data["user_position"]) <= int(row["Position"]) +4:
                            condensed_ptm_motif.append(row)
                            print(input_data["user_position"])
                            print(row["Position"])
                print(condensed_ptm_target)

                is_ptmTarget_df_empty = condensed_ptm_target.empty
                results["ptm_target"] = condensed_ptm_target
                results["ptm_target_motif"] = pd.DataFrame(condensed_ptm_motif)
                results["is_ptmTarget_df_empty"] = is_ptmTarget_df_empty
  
            # Additional PTM Text Info Query  --------
            additional_ptm = query_additionalPTMTextInfo(input_data)
            is_additionalPTM_df_empty = additional_ptm.empty

                







            # Data for JSmol section -------------------------------------------------------------------------------------------
            with connection.cursor() as cursor:
                cursor.execute("SELECT * FROM variants WHERE uniprot_id = %s", [uniprot_id])
                variant_all = cursor.fetchall()
                # Get column names
                col_names = [column[0] for column in cursor.description]
                # Post processing data
                variant_allDF = pd.DataFrame(variant_all, columns=col_names)
            
            # Topology data Query ------
            topologyQueryRes = query_topology(input_data)
            html_topology = topologyQueryRes[["Type", "Start", "End"]].to_html(classes="table table-bordered table-striped custom-table", index = False) # to be displayed on website 
            html_topology_with_id = html_topology.replace('<table', '<table id="topology-table"')

            # add topology information to results if it exists
            pdb_directory = "/project/data/alphafold/pdb/"
            pdb_file, sequence_file = select_pdb_file(uniprot_id, pdb_directory)

            if sequence_file != None:
                sequence = extract_protein_sequence(sequence_file)
                if (int(input_data['user_position']) <= len(sequence)) & (int(input_data['user_position']) >0):
                    wt_residue_fromSequenceFile = sequence[int(input_data['user_position']) - 1]
                    results["user_position_withinSequence"] = True
                    results["wt_residue_fromSequenceFile"] = wt_residue_fromSequenceFile
                else:
                    results["user_position_withinSequence"] = False
            else:
                sequence = None
            seq_ptm=[]
            if not is_ptmAll_empty:
                for i in ptm_all['ptm_start']:
                    seq_ptm.append(int(i))
                seq_ptm = list(dict.fromkeys(seq_ptm))
                results["seq_ptm"] = seq_ptm

            seq_var=[]
            for i in variant_all:
                seq_var.append(int(i[2]))
            seq_var = list(set(seq_var))
            results['seq_var'] = seq_var
    
            # Remove duplicate variants 
            from .Class import Segment
            from .Class import SegmentList 
            highlight =''
            segment_list = SegmentList()
            
            # for topology coloring button
            if topologyQueryRes["Type"].any():
                # add topology table 
                results["topology"] = html_topology_with_id

                for index, seg in topologyQueryRes.iterrows():
                    if upload_choice == 'upload' and 'upload_file' in request.FILES:
                        start = str(seg["Start"])
                        end = str(seg["End"])    
                        start2 = str(seg["Start"]-(int(start_for_pdb_match)-1))
                        end2 = str(seg["End"]-(int(start_for_pdb_match)-1))
                    else:
                        start = str(seg["Start"])
                        end = str(seg["End"])    
                        start2 = str(seg["Start"])
                        end2 = str(seg["End"])  
                    if seg["Type"] == 'Cytoplasmic':
                        highlight += "select within(group," +str(int(start2))+"-" +str(int(end2)) + '); color lightblue;'
                        segment_list.add_segment(int(start),int(end), "lightblue")
                    elif seg["Type"] == 'Transmembrane':
                        highlight += "select within(group," +str(int(start2))+"-" +str(int(end2)) + '); color aquamarine;'
                        segment_list.add_segment(int(start),int(end), "aquamarine")
                    elif seg["Type"] == 'Extracellular':
                        highlight += "select within(group," +str(int(start2))+"-" +str(int(end2)) + '); color blue;'
                        segment_list.add_segment(int(start),int(end), "blue")
            highlight += ' select ' +str(input_data['user_position'])+ '; color aqua;'   
            uniprot_residue_number = int(input_data['user_position'])

            # for sequence viewer
            from .Class import Spagetti
            from .Class import SpagettiList
            spagetti_list = SpagettiList()
            serialized_segments = [segment.to_dict() for segment in segment_list.get_segments()]
            
            if upload_choice == 'upload' and 'upload_file' in request.FILES:
                difference_in_position = start_for_pdb_match-1
            else:
                difference_in_position = 0
            template_residue_number = int(input_data["user_position"]) - difference_in_position
            results["residue_number"] = template_residue_number
            results["uniprot_residue_number"] = uniprot_residue_number

            # for coloring sequence viewer based on topology    
            if sequence != None:
                for index, residue in enumerate(sequence):
                    check = 0
                    if index == uniprot_residue_number-1:
                        spagetti_list.add_spagetti(residue, int(uniprot_residue_number),int(uniprot_residue_number)-int(difference_in_position), 'orange' )
                        check = 1
                    
                    else: 
                        for segment in serialized_segments:
                            if index+1 >= segment['start'] and index+1 <= segment['end']:
                                spagetti_list.add_spagetti(residue, index+1, index+1-int(difference_in_position),segment['color'] )
                                check = 1
                                break
                        
                        if check ==0:
                            spagetti_list.add_spagetti(residue, index+1, index+1-int(difference_in_position),'lightgrey' )
                serialized_spagettis = [spagetti.to_dict() for spagetti in spagetti_list.get_spagettis()]
                
            else:
                serialized_spagettis= None

            # for plddt coloring button
            from .utils import generate_jsmol_plddt_script
            plddt = generate_jsmol_plddt_script(sequence_file)
            results["pLDDT_Script"] = plddt


            # for distance slider button  
            from .utils import calculate_distance, extract_atom_info_and_min_distance, categorize_atoms_by_distance_adjusted
            if upload_choice == 'upload' and 'upload_file' in request.FILES:
                sequence_file_confirmed = public_path
            else: 
                sequence_file_confirmed = sequence_file
            
            if sequence_file_confirmed != None:
                atoms_info, atoms_info_DF = extract_atom_info_and_min_distance(sequence_file_confirmed, template_residue_number,Target_Chain,difference_in_position)
                other_ptms_df = pd.merge(condensed_ptm, atoms_info_DF, on="Position", how="left")
                other_ptms_df["min_distance_to_target"] = other_ptms_df["min_distance_to_variant"].fillna("NaN")
                other_ptms_df["chain_id"] = other_ptms_df["chain_id"].fillna("")
                other_ptms_df.drop(["atom_coordinates", "wt_check", "ptm_start"], axis=1, inplace=True)
                position_int =int(position)
                other_ptms_df = other_ptms_df[other_ptms_df["Position"] != position_int]
                group_1,group_2,group_3,group_4,group_5, group_6 = categorize_atoms_by_distance_adjusted(atoms_info, difference_in_position)
                
                results["atoms_info"] = atoms_info
                results["other_ptms_df"] = other_ptms_df
                results["G1"] = json.dumps(group_1)
                results["G2"] = json.dumps(group_2)
                results["G3"] = json.dumps(group_3)
                results["G4"] = json.dumps(group_4)
                results["G5"] = json.dumps(group_5)
                results["G6"] = json.dumps(group_6)
            else:
                results["other_ptms_df"] = condensed_ptm
            #-----------------------------------------------------------------------------------------------

            





            # Variant Query at target Position ------
            # Query to get all variants at user input position
            variantAtPosQueryRes = query_variantAtPos(input_data)
            print('variantAtPosQueryRes')
            print(variantAtPosQueryRes)

            from .utils import filter_and_add_dbsnp
            # List of dbs to retain
            allowed_dbs = {"dbSNP", "ClinGen", "gnomAD", "UniProt", "1000Genomes"}
            if variantAtPosQueryRes != None:
                results["variant_target"] = variantAtPosQueryRes
                # Apply filtering
                filtered_variantAtPosQueryRes = filter_and_add_dbsnp(variantAtPosQueryRes, allowed_dbs)
            if not is_ptmAll_empty:
                ptm_all_json = condensed_ptm.to_json(orient='records')
                request.session['ptm_all'] = (ptm_all_json)
            else:
                ptm_all_json = ptm_all.to_json(orient='records')
                request.session['ptm_all'] = (ptm_all_json)
            
            additional_ptmText_json = additional_ptm.to_json(orient='records')
            request.session['additional_ptmText'] = (additional_ptmText_json)
            
            # variant_target_json=variantAtPosQueryRes.to_json(orient='records')
            request.session['variant_target'] = (variantAtPosQueryRes)







            # PTM and Variant Co-Site #####-----
            if not is_ptmAll_empty:
                ptmVarCositeRes = pd.merge(condensed_ptm, variant_allDF, left_on='Position', right_on='position', how='inner')
                ptmVarCositeRes_condensed = ptmVarCositeRes.groupby(['Position', 'wt_residue', 'mutant', 'ptm_type'])['snp_id'].unique().reset_index()
                is_ptmVarCositeRes_df_empty = ptmVarCositeRes_condensed.empty
                results["ptmVarCosite"] = ptmVarCositeRes_condensed
                results['is_ptmVarCosite_df_empty'] = is_ptmVarCositeRes_df_empty
            if upload_choice == 'upload' and 'upload_file' in request.FILES:
                results["pdb_file"]= pdb_file_path2
                results["istemplateuploaded"] = True
            else:
                results["pdb_file"]= pdb_file

            print("ptmVarCosite")
            print(ptmVarCositeRes_condensed)


            





            # check differences between alphafold sequence and uniprot sequence for warning messages
            from .utils import check_differences
            warning_message = [] 
            
            # for testing if want to allow user to upload own pdb file:
            # serialized_spagettis =  
            # print('serialized_spagettis')
            # print(serialized_spagettis)

            if serialized_spagettis !=None:
                residue_by_position = {entry['position']: entry['residue'] for entry in serialized_spagettis}
                column_names_warning_df = ['Table', 'position', 'wt', 'residue', 'status']

                # Create an empty dataframe with the given column names
                warning_df = pd.DataFrame(columns=column_names_warning_df)
                
                # Add messages for variantAtPosQueryRes
                if variantAtPosQueryRes != None:
                    for variant in variantAtPosQueryRes:
                        position = variant['position']
                        wt = variant['wt']
                        residue = residue_by_position.get(position, None)
                        
                        if residue and wt != residue:
                            warning_df = warning_df.append({'Table': 'Query Variant Result', 'position': position, 'wt': wt, 'residue': residue, 'status': 'Difference found in AlphaFold Model.'}, ignore_index=True)
                        elif residue:
                            pass  # Do nothing when there is no difference
                        else:
                            warning_df = warning_df.append({'Table': 'Query Variant Result', 'position': position, 'wt': wt, 'residue': 'None', 'status': 'Position not found in AlphaFold Model.'}, ignore_index=True)
                    
                # Add messages for condensed_ptm_target
                if not condensed_ptm_target.empty:
                    differences_condensed_ptm_target = check_differences(condensed_ptm_target, residue_by_position)
                    for _, row in differences_condensed_ptm_target.iterrows():
                        warning_df = warning_df.append({'Table': 'Query PTM Results', 'position': row['Position'], 'wt': row['wt_residue'], 'residue': row['residue'], 'status': row['status']}, ignore_index=True)
                
                # Add messages for ptmVarCositeRes_condensed
                if not ptmVarCositeRes_condensed.empty:
                    differences_ptmVarCositeRes_condensed = check_differences(ptmVarCositeRes_condensed, residue_by_position)
                    for _, row in differences_ptmVarCositeRes_condensed.iterrows():
                        warning_df = warning_df.append({'Table': 'PTM & Variant Co-sites', 'position': row['Position'], 'wt': row['wt_residue'], 'residue': row['residue'], 'status': row['status']}, ignore_index=True)
                
                # Add messages for other_ptms_df
                if is_ptmAll_empty == False:
                    differences_other_ptms_df = check_differences(condensed_ptm, residue_by_position)
                    print(differences_other_ptms_df)
                    for _, row in differences_other_ptms_df.iterrows():
                        warning_df = warning_df.append({'Table': ' Additional PTMs', 'position': row['Position'], 'wt': row['wt_residue'], 'residue': row['residue'], 'status': row['status']}, ignore_index=True)

            # Combine the warning messages into a single string if desired
            merged_warning_df = (
                warning_df.groupby(['position', 'wt', 'residue', 'status'], as_index=False)
                .agg({'Table': lambda x: ', '.join(sorted(set(x)))})
            )
            # Start with a list instead of a string
            warning_message = []

            print("warning:")
            for index, row in merged_warning_df.iterrows():
                # Append formatted message to the list
                warning_message.append(
                    f"For table(s): {row['Table']}, at Position {row['position']}, UniProt reports WT = {row['wt']}, while AlphaFold reports Residue = {row['residue']}. Status: {row['status']}"
                )
                print(f"Row {index}: {row.to_dict()}")

            # Join all messages into a single string for display, if needed
            final_warning_message = "\n".join(warning_message)
            print(final_warning_message)









            # For downloading
            querysummary_download = {'Column1': ['uniprot_id', 'gene_name', 'residue_number', 'protein_name'], 'Column2': [uniprot_id, gene_name, uniprot_residue_number, protein_name]}
            request.session['querysummary_download'] = querysummary_download
            request.session['condensed_ptm_target'] = condensed_ptm_target.to_dict()
            ptmVarCosite_json = ptmVarCositeRes_condensed.copy()
            ptmVarCosite_json['snp_id'] = ptmVarCosite_json['snp_id'].apply(lambda x: str(x))
            request.session['ptmVarCositeRes_condensed'] = ptmVarCosite_json.to_dict()
            other_ptms_df_json = other_ptms_df.copy()
            other_ptms_df_json['evidence_code'] = other_ptms_df_json['evidence_code'].apply(lambda x: str(x))
            request.session['other_ptms_df'] = other_ptms_df_json.to_dict()
            atoms_info_json = atoms_info_DF.copy()
            atoms_info_json['atom_coordinates'] = atoms_info_json['atom_coordinates'].apply(lambda x: str(x))
            request.session['atoms_info_json'] = atoms_info_json.to_dict()
            request.session['serialized_spagettis'] = serialized_spagettis

            results.update({"uniprot_id": uniprot_id,
                       "gene_name" : gene_name,
                       "wt_residue": input_data["user_wt"],
                       "alt_gene_names": alt_gene_names,
                       "proteinName": protein_name,
                       "additional_ptmText": additional_ptm,
                       "is_additionalPTM_df_empty": is_additionalPTM_df_empty,
                       'highlighting': highlight,
                       'variant_all': variant_all,
                       'topologyQueryRes': topologyQueryRes,
                       'spagettis': serialized_spagettis,
                       'sequence':sequence,
                       'segments': serialized_segments,
                       'Target_Chain':Target_Chain,
                        'warning_message_combined': final_warning_message,
                       })                      
            return render(request, "dbQuery/results_v2.html", context = results)  
    else:
        form = QueryForm()

    return render(request, 'dbQuery/home.html', {'form': form})







def home(request):
    form = QueryForm()
    return render(request, 'dbQuery/home.html', {"form": form})

def base(request):
    return render(request, 'dbQuery/base.html')

def documentation(request):
    return render(request, 'dbQuery/documentation.html')

def dataset(request):
    return render(request, 'dbQuery/dataset.html')

def contact(request):
    return render(request, 'dbQuery/contact.html')

def formFill(request):
    return HttpResponse("Function executed successfully")

def download_proteins(request):
    return FileResponse(open('/project/data/VariantP/variant_data/canonical_list_website.txt', 'rb'), as_attachment=True, filename='canonical_proteins.txt')

from django.http import HttpResponse
from openpyxl import Workbook
from io import BytesIO

def download_excel(request):
    # Example data
    querysummary_download = request.session.get('querysummary_download', [])
    variant_target = request.session.get('variant_target', [])
    condensed_ptm_target = request.session.get('condensed_ptm_target', [])
    ptmVarCositeRes_condensed = request.session.get('ptmVarCositeRes_condensed', [])
    other_ptms_df = request.session.get('other_ptms_df', [])
    
    file_naming =querysummary_download['Column2'][0] + '_' +str(querysummary_download['Column2'][2]) 
    
   # Flatten variant_target (excluding some columns)
    flattened_variant_target = []
    for variant in variant_target:
        for info in variant.get('variantInfo', []):
            combined = {
                'position': variant.get('position'),
                'wt': variant.get('wt'),
                'mutant': variant.get('mutant'),
                'wt_check': variant.get('wt_check'),
                'accession_id': info.get('snp_id'),
                'db': ', '.join(info.get('db', [])),
                'clinical_significance': ', '.join(info.get('clinical_significance', [])),
                'conditions': ', '.join(info.get('conditions', [])),
                'disease': ', '.join(info.get('disease', [])),
                'review_status': ', '.join(info.get('review_status', [])),
            }
            flattened_variant_target.append(combined)

    
    # Data for separate sheets
    datasets = {
        'Query Summary': pd.DataFrame(querysummary_download),
        'Query Variant Result': pd.DataFrame(flattened_variant_target),
        'Query PTM Results': pd.DataFrame(condensed_ptm_target),
        'PTM & Variant Co-sites': pd.DataFrame(ptmVarCositeRes_condensed),
        'Additional PTMs': pd.DataFrame(other_ptms_df),
    }
    # Create an Excel writer
    with BytesIO() as buffer:
        with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
            # Write each dataset to a separate sheet
            for sheet_name, df in datasets.items():
                df.to_excel(writer, sheet_name=sheet_name, index=False)
        
        buffer.seek(0)
        response = HttpResponse(
            buffer,
            content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
        )
        response['Content-Disposition'] = f'attachment; filename="Query_Result_{file_naming}.xlsx"'
        return response



def home2(request):
    return render(request, 'dbQuery/home_copy.html')




