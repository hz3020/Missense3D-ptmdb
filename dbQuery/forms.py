from django import forms

class QueryForm(forms.Form):
    UNIPROT = forms.CharField(label="UniProt ID", min_length=6, required=False)
    GENENAME = forms.CharField(label="Gene Name", required=False)
    POSITION = forms.CharField(label="Residue Number", required=True)
    RESIDUE = forms.CharField(label="Wild-type amino acid", required=True)
    upload_choice = forms.ChoiceField(
        choices=[('upload', 'Upload a file'), ('no_upload', 'Proceed without upload')],
        widget=forms.RadioSelect,
        label="Do you want to upload a file?",
        required=True,
        initial='no_upload'
    )
    upload_file = forms.FileField(
        required=False,  # Only required if 'upload' is selected
        label="Upload a PDB file (optional)"
    )
    start_position = forms.IntegerField(
        required=False,  # Only required if 'upload' is selected
        label="Start Position",
        min_value=1
    )
    Target_Chain = forms.CharField(label="Target Chain", required=False)

    def clean(self):
        cleaned_data = super().clean()
        uniprot_id = cleaned_data.get("UNIPROT")
        gene_name = cleaned_data.get("GENENAME")

        if gene_name and uniprot_id:
            raise forms.ValidationError("Please enter either a UniProt ID or a gene name, not both.")
        elif not gene_name and not uniprot_id:
            raise forms.ValidationError("Please enter a UniProt ID or a gene name.")
        
        upload_choice = cleaned_data.get('upload_choice')
        upload_file = cleaned_data.get('upload_file')
        start_position = cleaned_data.get('start_position')
        Target_Chain = cleaned_data.get('Target_Chain')

        # Validate that start_position is provided if the user chooses to upload
        if upload_choice == 'upload':
            if not upload_file:
                raise forms.ValidationError("Please upload a file.")
            if start_position is None:
                raise forms.ValidationError("Please provide a start position for the uploaded file.")
            if Target_Chain is None:
                raise forms.ValidationError("Please provide a Target_Chain for the uploaded file.")
        return cleaned_data


        




        