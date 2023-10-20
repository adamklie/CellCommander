from PyInquirer import prompt


def check_marker_genes(adata, marker_genes):
    if isinstance(marker_genes, dict):
        # Check markers against adata.var
        marker_genes_in_data = dict()
        for ct, markers in marker_genes.items():
            markers_found = list()
            for marker in markers:
                if marker in adata.var.index:
                    markers_found.append(marker)
            marker_genes_in_data[ct] = markers_found

        # Remove any keys with lists of length 0
        marker_genes_in_data = {k: v for k, v in marker_genes_in_data.items() if len(v) > 0}
        
        return marker_genes_in_data
    

def get_user_input(keys, choices):
    user_data = {}
    for key in keys:
        questions = [
            {
                'type': 'list',
                'name': key,
                'message': f'Select value for cluster {key}',
                'choices': choices,
            }
        ]
        answer = prompt(questions)
        user_data[key] = answer[key]
    return user_data