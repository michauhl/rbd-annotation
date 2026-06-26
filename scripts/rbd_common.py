from typing import Optional
import requests


################################################################################

def fetch_interpro_json(interpro_id,
                        base_url="https://www.ebi.ac.uk/interpro/api/entry/interpro/"):

    """
    https://www.ebi.ac.uk/interpro/api/entry/interpro/IPR000504?format=json
    
    """
    url = f"{base_url}{interpro_id}?format=json"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        print(f"Failed to fetch data for {interpro_id} (HTTP {response.status_code})")
        return None


################################################################################

def get_interpro_parent_infos(interpro_ids_dic):
    """
    Given dictionary with Interpro IDs, try to get parent relationship, i.e., 
    get the parent Interpro ID + name. This could be the same ID, or a different one.
    
    """
    parent_infos_dic = {}
    
    for interpro_id in interpro_ids_dic:

        data = fetch_interpro_json(interpro_id)

        if not data:
            print("WARNING: Failed to fetch data for", interpro_id)
            continue  # Skip if data could not be fetched

        # Store data["metadata"]["hierarchy"]["accession"], if exists, else store None.
        parent_interpro_id = data.get("metadata", {}).get("hierarchy", {}).get("accession", "-")
        interpro_name = data.get("metadata", {}).get("name", {}).get("name", "-")
        interpro_type = data.get("metadata", {}).get("type", "-")
        interpro_short = data.get("metadata", {}).get("name", {}).get("short", "-")

        rbd_infos = InterproRBDInfo(interpro_id, interpro_name, interpro_type, interpro_short,
                                    interpro_id, interpro_name, interpro_type, interpro_short)

        if interpro_id != parent_interpro_id:
            parent_data = fetch_interpro_json(parent_interpro_id)
            if parent_data:
                parent_interpro_name = parent_data.get("metadata", {}).get("name", {}).get("name", "-")
                parent_interpro_type = parent_data.get("metadata", {}).get("type", "-")
                parent_interpro_short = parent_data.get("metadata", {}).get("name", {}).get("short", "-")
                rbd_infos.parent_id = parent_interpro_id
                rbd_infos.parent_name = parent_interpro_name
                rbd_infos.parent_type = parent_interpro_type
                rbd_infos.parent_short = parent_interpro_short

        parent_infos_dic[interpro_id] = rbd_infos

    return parent_infos_dic


################################################################################

class InterproRBDInfo:
    """
    Stores Interpro RBD (RNA-binding domain) info, including parent information.  

    """
    def __init__(self,
                 interpro_id: str,
                 interpro_name: str,
                 interpro_type: str,
                 interpro_short: str,
                 parent_id: Optional[str] = None,
                 parent_name: Optional[str] = None,
                 parent_type: Optional[str] = None,
                 parent_short: Optional[str] = None) -> None:
        
        self.interpro_id = interpro_id
        self.interpro_name = interpro_name
        self.interpro_type = interpro_type
        self.interpro_short = interpro_short
        self.parent_id = parent_id
        self.parent_name = parent_name
        self.parent_type = parent_type
        self.parent_short = parent_short


################################################################################



