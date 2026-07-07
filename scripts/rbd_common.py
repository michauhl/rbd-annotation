from typing import Optional
import requests


################################################################################

def read_fasta_file(fasta_file):
    """
    Read a FASTA file into a dictionary with sequence ID -> sequence.

    The sequence ID is the first whitespace-delimited token after ">".
    """
    seqs_dic = {}
    seq_id = None
    seq_chunks = []

    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_id is not None:
                    assert seq_id not in seqs_dic, f"Duplicate FASTA sequence ID: {seq_id}"
                    seqs_dic[seq_id] = "".join(seq_chunks)
                seq_id = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)

    if seq_id is not None:
        assert seq_id not in seqs_dic, f"Duplicate FASTA sequence ID: {seq_id}"
        seqs_dic[seq_id] = "".join(seq_chunks)

    return seqs_dic


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

def refine_intervals(intervals):
    """
    Filter list of start, end pairs to remove intervals that are fully 
    contained within another interval.


    >>> l = [[1, 10], [2, 5], [3, 6], [8, 9], [12, 15], [13, 14]]
    >>> refine_intervals(l)
    [[1, 10], [12, 15]]
    >>> l = []
    >>> refine_intervals(l)
    []
    >>> l = [[1, 10]]
    >>> refine_intervals(l)
    [[1, 10]]

    """

    # Sort by start coordinate, then by end coordinate in descending order.
    intervals.sort(key=lambda x: (x[0], -x[1]))

    refined = []
    for interval in intervals:
        # Add to refined list if it's not fully contained within the last added interval.
        if not refined or interval[1] > refined[-1][1]:
            refined.append(interval)
    
    return refined


################################################################################

def merge_intervals(intervals, overlap_threshold=0.6):

    """
    Merge overlapping intervals, and let the new start and end positions be the 
    minimum and maximum start and end positions of the overlapping intervals, 
    respectively.

    >>> intervals = [[100, 200], [120, 220]]
    >>> merge_intervals(intervals)
    [[100, 220]]
    >>> merge_intervals(intervals, overlap_threshold=0.9)
    [[100, 200], [120, 220]]
    >>> intervals = [[100, 200], [200, 300]]
    >>> merge_intervals(intervals)
    [[100, 200], [200, 300]]
    >>> intervals = [[100, 200]]
    >>> merge_intervals(intervals)
    [[100, 200]]
    >>> intervals = []
    >>> merge_intervals(intervals)
    []

    """

    # Sort intervals by start position.
    intervals.sort()
    
    merged = []
    for start, end in intervals:
        if not merged:
            merged.append([start, end])
        else:
            last_start, last_end = merged[-1]
            
            # Compute overlap
            overlap = max(0, min(end, last_end) - max(start, last_start))
            interval_length = min(end - start, last_end - last_start)  # The smaller interval
            
            if overlap / interval_length >= overlap_threshold:  # If overlap is >= 60%
                # Merge intervals
                merged[-1] = [min(last_start, start), max(last_end, end)]
            else:
                merged.append([start, end])
    
    return merged


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


