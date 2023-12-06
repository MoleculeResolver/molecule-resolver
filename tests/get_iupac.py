import json
import requests
from lxml import etree
from copy import deepcopy
from tqdm import tqdm

def get_iupac(smiles: str):
    html_doc = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/" + smiles + "/record/XML")
    if html_doc.status_code != 200:
        return 
    root = etree.XML(html_doc.text)

    iupac_elements = root.findall(".//{*}PC-Urn_label")
    for e in iupac_elements:
        if "IUPAC Name" == e.text:
            urn = e.getparent()
            iupac_name_type = urn.find(".//{*}PC-Urn_name").text

            info_data = urn.getparent().getparent()
            iupac_name = info_data.find(".//{*}PC-InfoData_value_sval").text

            if iupac_name_type == "Preferred":
                return iupac_name

def main():     
    with open("benchmark_component_molecules.json", "r") as f:
        benchmark = json.load(f)

    new_benchmark = {}
    for name, data in tqdm(benchmark.items()):
        iupac_name = get_iupac(data["SMILES"])
        if iupac_name is not None:
            data["iupac_name"] = iupac_name
            new_benchmark[name] = data

    with open("benchmark_component_molecules_iupac.json", "w") as f:
        json.dump(new_benchmark, f, indent=4)

if __name__ == "__main__":
    main()

