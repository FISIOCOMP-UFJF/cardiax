import xml.etree.ElementTree as ET
from collections import defaultdict

# ============================================================
# Arquivos
# ============================================================

ARQUIVO_ENTRADA = "malha_antiga.xml"
ARQUIVO_SAIDA = "malha_nova.xml"

# ============================================================
# Leitura do XML
# ============================================================

with open(ARQUIVO_ENTRADA, "r", encoding="utf-8") as f:
    conteudo = f.read()

root = ET.fromstring(f"<root>{conteudo}</root>")

mesh = root.find("mesh")

nodes_tag = mesh.find("nodes")
elements_tag = mesh.find("elements")
element_data_tag = mesh.find("element_data")

if element_data_tag is None:
    raise RuntimeError("A tag <element_data> não foi encontrada.")

# ============================================================
# Acumuladores
# ============================================================

soma_lat = defaultdict(float)
contagem = defaultdict(int)

# ============================================================
# Percorre elementos
# ============================================================

for elemento in element_data_tag.findall("element"):

    eid = int(elemento.get("id"))

    eikonal_tag = elemento.find("eikonal")

    if eikonal_tag is None:
        continue

    lat = float(eikonal_tag.text)

    elem_xml = elements_tag.find(f"./element[@id='{eid}']")

    if elem_xml is None:
        continue

    nos = [
        int(elem_xml.get("v0")),
        int(elem_xml.get("v1")),
        int(elem_xml.get("v2")),
        int(elem_xml.get("v3")),
    ]

    for nid in nos:
        soma_lat[nid] += lat
        contagem[nid] += 1

# ============================================================
# Remove eikonal antigo dos elementos
# ============================================================

for elemento in element_data_tag.findall("element"):

    eik = elemento.find("eikonal")

    if eik is not None:
        elemento.remove(eik)

# ============================================================
# Cria nova seção <eikonal>
# ============================================================

eikonal_xml = ET.Element("eikonal")
eikonal_xml.set("size", nodes_tag.get("size"))

n_nos = int(nodes_tag.get("size"))

for nid in range(n_nos):

    if contagem[nid] > 0:
        lat = soma_lat[nid] / contagem[nid]
    else:
        lat = 0.0

    node = ET.SubElement(eikonal_xml, "node")
    node.set("id", str(nid))
    node.set("lat", f"{lat:.6f}")

# ============================================================
# Insere após <element_data>
# ============================================================

indice = list(mesh).index(element_data_tag)

mesh.insert(indice + 1, eikonal_xml)

# ============================================================
# Salva mantendo elasticity e mesh
# ============================================================

with open(ARQUIVO_SAIDA, "w", encoding="utf-8") as f:

    for child in root:
        f.write(ET.tostring(child, encoding="unicode"))
        f.write("\n")

print("Conversão concluída!")
print(f"Arquivo salvo em: {ARQUIVO_SAIDA}")
