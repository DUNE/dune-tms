import xml.etree.ElementTree as ET

src = "nd_hall_with_lar_tms_nosand_CLEAN_NEW_GEOMETRY.gdml"
dst = "nd_hall_80layers_MAPPED.gdml"
mapper = "Mapper_global_DIAGNOSTIC.txt"

steel_volumes = {
    "thinvolTMS",
    "thinvol2TMS",
    "thickvolTMS",
    "thickvol2TMS",
    "doublevolTMS",
    "doublevol2TMS",
}

tree = ET.parse(src)
root = tree.getroot()

found_tms = False
removed = 0

for volume in root.iter():
    if volume.tag.split("}")[-1] != "volume":
        continue

    name = volume.get("name")

    # Remove BField/ArbBField from the six reused steel logical volumes
    # and from volTMS before inserting the one parent field.
    if name in steel_volumes or name == "volTMS":
        for child in list(volume):
            if child.tag.split("}")[-1] != "auxiliary":
                continue

            field_type = child.get("auxtype") or child.get("auxType")

            if field_type in {"BField", "ArbBField"}:
                volume.remove(child)
                removed += 1

    if name == "volTMS":
        found_tms = True

        namespace = ""
        if volume.tag.startswith("{"):
            namespace = volume.tag.split("}", 1)[0] + "}"

        ET.SubElement(
            volume,
            namespace + "auxiliary",
            {
                "auxtype": "ArbBField",
                "auxvalue": mapper,
            },
        )

if not found_tms:
    raise RuntimeError("volTMS was not found")

tree.write(dst, encoding="utf-8", xml_declaration=True)

print("Removed old field auxiliaries:", removed)
print("Inserted ONE ArbBField on volTMS")
print("Output:", dst)
