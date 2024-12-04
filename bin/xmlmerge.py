#!/usr/bin/env python
import xml.etree.ElementTree as ET
import sys

def merge_xml_files(input_files, output_file):
    if not input_files:
        print("No XML files provided.")
        return

    # Parse the first XML file and get the root element
    first_file = input_files[0]
    tree = ET.parse(first_file)
    root = tree.getroot()

    # Iterate over the remaining XML files and append their content to the root
    for xml_file in input_files[1:]:
        tree = ET.parse(xml_file)
        root.extend(tree.getroot())

    # Write the merged content to the output file
    tree = ET.ElementTree(root)
    tree.write(output_file, encoding='utf-8', xml_declaration=True)
    print(f"Merged XML file saved as {output_file}")

if __name__ == "__main__":
    input_files = sys.argv[1:-1]
    output_file = sys.argv[-1]
    merge_xml_files(input_files, output_file)