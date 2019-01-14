import sys
import os
import xml.etree.ElementTree as ET
import json
from json import JSONEncoder


class MarkedList:
    list = None

    def __init__(self, l):
        self.list = l


class CustomJSONEncoder(JSONEncoder):
    def default(self, o):
        if isinstance(o, MarkedList):
            return f'##<{o.list}>##'

        return json.JSONEncoder.default(self, o)


if len(sys.argv) != 3:
    print('Incorrect number of arguments', file=sys.stderr)
    sys.exit(1)

xmlfile = sys.argv[1]
outdir = sys.argv[2]

tree = ET.parse(xmlfile)
root = tree.getroot()

for idx, s in enumerate(root):
    data = {}
    metadata = {}
    common = {'names': MarkedList(['kappa'])}
    atom = {'names': MarkedList(['A', 'B'])}
    name = s.attrib['Name']

    metadata['method'] = 'EEM'
    metadata['name'] = name
    metadata['publication'] = ''

    data['metadata'] = metadata
    data['common'] = common
    data['atom'] = atom

    uc = s.find('UnitConversion')
    kappa_factor = float(uc.attrib['KappaFactor'])
    ab_factor = float(uc.attrib['ABFactor'])

    parameters = s.find('Parameters')
    common['values'] = MarkedList([float(parameters.attrib['Kappa']) * kappa_factor])

    all_values = []
    for e in parameters.findall('Element'):
        element_name = e.attrib['Name']
        bond_info = []
        for b in e.findall('Bond'):
            bond_info.append((b.attrib['Type'], float(b.attrib['A']), float(b.attrib['B'])))

        for t, a, b in bond_info:
            all_values.append(
                {'key': MarkedList([element_name, "hbo", t]), 'value': MarkedList([a * ab_factor, b * ab_factor])})

    atom['data'] = all_values

    with open(os.path.join(outdir, f'EEM_{idx:02d}_{name}.json'), 'w') as f:
        output = json.dumps(data, indent=2, cls=CustomJSONEncoder)
        output = output.replace('"##<', "").replace('>##"', "").replace('\'', '"')
        f.write(output)
