import glob
import re
import os
import urllib.request
import urllib.urlretrieve

cath_id_re = re.compile("(\d.{3}\w\d{2})")
for name in glob.glob('*.py'):
    with open(name) as fh:
        contents = fh.readlines()
        for line in contents:
            if "knowns" in line:
                result = cath_id_re.search(line)
                if result:
                    for cath_id in result.groups():
                        if not os.path.isfile("/data/cath_data/dompdb"+cath_id):
                            urllib.urlretrieve ("http://www.cathdb.info/version/v4_2_0/api/rest/id/"+cath_id+".pdb", cath_id)
