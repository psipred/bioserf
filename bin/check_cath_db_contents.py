import glob
import re

cath_id_re = re.compile("(\d.{3}\w\d{2})")
for name in glob.glob('*.py'):
    with open(name) as fh:
        contents = fh.readlines()
        for line in contents:
            if "knowns" in line:
                print(line)
                result = cath_id_re.search(line)
                if result:
                    for cath_id in result.groups():
                        print(cath_id)
