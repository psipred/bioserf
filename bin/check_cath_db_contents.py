import glob
import re

cath_id_re = re.compile("\d.{3}\w\d{2}")
for name in glob.glob('*.py'):
    print(name)
    with open(name) as fh:
        contents = fh.readlines()
        result = cath_id_re.match(contents)
        if result:
            for match in result:
                print(match)
