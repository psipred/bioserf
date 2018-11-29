import glob
import re
import shutil

pdb_dir = '/scratch0/NOT_BACKED_UP/dbuchan/genome3d/home/camp/buchand/' \
          'working/genome3d/parse_cath_domth_output/*.pdb'

# REMARK GENOME3D UNIPROT_ID sp|Q6UWF7|NXPE4
# REMARK GENOME3D UNIPROT_MD5 .0d7e5e18b2d0740b3737f2376ad64e3d
id_pattern = re.compile("REMARK GENOME3D UNIPROT_ID .{2}\|(.+)\|.+")
md_pattern = re.compile("REMARK GENOME3D UNIPROT_MD5 \.(.+)")

for file in glob.glob(pdb_dir):
    print(file)
    output_str = ''
    with open(file) as pdb_file:
        for line in pdb_file:
            id_result = re.match(id_pattern, line)
            md_result = re.match(md_pattern, line)
            if(id_result):
                line = "REMARK GENOME3D UNIPROT_ID "+id_result.group(1)+"\n"
            if(md_result):
                line = "REMARK GENOME3D UNIPROT_MD5 "+md_result.group(1)+"\n"
            output_str += line
    output_file = open('tmp.pdb', 'w')
    output_file.write(output_str)
    shutil.copy2('tmp.pdb', file)
