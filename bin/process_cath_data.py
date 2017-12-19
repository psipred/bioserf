import sys
import glob

# we read in the cath domain list and the cath domall file and produce
# an annotated version of the domain list
#
# For the sake of sanity we only handle contiguous chains whose starts
# as positive integers


class Vividict(dict):
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value


def read_domall(domall_file):
    domains = Vividict()
    with open(domall_file, 'r') as dfh:
        for line in dfh:
            domain_count = 0
            line = line.rstrip()
            if line.startswith('#'):
                continue
            else:
                # print(line)
                leader = line[0:15]
                data = line[15:]
                leader_entries = leader.split()
                total_domains = int(leader_entries[1].lstrip('D'))
                while domain_count != total_domains:
                    seg_number = int(data[0])
                    entry_length = (seg_number*19)+3
                    domain_entry = data[0:entry_length]
                    data = data[entry_length:]
                    # print(domain_entry)
                    seg_count = 0
                    total_segs = int(domain_entry[0])
                    domain_count += 1
                    while seg_count != total_segs:
                        seg_data = domain_entry[1:]
                        start = int(seg_data[3+(seg_count*19):9+(seg_count*19)])
                        stop = int(seg_data[12+(seg_count*19):17+(seg_count*19)])
                        seg_count += 1
                        if start < 1:
                            start = 1
                        domains[leader_entries[0]][domain_count][seg_count]['start']=start
                        domains[leader_entries[0]][domain_count][seg_count]['stop']=stop
                        domains[leader_entries[0]][domain_count][seg_count]['length']=stop-start+1


    return(domains)


def read_domain_list(domall, domain_file):
    with open(domain_file, 'r') as dfh:
        for line in dfh:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            entries = line.split()
            # if entries[0][0:5] != '4q43F':
            #     continue
            if entries[0][0:5] in domall.keys():
                # print(domall[entries[0][0:5]])
                if entries[0][5:7] == '00':
                    if len(domall[entries[0][0:5]][1]) == 1:
                        print(line+" " +
                              str(domall[entries[0][0:5]][1][1]['start'])+" " +
                              str(domall[entries[0][0:5]][1][1]['stop'])+" " +
                              str(domall[entries[0][0:5]][1][1]['length'])
                              )
                    # handle single domain chain
                else:
                    i = int(entries[0][5:7])
                    if len(domall[entries[0][0:5]][i]) == 1:
                        print(line+" " +
                              str(domall[entries[0][0:5]][i][1]['start'])+" " +
                              str(domall[entries[0][0:5]][i][1]['stop'])+" " +
                              str(domall[entries[0][0:5]][1][1]['length'])
                              )


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i


def read_tdb(domall, tdb_list):
    tdb_files = glob.glob(tdb_list+"/*.tdb")
    for tdb in tdb_files:
        domid = tdb[-11:-4]
        chainid = tdb[-11:-6]
        domain_number = int(tdb[-6:-4])
        domain_len = file_len(tdb)
        if domain_number == 0:
            domain_number == 1
        # print(domid, chainid, domain_number)
        if chainid in domall:
            if domain_number in domall[chainid]:
                pass
            else:
                domall[chainid][domain_number][1]['start'] = '-'
                domall[chainid][domain_number][1]['stop'] = '-'
                domall[chainid][domain_number][1]['length'] = domain_len


domain_list = sys.argv[1]
domall_list = sys.argv[2]
tdb_list = sys.argv[3]

domall = read_domall(domall_list)
#print(domall)
domall = read_tdb(domall, tdb_list)
print(domall)
# read_domain_list(domall, domain_list)
