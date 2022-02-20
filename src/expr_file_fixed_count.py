'''
expr_file_fixed_count.py
Generate expression file for simulation

@author Jorge Mestre Tomas (jormart2@alumni.uv.es)
@date 19/01/2020
'''

def create_expr_file(f_cat: str, f_del: str, n_trans: int, coverage: int, output: str):

    deleted = set()
    with open(f_del, 'r') as del_in:
        skip = del_in.readline()
        for line in del_in:
            line = line.split()
            deleted.add(line[0])
    del_in.close()


    tot_trans = deleted
    with open(f_cat, 'r') as cat_in:
        skip = cat_in.readline()
        for line in cat_in:
            trans = line.split()[0]
            if trans not in deleted:
                tot_trans.add(trans)
            if len(tot_trans) == n_trans:
                break

    cat_in.close()
    
    f_out = open(output, 'w')
    f_out.write('target_id\test_counts\ttpm\n')
    for trans in tot_trans:
        f_out.write(trans + '\t' + coverage + '\t' + coverage + '\n')
    f_out.close()