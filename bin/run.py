# %%
import sys
import os
import subprocess
import yaml

def main():
    inp = os.path.dirname(__file__) + '/config.yaml'
    default_config = ''
    help = f'''
        One command executing mode

        Usage: 
            phasiHunter run [-i] [config file]
            phasiHunter run -d

        option:
            -i: yaml format config file
            -d: using the default config, defalut config file is {inp}
            -h: print help information
        
        WARNIG: make sure choose the correct config file before run this command
    '''

    version = '''
    '''

    if len(sys.argv) == 1:
        print(help)
        sys.exit()

    for i in range(1, len(sys.argv)):
        if sys.argv[i] == '-i':
            inp = sys.argv[i+1]
        elif sys.argv[i] == '-d':
            default_config = inp
        elif sys.argv[i] == '-h':
            print(help)
            sys.exit()

    print('Parse yaml file')
    if default_config != '':
        with open(default_config, 'r') as fn:
            data = yaml.load(fn, Loader=yaml.FullLoader)
    else:
        with open(inp, 'r') as fn:
            data = yaml.load(fn, Loader=yaml.FullLoader)
        
    Runing_module = data['Runing_module']

    module_list = list()
    for k,v in Runing_module.items():
        if type(v) != dict:
            if Runing_module[k] == 'y':
                module_list.append(k)
        if type(v) == dict:
            for i in Runing_module[k]:
                if Runing_module[k][i] == 'y':
                    if k in module_list:
                        pass
                    else:
                        module_list.append(k)
                    module_list.append(i)

    string = ''
    for i in (module_list):
        if i == 'target' or i == 'initiator' or i== 'deg' or i == 'phasiRNA_target' or i == 'phasiRNA_deg':
            string += '    '
            string += i + '\n'
        else:
            string += '  '
            string += i + '\n'

    print('Runing Module:  ')
    print(string)


    if 'preprocess' in module_list:
        Preprocess(data)
    if 'phase' in module_list:
        Phase(data)
    if 'integration' in module_list:
        Integration(data)
    if 'visualization' in module_list:
        Visualization(data)
    if 'initiator_prediction_and_verification' in module_list and 'target' in module_list:
        Target(data)
    if 'initiator_prediction_and_verification' in module_list and 'initiator' in module_list:
        Initiator(data)
    if 'initiator_prediction_and_verification' in module_list and 'deg' in module_list:
        Deg(data)
    if 'phasiRNA_target_prediction_and_verification' in module_list and 'phasiRNA_target' in module_list:
        phasiRNA_target(data)
    if 'phasiRNA_target_prediction_and_verification' in module_list and 'phasiRNA_deg' in module_list:
        phasiRNA_deg(data)
    print('Analysis finished')


def Preprocess(data):
    print('-----------------------------------------------------')
    print('Run preprocess module')
    # build index if index not exist
    print('Detect index ...')
    index_list = []
    preprocess_subprocess = []
    if data['preprocess']['index'] != None and len(data['preprocess']['index']) > 0:
        print('Index exsits')
        if data['preprocess']['reference_fasta'] != None and len(data['preprocess']['reference_fasta']) > 0:
            print('Reference sequence also exsit, which will be ignored')
        for i in data['preprocess']['index']:
            index_list.append(i)
    else:
        print('No index exist, build index with reference_fasta')
        cmd = f"[ -e index ] || mkdir index"
        os.system(cmd)
        for i in range(0, len(data['preprocess']['reference_fasta'])):
            bindex=f"index/{data['preprocess']['reference_fasta'][i].split('.')[0].split('/')[-1]}"
            index_list.append(bindex)
            cmd = f"bowtie-build {data['preprocess']['reference_fasta'][i]} {bindex}"
            print(cmd)
            process = subprocess.Popen(cmd, shell=True)
            preprocess_subprocess.append(process)

        for process in preprocess_subprocess:
            process.wait()

        all_finished = all([process.poll() is not None for process in preprocess_subprocess])

        if all_finished:
            print("All index build finished")


    cmd = f"phasiHunter preprocess -m {data['preprocess']['mode']} -i {data['preprocess']['inputfile']} -mi {data['preprocess']['minimal_sRNA_length_cutoff']} -ma {data['preprocess']['maxmial_sRNA_length_cutoff']} -e {data['preprocess']['sRNA_expression_cutoff']} -n {data['preprocess']['library_normalization_base']} -o {data['preprocess']['outfile_name'][0]} -in {index_list[0]}"
    print(cmd)
    os.system(cmd)

    for i in range(1, len(index_list)):
        cmd = f"phasiHunter preprocess -m m -i {data['preprocess']['inputfile'].split('.')[0].split('./')[-1]}_trimmed_format_filter.fa -mi {data['preprocess']['minimal_sRNA_length_cutoff']} -ma {data['preprocess']['maxmial_sRNA_length_cutoff']} -e {data['preprocess']['sRNA_expression_cutoff']} -n {data['preprocess']['library_normalization_base']} -o {data['preprocess']['outfile_name'][i]} -in {index_list[i]}"
        print(cmd)
        os.system(cmd)


def Phase(data):
    print('-----------------------------------------------------')
    print('Run phase module')
    phase_parameter = data['phase']
    _cm = phase_parameter['mapped_cdna_file']
    _c = phase_parameter['cdna_fasta']
    _gm = phase_parameter['mapped_gdna_file']
    _g = phase_parameter['gdna_fasta']
    _fm = phase_parameter['mapped_flnc_file']
    _f = phase_parameter['flnc_fasta']
    _fa = phase_parameter['sRNA_fa']
    _a = phase_parameter['allsiRNA_cluster_output']
    _o = phase_parameter['phasiRNA_cluster_output']
    _me = phase_parameter['phasiRNA_prediction_method']
    _il = phase_parameter['phasiRNA_cluster_island']
    _pl = phase_parameter['phase_length']
    _pn = phase_parameter['phase_number_cutoff']
    _mh = phase_parameter['bowtie_max_hits_cutoff']
    _j = phase_parameter['parallel_cores']
    _pv = phase_parameter['pvalue_cutoff']
    _ps = phase_parameter['phase_score_cutoff']
    _pr = phase_parameter['phase_ratio_cutoff']
    _cl = phase_parameter['delete_index']

    cmd = f"phasiHunter phase -cm {_cm} -c {_c} -gm {_gm} -g {_g} -fm {_fm} -f {_f} -fa {_fa} -a {_a} -o {_o} -me {_me} -il {_il} -pl {_pl} -pn {_pn} -mh {_mh} -j {_j} -pv {_pv} -ps {_ps} -pr {_pr} -cl {_cl}" 

    print(cmd)
    os.system(cmd)

def Integration(data):
    print('-----------------------------------------------------')
    print('Run integration module')
    integration_parameter = data['integration']

    _io = integration_parameter['o_inputfile']
    _ia = integration_parameter['a_inputfile']
    _an = integration_parameter['gff3']
    _g = integration_parameter['gdna_based_PHAS_Loci']
    _o = integration_parameter['integration_phasiRNA_cluster']
    _a = integration_parameter['integration_allsiRNA_cluster']
    _s = integration_parameter['integration_summary']
    _po = integration_parameter['integration_PHAS_Loci_info']
    _j = integration_parameter['parallel_cores']
    _pn = integration_parameter['phase_number_cutoff']
    _pl = integration_parameter['phase_length']
    _pv = integration_parameter['pvalue_cutoff']
    _il = integration_parameter['phasiRNA_cluster_island']
    _dp = integration_parameter['discard_only_P_method_result']
    _fn = integration_parameter['flnc_annotation_file']

    cmd = f"phasiHunter integration -io {_io} -ia {_ia} -an {_an} -g {_g} -o {_o} -a {_a} -s {_s} -po {_po} -j {_j} -pn {_pn} -pl {_pl} -pv {_pv} -il {_il} -dp {_dp} {_fn}" 

    print(cmd)
    os.system(cmd)

def Visualization(data):
    print('-----------------------------------------------------')
    print('Run visualization module')

    visulization_parameters = data['visulization']

    _io = visulization_parameters['o_inputfile']
    _ia = visulization_parameters['a_inputfile']
    _ip = visulization_parameters['p_inputfile']
    _a = visulization_parameters['output_alignment_file']
    _o = visulization_parameters['output_phasiRNA_fa']
    _p = visulization_parameters['output_PHAS_fa']
    _pl = visulization_parameters['phase_length']
    _m = visulization_parameters['Y_axis']
    _c = visulization_parameters['cdna_fasta']
    _g = visulization_parameters['gdna_fasta']
    _f = visulization_parameters['flnc_fasta']
    _pc = visulization_parameters['plot_cdna_based_phasiRNA_cluster']
    _pg = visulization_parameters['plot_gdna_based_phasiRNA_cluster']
    _pf = visulization_parameters['plot_flnc_based_phasiRNA_cluster']

    cmd = f"phasiHunter visulization -io {_io} -ia {_ia} -ip {_ip} -a {_a} -o {_o} -p {_p} -pl {_pl} -m {_m} -c {_c} -g {_g} -f {_f} -pc {_pc} -pg {_pg} -pf {_pf}" 

    print(cmd)
    os.system(cmd)

def Target(data):
    print('=====================================================')
    print('Initiator prediction and verification')

    print('-----------------------------------------------------')
    print('Run target module for miRNA target prediction')
    target_paramter = data['target']

    _q = target_paramter['query_fa']
    _b = target_paramter['subject_fa']
    _o = target_paramter['output']
    # _M = target_paramter['total_misp']
    # if _M == None:
    #     _M = 'off'
    # _m = target_paramter['seed_misp']
    # if _m == None:
    #     _m = 'off'
    # _f = target_paramter['score']
    # if _f == None:
    #     _f = 4
    # _I = target_paramter['mimics']
    # if _I == None:
    #     _I = 'off'
    # _i = target_paramter['mimics_str']
    # if _i == None:
    #     _i = 0
    _T = target_paramter['threads']
    if _T == None:
        _T = 1

    # cmd = f"phasiHunter target -q {_q} -b {_b} -o {_o} -M {_M} -m {_m} -f {_f} -I {_I} -i {_i} -T {_T}" 
    cmd = f"phasiHunter target -q {_q} -b {_b} -o {_o} -T {_T}" 

    print(cmd)
    os.system(cmd)

def Initiator(data):
    print('-----------------------------------------------------')
    print('Run initiator module')
    print('Adjust miRNA_target file format')
    initiator_parameter = data['initiator']

    _i = initiator_parameter['i_input_file']
    _j = initiator_parameter['j_input_file']
    _ip = initiator_parameter['p_input_file']
    _pd = initiator_parameter['sRNA_distance']
    _pl = initiator_parameter['phase_length']
    _ps = initiator_parameter['cleavage_shift']
    _o = initiator_parameter['outputfile']

    cmd = f"phasiHunter initiator -i {_i} -j {_j} -ip {_ip} -pd {_pd} -pl {_pl} -ps {_ps} -o {_o}" 

    print(cmd)
    os.system(cmd)

def Deg(data):
    print('-----------------------------------------------------')
    print('Run deg module for phase initiator verification')
    deg_parameter = data['deg']

    _i = deg_parameter['inputfile']
    _q = deg_parameter['query_fa']
    _j = deg_parameter['STI_result']
    _t = deg_parameter['transcript_fa']
    _o = deg_parameter['output']
    _s = deg_parameter['shift']
    _m = deg_parameter['minum_deg_abun']
    _p = deg_parameter['T_plot']
    _in = deg_parameter['initiator']
    _pl = deg_parameter['plot_categories']
    _pf = deg_parameter['plot_folder']
    __lib = deg_parameter['library']
    _less = deg_parameter['less']

    if _less == 'y':
        _less = '-less'
    else:
        _less = ''

    preprocess_subprocess = []
    if len(__lib) != len(_i):
        print('There are not enough library name, ignore add library column')
        for i in range(0, len(_i)):
            cmd = f"phasiHunter deg -i {_i[i]} -q {_q} -j {_j} -t {_t} -o {_o[i]} -s {_s} -m {_m} -p {_p} -in {_in} -pl {_pl} -pf {_pf} {_less}"
            print(cmd)
            process = subprocess.Popen(cmd, shell=True)
            preprocess_subprocess.append(process)
    else:
        for i in range(0, len(_i)):
            cmd = f"phasiHunter deg -i {_i[i]} -q {_q} -j {_j} -t {_t} -o {_o[i]} -s {_s} -m {_m} -p {_p} -in {_in} -pl {_pl} -pf {_pf} --lib {__lib[i]} {_less}"
            print(cmd)
            process = subprocess.Popen(cmd, shell=True)
            preprocess_subprocess.append(process)

    for process in preprocess_subprocess:
        process.wait()

    all_finished = all([process.poll() is not None for process in preprocess_subprocess])

    if all_finished:
        print("All deg module analysis finished")

def phasiRNA_target(data):
    print('=====================================================')
    print('PhasiRNA target prediction and verification')
    print('-----------------------------------------------------')
    print('Run target module for phasiRNA target prediction')
    phasiRNA_target_parameter = data['phasiRNA_target']

    _q = phasiRNA_target_parameter['query_fa']
    _b = phasiRNA_target_parameter['subject_fa']
    _o = phasiRNA_target_parameter['output']
    # _M = phasiRNA_target_parameter['total_misp']
    # if _M == None:
    #     _M = 'off'
    # _m = phasiRNA_target_parameter['seed_misp']
    # if _m == None:
    #     _m = 'off'
    # _f = phasiRNA_target_parameter['score']
    # if _f == None:
    #     _f = 4
    # _I = phasiRNA_target_parameter['mimics']
    # if _I == None:
    #     _I = 'off'
    # _i = phasiRNA_target_parameter['mimics_str']
    # if _i == None:
    #     _i = 0
    _T = phasiRNA_target_parameter['threads']
    if _T == None:
        _T = 1

    # cmd = f"phasiHunter target -q {_q} -b {_b} -o {_o} -M {_M} -m {_m} -f {_f} -I {_I} -i {_i} -T {_T}" 
    cmd = f"phasiHunter target -q {_q} -b {_b} -o {_o} -T {_T}" 

    print(cmd)
    os.system(cmd)

def phasiRNA_deg(data):
    print('-----------------------------------------------------')
    print('Run deg module for phasiRNA target verification')

    phasiRNA_deg_parameter = data['phasiRNA_deg']

    _i = phasiRNA_deg_parameter['inputfile']
    _q = phasiRNA_deg_parameter['query_fa']
    _j = phasiRNA_deg_parameter['STI_result']
    _t = phasiRNA_deg_parameter['transcript_fa']
    _o = phasiRNA_deg_parameter['output']
    _s = phasiRNA_deg_parameter['shift']
    _m = phasiRNA_deg_parameter['minum_deg_abun']
    _p = phasiRNA_deg_parameter['T_plot']
    _in = phasiRNA_deg_parameter['initiator']
    _pl = phasiRNA_deg_parameter['plot_categories']
    _pf = phasiRNA_deg_parameter['plot_folder']
    __lib = phasiRNA_deg_parameter['library']
    _less = phasiRNA_deg_parameter['less']

    if _less == 'y':
        _less = '-less'
    else:
        _less = ''

    if _in == 'y':
        print('WARNIG: -in must no when executing phasiRNA target verification')
        sys.exit()

    preprocess_subprocess = []
    if len(__lib) != len(_i):
        print('There are not enough library name, ignore add library column')
        for i in range(0, len(_i)):
            cmd = f"phasiHunter deg -i {_i[i]} -q {_q} -j {_j} -t {_t} -o {_o[i]} -s {_s} -m {_m} -p {_p} -in {_in} -pl {_pl} -pf {_pf} {_less}"
            print(cmd)
            process = subprocess.Popen(cmd, shell=True)
            preprocess_subprocess.append(process)
    else:
        for i in range(0, len(_i)):
            cmd = f"phasiHunter deg -i {_i[i]} -q {_q} -j {_j} -t {_t} -o {_o[i]} -s {_s} -m {_m} -p {_p} -in {_in} -pl {_pl} -pf {_pf} --lib {__lib[i]} {_less}"
            print(cmd)
            process = subprocess.Popen(cmd, shell=True)
            preprocess_subprocess.append(process)

    for process in preprocess_subprocess:
        process.wait()

    all_finished = all([process.poll() is not None for process in preprocess_subprocess])

    if all_finished:
        print("All deg module analysis finished")

if __name__ == '__main__':
    main()