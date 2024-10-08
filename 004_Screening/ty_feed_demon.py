import sys, os, time
import subprocess, glob, json
from Daemon import Daemon
import openbabel as ob
# import pickle
from DFT_function import *
# from check import *
class job(Daemon):
    def run(self):
        while True:
            max_jobs = 60
            
            ### oxo
            dft_recode = '/home/u7046082/auto_submit/metal-oxo_already_submit.csv'
            with open('/home/u7046082/auto_submit/metal-oxo_need_submit.csv','r') as f1:
                files = f1.read().splitlines()            
            with open(dft_recode,'r') as f:
                submited_dft = f.read().splitlines()
            unsubmit_list = list(set(files)-set(submited_dft))
            unsubmit_list_sorted = sorted(unsubmit_list)
            next_submit = unsubmit_list_sorted[0]
            cmd = 'squeue -u u7046082'
            cmds = cmd.split(' ')
            process = subprocess.Popen(cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = process.communicate()
            num_jobs = len(out.decode('utf-8').split('\n')) - 2
            if num_jobs < max_jobs:
                DFT_sp('/home/u7046082/VHTS/scripts/manager_utils/manager-final.py','/home/u7046082/xyz_complex_oxo/' + next_submit)
                with open(dft_recode,'a') as f5:
                    f5.write(next_submit+'\n')    
            time.sleep(180)
            DFT_penta_need_submit('/work/u7046082/dft_data_management/done','/home/u7046082/auto_submit/meatl-penta_need_submit.csv')
            time.sleep(200)
            
if __name__ == '__main__':
    daemon = job('/home/u7046082/auto_submit/metal-oxo.pid')
    if len(sys.argv) == 2:
        if 'start' == sys.argv[1]:
            print('daemon start')
            daemon.start()
        elif 'stop' == sys.argv[1]:
            daemon.stop()
        elif 'restart' == sys.argv[1]:
            daemon.restart()
        else:
            print("Unknown command")
            sys.exit(2)

def DFT_sp(manager, xyz, outf='/work/u7046082/DFT_complex_oxo/', details=None):
    inputs = '/home/u7046082/VHTS/scripts/orca/orca5_tpss_d4_lanl2dz_631gpol_opt.inp'    
    charge = xyz.split('.')[0].split('_')[-4]
    multi = xyz.split('.')[0].split('_')[-1]
    metal = xyz.split('.')[0].split('_')[-2]    
    cmd = 'python %s -i %s -x %s -o %s -c %s -m %s -d %s' % (manager, inputs, xyz, outf, charge, multi, metal)
    cmds = cmd.split(' ')
    p = subprocess.Popen(cmds, stdin=None, stdout=None, stderr=None, close_fds=True)
    time.sleep(10)









