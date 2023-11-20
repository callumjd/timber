# timber

import os
import glob
import numpy as np
from .ligprep_tools import check_file

##############################################################################

class OnlineAvVar(object):
    def __init__(self, store_data = True):
        self.step = 0
        self.mean = 0.0
        self.M2 = 0.0
        self.store = store_data
        self.data = []

    def accumulate(self, x):
        self.step += 1
        delta = x - self.mean
        self.mean += delta / self.step
        self.M2 += delta * (x - self.mean)
        if self.store:
            self.data.append(x)

    def get_variance(self):
        return self.M2 / (self.step - 1)

    def get_stat(self):
        return self.mean, math.sqrt(self.M2 / (self.step - 1))

def run_analysis(df):

    # analysis all dir lig1->lig2
    for index,row in df.iterrows():
        name1_col=list(df.columns)[0]
        name1=row[name1_col]
        # rbfe
        if len(list(df.columns))==2:
            name2_col=list(df.columns)[1]
            name2=row[name2_col]
            pair_dir=name1+'~'+name2
        # absolute
        else:
            pair_dir=name1

    os.chdir(pair_dir)

    # get protocol and ti_repeats
    ti_repeats=len(glob.glob('complex/*_rep*'))
    protocol=glob.glob('complex/*_rep*')[0].split('_rep')[0].replace('complex/','')

    # run the analysis
    if criteria_analysis():
        output_dG=perform_analysis(protocol,ti_repeats)
        # print dG results
        if len(output_dG)>1:
            output_dG=np.array(output_dG)
            print('%s calc dG kcal/mol: %lf +/- %lf\n' % (pair_dir,np.average(output_dG,axis=0),np.std(output_dG,axis=0)))
        elif len(output_dG)==1:
            print('%s calc dG kcal/mol: %lf\n' % (pair_dir,output_dG[0]))
        else:
            print('Error: calc dG failed for %s\n' % (pair_dir))

    else:
        print('Error: cannot run analysis %s\n' % (pair_dir))
        print('Checking for failed prod runs\n')
        resubmit_gpu_error()

    # pair dir
    os.chdir('../')

def perform_analysis(protocol,ti_repeats):

    output_dG=[]

    # analytical solvent restraint for absolute runs
    if protocol in ['absolute','absolute-three-step']:
        solv_rest=get_analytical_rest('disang.RST')
    else:
        solv_rest=0.0

    complex_stage_dict={}
    complex_lambda_dict={}
    complex_free_enegy={}

    solvent_stage_dict={}
    solvent_lambda_dict={}
    solvent_free_enegy={}

    for rep in range(1,ti_repeats+1):

        # complex
        os.chdir('complex/%s_rep%d' % (protocol,rep))

        for stage_dir in glob.glob('*/'):
            stage=stage_dir.strip('/')
            os.chdir(stage_dir)

            dvdl_list=get_dvdl('lambda_*/')
            complex_stage_dict.update({stage:dvdl_list})

            lambda_list=get_lambda('lambda_*/')
            complex_lambda_dict.update({stage:lambda_list})

            complex_free_enegy.update({stage:FE_DVDL(complex_stage_dict,complex_lambda_dict)})

            write_dvdl_data('DVDL.dat',dvdl_list,lambda_list)

            os.chdir('../')

        # complex dir
        os.chdir('../../')

        # solvent
        os.chdir('solvent/%s_rep%d' % (protocol,rep))

        for stage_dir in glob.glob('*/'):
            stage=stage_dir.strip('/')
            os.chdir(stage_dir)

            dvdl_list=get_dvdl('lambda_*/')
            solvent_stage_dict.update({stage:dvdl_list})

            lambda_list=get_lambda('lambda_*/')
            solvent_lambda_dict.update({stage:lambda_list})

            solvent_free_enegy.update({stage:FE_DVDL(solvent_stage_dict,solvent_lambda_dict)})

            write_dvdl_data('DVDL.dat',dvdl_list,lambda_list)

            os.chdir('../')

        # solvent dir
        os.chdir('../../')

        # check if we can run convergence 
        convergence=True
        dG_array=np.zeros(len(complex_free_enegy['complex_ligands']))

        for stage,fe_list in complex_free_enegy.items():
            if len(fe_list)!=len(dG_array):
                convergence=False
                dG_array=np.zeros(1)
                break

        for stage,fe_list in solvent_free_enegy.items():
            if len(fe_list)!=len(dG_array):
                convergence=False
                dG_array=np.zeros(1)
                break

        # final free energy array
        for stage,fe_list in complex_free_enegy.items():
            if convergence:
                dG_array+=np.array(fe_list)
            else:
                dG_array[0]+=fe_list[-1]

        # TEST signs for absolute
        for stage,fe_list in solvent_free_enegy.items():
            if convergence:
                dG_array-=(np.array(fe_list)-solv_rest)
            else:
                dG_array[0]-=(fe_list[-1]-solv_rest)

        output_dG.append(list(dG_array)[-1])

        # write convergence and contributions
        write_ti_convergence(protocol,rep,dG_array,complex_free_enegy,solvent_free_enegy,solv_rest)

    return output_dG

def write_dvdl_data(dvdl_filename,dvdl_list,lambda_list):

    with open(dvdl_filename,'w') as f:
        for i in range(0,len(lambda_list)):
            f.write('%lf %lf\n' % (lambda_list[i],np.mean(dvdl_list[i])))

def write_ti_convergence(protocol,rep,dG_array,complex_free_enegy,solvent_free_enegy,solv_rest):

    with open('convergence_%s_rep%d.dat' % (protocol,rep),'w') as f:
        for i in range(0,len(dG_array)):
            f.write('%d %lf\n' % (i,dG_array[i]))

    with open('breakdown_%s_rep%d.dat' % (protocol,rep),'w') as f:
        for stage in sorted(list(complex_free_enegy.keys())):
            f.write('%s %lf\n' % (stage,complex_free_enegy[stage][-1]))

        for stage in sorted(list(solvent_free_enegy.keys())):
            f.write('%s %lf\n' % (stage,solvent_free_enegy[stage][-1]))
        if abs(solv_rest)>0:
            f.write('solvent_restraint %lf\n' % (solv_rest))

def FE_DVDL(stage_dict,lambda_dict):

    fe_stage_dict={}
    for stage,dvdl_list in stage_dict.items():
        lambda_values=lambda_dict[stage]

        weights=None

        # gaussian quadrature weights
        if lambda_values==[0.1127,0.5,0.88729]:
            weights=np.array([0.27777,0.44444,0.27777])
        elif lambda_values==[0.04691,0.23076,0.5,0.76923,0.95308]:
            weights=np.array([0.11846,0.23931,0.28444,0.23931,0.11846])
        elif lambda_values==[0.02544,0.12923,0.29707,0.5,0.70292,0.87076,0.97455]:
            weights=np.array([0.06474,0.13985,0.19091,0.20897,0.19091,0.13985,0.06474])
        elif lambda_values==[0.01592,0.08198,0.19331,0.33787,0.5,0.66213,0.80669,0.91802,0.98408]:
            weights=np.array([0.04064, 0.09032, 0.13031, 0.15617, 0.16512, 0.15617, 0.13031, 0.09032, 0.04064])
        elif lambda_values==[0.00922,0.04794,0.11505,0.20634,0.31608,0.43738,0.56262,0.68392,0.79366,0.88495,0.95206,0.99078]:
            weights=np.array([0.02359,0.05347,0.08004,0.10158,0.11675,0.12457,0.12457,0.11675,0.10158,0.08004,0.05347,0.02359])

        dG_list=[]

        # if all windows have the same number of dvdl values, calculate dG convergence
        conv_test=0
        dvdl_len=len(dvdl_list[0])
        for dvdl_windows in dvdl_list:
            if len(dvdl_windows)==dvdl_len:
                conv_test+=1

        if conv_test==len(dvdl_list):
            convergence=True
        else:
            convergence=False

        if convergence:

            for conv in range(1,dvdl_len+1):
                dvdl_frac=[]
                for i in range(0,len(lambda_values)):
                    dvdl_frac.append(np.mean(dvdl_list[i][0:conv]))

                if isinstance(weights,np.ndarray):
                    run_sum=0
                    for gs in range(0,len(lambda_values)):
                        run_sum+=(dvdl_frac[gs]*weights[gs])
                else:
                    run_sum=np.trapz(dvdl_frac,x=lambda_values)

                dG_list.append(run_sum)

        else:

            dvdl_frac=[]
            for i in range(0,len(lambda_values)):
                dvdl_frac.append(np.mean(dvdl_list[i]))

            if isinstance(weights,np.ndarray):
                run_sum=0
                for gs in range(0,len(lambda_values)):
                    run_sum+=(dvdl_frac[gs]*weights[gs])
            else:
                run_sum=np.trapz(dvdl_frac,x=lambda_values)

            dG_list.append(run_sum)

    return dG_list

def get_dvdl(dir_str):
    
    # list of dvdl lists for each lambda window
    dvdl_list=[]

    file_list=glob.glob(dir_str)
    file_list=sorted(file_list)

    for lambda_dir in file_list:
        if check_file(lambda_dir+'prod.en'):
            dVdl = OnlineAvVar()

            with open(lambda_dir+'prod.en', 'r') as en_file:
                for line in en_file:
                    if line.startswith('L9') and not 'dV/dlambda' in line:
                        dVdl.accumulate(float(line.split()[5]))

            # check for extended prod runs
            if len(glob.glob(lambda_dir+'extend*.en'))>0:
                for extend in sorted(glob.glob(lambda_dir+'extend*.en')):
                    if check_file(extend):

                        with open(extend, 'r') as en_file:
                            for line in en_file:
                                if line.startswith('L9') and not 'dV/dlambda' in line:
                                    dVdl.accumulate(float(line.split()[5]))
                        
        dvdl_list.append(dVdl.data)

    return dvdl_list

def get_lambda(dir_str):

    # list of lambda values
    lambda_list=[]

    file_list=glob.glob(dir_str)
    file_list=sorted(file_list)

    for window in file_list:
        if check_file(window+'prod.in'):
            with open(window+'prod.in','r') as f:
                for line in f:
                    if 'clambda =' in line:
                        val=float(line.split('clambda =')[1].split(',')[0])
                        lambda_list.append(val)

    return lambda_list

def criteria_analysis():

    output=[]

    dir_list=glob.glob('complex/*_rep*/*/lambda_*/')+glob.glob('solvent/*_rep*/*/lambda_*/')

    for mydir in dir_list:
        output.append(check_file(mydir+'prod.en'))

    return all(output)

def resubmit_gpu_error():

    # TESTING REQUIRED
    dir_list=glob.glob('complex/*_rep*/*/')+glob.glob('solvent/*_rep*/*/')

    cwd=os.getcwd()
    for mydir in dir_list:
        os.chdir(mydir)

        runfiles=glob.glob('TI.o*')
        runfiles=sorted(runfiles, key=lambda x: int(x.split('.')[-1]),reverse=False)
        to_resubmit=[]
        for run in runfiles:
            with open(run,'r') as f:
                data=f.readlines()
            if len(data)>0 and 'Error' in data[-1]:
                run_val=int(run.split('.')[-1])-1
                to_resubmit.append(run_val)
            elif len(data)>0 and 'STOP PMEMD Terminated Abnormally!' in data[-1]:
                run_val=int(run.split('.')[-1])-1
                to_resubmit.append(run_val)
            elif len(data)>0 and 'Segmentation fault' in data[-1]:
                run_val=int(run.split('.')[-1])-1
                to_resubmit.append(run_val)
            elif len(data)>0 and 'Error: an illegal memory access' in data[-1]:
                run_val=int(run.split('.')[-1])-1
                to_resubmit.append(run_val)
            elif len(data)>0 and 'cudaMemcpy GpuBuffer::Download failed' in data[-1]:
                run_val=int(run.split('.')[-1])-1
                to_resubmit.append(run_val)

        if len(to_resubmit)>0:
            print('Re-submitting %d runs\n' % (len(to_resubmit)))
            if check_file('dir_list.dat'):
                with open('dir_list.dat','r') as f:
                    dat_list=f.readlines()
                os.system('rm dir_list.dat')
                if len(to_resubmit)<=len(dat_list):
                    with open('dir_list.dat','w') as f_out:
                        for val in to_resubmit:
                            f_out.write('%s' % (dat_list[val]))

                    os.system('rm TI.o*')
                    #os.system('qsub -t 1-%d run_prod.sh' % (len(to_resubmit)))
                else:
                    print('Error with resubmission!\n')

        os.chdir(cwd)

def get_analytical_rest(disang_file,temp=300):

    K=1.9872 * 1e-3   # Gas constant in kcal/mol/K
    V=1.66            # standard volume in nm^3
    T=float(temp)     # Temperature in Kelvin

    with open(disang_file,'r') as f:
        f.readline()
        data=f.readlines()

    for val in data[0].split(','):
        if 'r2' in val:
            r0=0.1*float(val.strip().strip('r2=').strip())
        elif 'rk2' in val:
            K_r=float(val.strip().strip('rk2=').strip())

    for val in data[1].split(','):
        if 'r2' in val:
            thA=float(val.strip().strip('r2=').strip())
        elif 'rk2' in val:
            K_thA=float(val.strip().strip('rk2=').strip())

    for val in data[2].split(','):
        if 'r2' in val:
            thB=float(val.strip().strip('r2=').strip())
        elif 'rk2' in val:
            K_thB=float(val.strip().strip('rk2=').strip())

    for val in data[3].split(','):
        if 'rk2' in val:
            K_phiA=float(val.strip().strip('rk2=').strip())

    for val in data[4].split(','):
        if 'rk2' in val:
            K_phiB=float(val.strip().strip('rk2=').strip())

    for val in data[5].split(','):
        if 'rk2' in val:
            K_phiC=float(val.strip().strip('rk2=').strip())

    # BORESCH FORMULA

    thA = math.radians(thA)
    thB = math.radians(thB)

    arg =(
    (8.0 * math.pi**2.0 * V) / (r0**2.0 * math.sin(thA) * math.sin(thB))
    *
    (
        ( (K_r * K_thA * K_thB * K_phiA * K_phiB * K_phiC)**0.5 ) / ( (2.0 * math.pi * K * T)**(3.0) )
        )
    )

    dG = - K * T * math.log(arg)

    return -dG

