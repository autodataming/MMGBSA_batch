#!/usr/bin/env python3


# Function: batch MMGBSA calculations of protein (chain A) and peptide (chain B)

## input: A1.cif A2.cif A3.cif


from glob import glob
import gemmi
import os
import subprocess
import shutil
import re
from pathlib import Path


def extract_residues(topology_file, output_file, offset=0):
    """使用grep和awk提取残基信息"""
    cmd = f'grep "; residue" {topology_file} | awk \'{{print ${3+offset},$4}}\''

    cmd = f'grep "; residue" {topology_file} | awk \'{{print ${3}+{offset},$4}}\''


    try:
        result = subprocess.run(cmd, shell=True, check=True,
                            capture_output=True, text=True)
        with open(output_file, 'w') as f:
            f.write(result.stdout)
        print(f"已提取残基信息到 {output_file}")
        return len(result.stdout.splitlines())
    except subprocess.CalledProcessError as e:
        print(f"错误: 提取残基失败 - {e.stderr}")
        return 0

def create_index(md_file, chain_a_file, chain_b_file, output_file):
    """创建GROMACS索引文件"""
    # 读取残基编号
    with open(chain_a_file) as f:
        res_a = [line.split()[0] for line in f]
    with open(chain_b_file) as f:
        res_b = [line.split()[0] for line in f]

    # 生成make_ndx输入
    input_script = f"""
    ri {' '.join(res_a)}
    name 17 Chain_A
    ri {' '.join(res_b)}
    name 18 Chain_B
    q
    """
    # 写入临时文件
    temp_input = "make_ndx_input.txt"
    with open(temp_input, 'w') as f:
        f.write(input_script)

    # 运行make_ndx
    cmd = f"gmx make_ndx -f {md_file} -o {output_file} < {temp_input}"
    try:
        subprocess.run(cmd, shell=True, check=True)
        print(f"已创建索引文件 {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"错误: 创建索引失败 - {e.stderr}")
    finally:
        if os.path.exists(temp_input):
            os.remove(temp_input)






def extract_gbsa_energy(file_path):
    """
    从 FINAL_RESULTS_MMPBSA.dat 文件中提取 ΔTOTAL 能量值及其标准差（SD）

    参数:
        file_path (str): 文件路径

    返回:
        tuple: (ΔTOTAL 值, SD) 或 (None, None) 如果未找到
    """
    try:
        with open(file_path, 'r') as file:
            content = file.read()

            # 使用正则表达式匹配 ΔTOTAL 行
            pattern = r"ΔTOTAL\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)"
            match = re.search(pattern, content)

            if match:
                delta_total = float(match.group(1))  # ΔTOTAL 平均值
                sd_prop = float(match.group(2))      # SD(Prop.)
                sd = float(match.group(3))           # SD
                sem_prop = float(match.group(4))     # SEM(Prop.)
                sem = float(match.group(5))          # SEM

                # 返回 ΔTOTAL 和 SD（标准偏差）
                return delta_total, sd
            else:
                return None, None

    except FileNotFoundError:
        print(f"错误: 文件 {file_path} 未找到")
        return None, None
    except Exception as e:
        print(f"处理文件时出错: {e}")
        return None, None



cifs = glob("*cif")

for cif in cifs:


    # 读取 CIF 文件
    cif_obj = gemmi.read_structure(cif)
    outf = cif.replace('.cif','.pdb')

    # 保存为 PDB 文件
    cif_obj.write_pdb(outf)



pdbfs = glob("*pdb")


fw = open("GBSA_results.txt",'w')

for pdbf in pdbfs:
    pdb_name = os.path.basename(pdbf)
    dir_name = os.path.splitext(os.path.basename(pdbf))[0]
    base_name = os.path.splitext(os.path.basename(pdbf))[0]
    # 创建文件夹（如果不存在）
    os.makedirs(dir_name, exist_ok=True)  # Python 3.2+ 推荐
    print(f"文件夹 '{dir_name}' 已创建或已存在")

    # 复制文件到文件夹
    dst_path = os.path.join(dir_name, os.path.basename(pdbf))
    shutil.copy(pdbf, dst_path)
    print(f"文件 '{pdbf}' 已复制到 '{dst_path}'")


    # 切换到新文件夹（可选）
    os.chdir(dir_name)
    print(f"当前工作目录已切换到: {os.getcwd()}")

    # 在这里执行需要在新目录下运行的代码（如 gmx 命令）
    # 例如：os.system("gmx pdb2gmx -f A0.pdb -o processed.gro")




    # pdb to gro

    command = "gmx pdb2gmx -f %s.pdb -o %s.gro -ignh -ter -ff amber03 -water tip3p"%(base_name,base_name)

    # 执行命令（捕获输出和错误）
    result = subprocess.run(command, shell=True, check=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    print("命令执行成功！")
    print("输出内容：", result.stdout)
    if result.stderr:
        print("错误信息：", result.stderr)








    # cp 9TEM.gro complex.gro
    grof = "%s.gro"%base_name
    newf="complex.gro"
    shutil.copy(grof, newf)
    print(f"文件 '{grof}' 已复制到 '{newf}'")



    print("pdb2gmx 已成功完成！")








    # add box
    command = "gmx editconf -f complex.gro -o newbox.gro -c -bt cubic -d 1.0"


    # 执行命令（捕获输出和错误）
    result = subprocess.run(command, shell=True, check=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    print("命令执行成功！")
    print("输出内容：", result.stdout)
    if result.stderr:
        print("错误信息：", result.stderr)

    print(" 盒子已成功加完！")




    # add water
    command = "gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro"

    # 执行命令（捕获输出和错误）
    result = subprocess.run(command, shell=True, check=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    print("命令执行成功！")
    print("输出内容：", result.stdout)
    if result.stderr:
        print("错误信息：", result.stderr)

    print(" 水已成功加完！")





    # 定义要写入的内容
    mdp_content = """
    ; LINES STARTING WITH ';' ARE COMMENTS
    title               = Minimization      ; Title of run

    ; Parameters describing what to do, when to stop and what to save
    integrator          = steep             ; Algorithm (steep = steepest descent minimization)
    emtol               = 1000.0    ; Stop minimization when the maximum force < 10.0 kJ/mol
    emstep          = 0.01      ; Energy step size
    nsteps              = 50000             ; Maximum number of (minimization) steps to perform

    ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
    nstlist             = 1             ; Frequency to update the neighbor list and long range forces
    cutoff-scheme   = Verlet
    ns_type             = grid              ; Method to determine neighbor list (simple, grid)
    rlist               = 1.0               ; Cut-off for making neighbor list (short range forces)
    coulombtype         = cutoff    ; Treatment of long range electrostatic interactions
    rcoulomb            = 1.0               ; long range electrostatic cut-off
    rvdw                = 1.0               ; long range Van der Waals cut-off
    pbc             = xyz           ; Periodic Boundary Conditions
    """

    # 写入文件
    with open('ions.mdp', 'w') as f:
        f.write(mdp_content.strip())  # strip() 移除开头和结尾的多余空行

    print("ions.mdp 文件已成功生成！")


    # add ions
    command = "gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr"

    # 执行命令（捕获输出和错误）
    result = subprocess.run(command, shell=True, check=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    print("命令执行成功！")
    print("输出内容：", result.stdout)
    if result.stderr:
        print("错误信息：", result.stderr)


    command = 'echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top       -pname NA -nname CL -neutral -conc 0.15'
    # 执行命令（捕获输出和错误）
    result = subprocess.run(command, shell=True, check=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    print("命令执行成功！")
    print("输出内容：", result.stdout)
    if result.stderr:
        print("错误信息：", result.stderr)



    print(" 离子已成功加完！")


    # 能量最小化

    # 定义要写入的内容
    mdp_content = """
    ; LINES STARTING WITH ';' ARE COMMENTS
    title               = Minimization      ; Title of run

    ; Parameters describing what to do, when to stop and what to save
    integrator          = steep             ; Algorithm (steep = steepest descent minimization)
    emtol               = 1000.0    ; Stop minimization when the maximum force < 10.0 kJ/mol
    emstep          = 0.01      ; Energy step size
    nsteps              = 50000             ; Maximum number of (minimization) steps to perform

    ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
    nstlist             = 1                 ; Frequency to update the neighbor list and long range forces
    cutoff-scheme   = Verlet
    ns_type             = grid                  ; Method to determine neighbor list (simple, grid)
    rlist               = 1.2                   ; Cut-off for making neighbor list (short range forces)
    coulombtype         = PME                   ; Treatment of long range electrostatic interactions
    rcoulomb            = 1.2                   ; long range electrostatic cut-off
    vdwtype         = cutoff
    vdw-modifier    = force-switch
    rvdw-switch     = 1.0
    rvdw                = 1.2                   ; long range Van der Waals cut-off
    pbc             = xyz               ; Periodic Boundary Conditions
    DispCorr        = no
    """

    # 写入文件
    with open('em.mdp', 'w') as f:
        f.write(mdp_content.strip())  # strip() 移除开头和结尾的多余空行

    print("em.mdp 文件已成功生成！")



    command = "gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr"

   # 执行命令（捕获输出和错误）
    result = subprocess.run(command, shell=True, check=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    print("命令执行成功！")
    print("输出内容：", result.stdout)
    if result.stderr:
        print("错误信息：", result.stderr)


    command = "gmx mdrun -v -deffnm em -nt 4 "

   # 执行命令（捕获输出和错误）
    result = subprocess.run(command, shell=True, check=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    print("命令执行成功！")
    print("输出内容：", result.stdout)
    if result.stderr:
        print("错误信息：", result.stderr)

    print("em 已经完成！")

    # nvt 500 ps

    nsteps = 250000
    tcgrps = "Protein Water_and_ions"
    temperature = 310

    mdp_content = f"""title                   = Protein-ligand complex NVT equilibration
    define                  = -DPOSRES  ; position restrain the protein and ligand
    ; Run parameters
    integrator              = md        ; leap-frog integrator
    nsteps                  = {nsteps}    ; 2 * 50000 = 100 ps
    dt                      = 0.002     ; 2 fs
    ; Output control
    nstenergy               = 500   ; save energies every 1.0 ps
    nstlog                  = 500   ; update log file every 1.0 ps
    nstxout-compressed      = 500   ; save coordinates every 1.0 ps
    ; Bond parameters
    continuation            = no        ; first dynamics run
    constraint_algorithm    = lincs     ; holonomic constraints
    constraints             = h-bonds   ; bonds to H are constrained
    lincs_iter              = 1         ; accuracy of LINCS
    lincs_order             = 4         ; also related to accuracy
    ; Neighbor searching and vdW
    cutoff-scheme           = Verlet
    ns_type                 = grid      ; search neighboring grid cells
    nstlist                 = 20        ; largely irrelevant with Verlet
    rlist                   = 1.2
    vdwtype                 = cutoff
    vdw-modifier            = force-switch
    rvdw-switch             = 1.0
    rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
    ; Electrostatics
    coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
    rcoulomb                = 1.2       ; short-range electrostatic cutoff (in nm)
    pme_order               = 4         ; cubic interpolation
    fourierspacing          = 0.16      ; grid spacing for FFT
    ; Temperature coupling
    tcoupl                  = V-rescale                     ; modified Berendsen thermostat
    tc-grps                 = {tcgrps}    ; two coupling groups - more accurate
    tau_t                   = 0.1   0.1                     ; time constant, in ps
    ref_t                   = {temperature}   {temperature}                     ; reference temperature, one for each group, in K
    ; Pressure coupling
    pcoupl                  = no        ; no pressure coupling in NVT
    ; Periodic boundary conditions
    pbc                     = xyz       ; 3-D PBC
    ; Dispersion correction is not used for proteins with the C36 additive FF
    DispCorr                = no
    ; Velocity generation
    gen_vel                 = yes       ; assign velocities from Maxwell distribution
    gen_temp                = {temperature}       ; temperature for Maxwell distribution
    gen_seed                = -1        ; generate a random seed
    """


    # 写入文件
    with open('nvt.mdp', 'w') as f:
        f.write(mdp_content.strip())  # strip() 移除开头和结尾的多余空行

    print("nvt.mdp 文件已成功生成！")

    command = "gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top  -o nvt.tpr;gmx mdrun -deffnm nvt -v -nt 16"



   # 执行命令（捕获输出和错误）
    result = subprocess.run(command, shell=True, check=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    print("命令执行成功！")
    print("输出内容：", result.stdout)
    if result.stderr:
        print("错误信息：", result.stderr)

    print("nvt  已经完成！")


    # npt


    # npt 500 ps

    nsteps = 250000
    tcgrps = "Protein Water_and_ions"

    mdp_content = f"""title                   = Protein-peptide complex NPT equilibration
    define                  = -DPOSRES  ; position restrain the protein and ligand
    ; Run parameters
    integrator              = md        ; leap-frog integrator
    nsteps                  = {nsteps}     ; 2 * 50000 = 100 ps
    dt                      = 0.002     ; 2 fs
    ; Output control
    nstenergy               = 500       ; save energies every 1.0 ps
    nstlog                  = 500       ; update log file every 1.0 ps
    nstxout-compressed      = 500       ; save coordinates every 1.0 ps
    ; Bond parameters
    continuation            = yes       ; continuing from NVT
    constraint_algorithm    = lincs     ; holonomic constraints
    constraints             = h-bonds   ; bonds to H are constrained
    lincs_iter              = 1         ; accuracy of LINCS
    lincs_order             = 4         ; also related to accuracy
    ; Neighbor searching and vdW
    cutoff-scheme           = Verlet
    ns_type                 = grid      ; search neighboring grid cells
    nstlist                 = 20        ; largely irrelevant with Verlet
    rlist                   = 1.2
    vdwtype                 = cutoff
    vdw-modifier            = force-switch
    rvdw-switch             = 1.0
    rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
    ; Electrostatics
    coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
    rcoulomb                = 1.2
    pme_order               = 4         ; cubic interpolation
    fourierspacing          = 0.16      ; grid spacing for FFT
    ; Temperature coupling
    tcoupl                  = V-rescale                     ; modified Berendsen thermostat
    tc-grps                 = {tcgrps}    ; two coupling groups - more accurate
    tau_t                   = 0.1   0.1                     ; time constant, in ps
    ref_t                   = {temperature}   {temperature}                     ; reference temperature, one for each group, in K
    ; Pressure coupling
    pcoupl                  = C-rescale                      ; pressure coupling is on for NPT
    pcoupltype              = isotropic                     ; uniform scaling of box vectors
    tau_p                   = 2.0                           ; time constant, in ps
    ref_p                   = 1.0                           ; reference pressure, in bar
    compressibility         = 4.5e-5                        ; isothermal compressibility of water, bar^-1
    refcoord_scaling        = com
    ; Periodic boundary conditions
    pbc                     = xyz       ; 3-D PBC
    ; Dispersion correction is not used for proteins with the C36 additive FF
    DispCorr                = no
    ; Velocity generation
    gen_vel                 = no        ; velocity generation off after NVT
    """





    # 写入文件
    with open('npt.mdp', 'w') as f:
        f.write(mdp_content.strip())  # strip() 移除开头和结尾的多余空行

    print("npt.mdp 文件已成功生成！")





    command = "gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top  -o npt.tpr;gmx mdrun -deffnm npt -v -nt 16"

   # 执行命令（捕获输出和错误）
    result = subprocess.run(command, shell=True, check=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    print("命令执行成功！")
    print("输出内容：", result.stdout)
    if result.stderr:
        print("错误信息：", result.stderr)

    print("npt  已经完成！")




    # md


    nsteps = 250000
    tcgrps = "Protein Water_and_ions"
    temperature =310
    mdp_content = f"""title                   = Protein-ligand complex MD simulation
    ; Run parameters
    integrator              = md        ; leap-frog integrator
    nsteps                  = {nsteps}   ; 2 * 5000000 = 10000 ps (10 ns)
    dt                      = 0.002     ; 2 fs
    ; Output control
    nstenergy               = 5000      ; save energies every 10.0 ps
    nstlog                  = 5000      ; update log file every 10.0 ps
    nstxout-compressed      = 5000      ; save coordinates every 10.0 ps
    ; Bond parameters
    continuation            = yes       ; continuing from NPT
    constraint_algorithm    = lincs     ; holonomic constraints
    constraints             = h-bonds   ; bonds to H are constrained
    lincs_iter              = 1         ; accuracy of LINCS
    lincs_order             = 4         ; also related to accuracy
    ; Neighbor searching and vdW
    cutoff-scheme           = Verlet
    ns_type                 = grid      ; search neighboring grid cells
    nstlist                 = 20        ; largely irrelevant with Verlet
    rlist                   = 1.2
    vdwtype                 = cutoff
    vdw-modifier            = force-switch
    rvdw-switch             = 1.0
    rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
    ; Electrostatics
    coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
    rcoulomb                = 1.2
    pme_order               = 4         ; cubic interpolation
    fourierspacing          = 0.16      ; grid spacing for FFT
    ; Temperature coupling
    tcoupl                  = V-rescale                     ; modified Berendsen thermostat
    tc-grps                 = {tcgrps}    ; two coupling groups - more accurate
    tau_t                   = 0.1   0.1                     ; time constant, in ps
    ref_t                   = {temperature}  {temperature}                     ; reference temperature, one for each group, in K
    ; Pressure coupling
    pcoupl                  = Parrinello-Rahman             ; pressure coupling is on for NPT
    pcoupltype              = isotropic                     ; uniform scaling of box vectors
    tau_p                   = 2.0                           ; time constant, in ps
    ref_p                   = 1.0                           ; reference pressure, in bar
    compressibility         = 4.5e-5                        ; isothermal compressibility of water, bar^-1
    ; Periodic boundary conditions
    pbc                     = xyz       ; 3-D PBC
    ; Dispersion correction is not used for proteins with the C36 additive FF
    DispCorr                = no
    ; Velocity generation
    gen_vel                 = no        ; continuing from NPT equilibration
    """





    # 写入文件
    with open('md.mdp', 'w') as f:
        f.write(mdp_content.strip())  # strip() 移除开头和结尾的多余空行

    print("md.mdp 文件已成功生成！")







    command = "gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top  -o md_500ps.tpr; gmx mdrun -deffnm md_500ps -v  -nt 16"


   # 执行命令（捕获输出和错误）
    result = subprocess.run(command, shell=True, check=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    print("命令执行成功！")
    print("输出内容：", result.stdout)
    if result.stderr:
        print("错误信息：", result.stderr)

    print("md  已经完成！")


    # 消除轨迹平动 PBC
    command = 'echo "Protein" | gmx trjconv -s md_500ps.tpr -f md_500ps.xtc -o md_500ps_nojump.xtc -pbc nojump'

    result = subprocess.run(command, shell=True, check=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    print("命令执行成功！")
    print("输出内容：", result.stdout)
    if result.stderr:
        print("错误信息：", result.stderr)

    print("PBC 周期消除 已经完成！")




    command='echo -e "Protein\nProtein" | gmx trjconv -f md_500ps_nojump.xtc   -s md_500ps.tpr  -fit rot+trans -o md_500ps_nojump_fit.xtc '


    result = subprocess.run(command, shell=True, check=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    print("命令执行成功！")
    print("输出内容：", result.stdout)
    if result.stderr:
        print("错误信息：", result.stderr)

    print("平动 消除 已经完成！")


    ## MMGBSA
    ## chain A 和 chain B 的编号，PDB 文件和gro list 对应起来


    # 文件路径设置
    current_dir = Path.cwd()
    chain_a_top = current_dir / "topol_Protein_chain_A.itp"
    chain_b_top = current_dir / "topol_Protein_chain_B.itp"
    md_file = current_dir / "md_500ps.gro"
    index_file = current_dir / "index.ndx"

    # 1. 提取Chain A残基
    chain_a_file = current_dir / "chain_A_residues.dat"
    chain_a_count = extract_residues(chain_a_top, chain_a_file)

    if chain_a_count == 0:
        print("错误: 未提取到Chain A残基")
        raise "ERROR"

    # 2. 提取Chain B残基(自动偏移)
    chain_b_file = current_dir / "chain_B_residues.dat"
    extract_residues(chain_b_top, chain_b_file, offset=chain_a_count)

    # 3. 创建索引文件
    create_index(md_file, chain_a_file, chain_b_file, index_file)

    #raise "stop" 

    intdiel = 3
    temperature =310
    mmpbsa_content = f"""
    Input file generated by gmx_MMPBSA (1.6.4)
    Be careful with the variables you modify, some can have severe consequences on the results you obtain.

    # General namelist variables
    &general
    sys_name             = "Protein peptide GBSA calculation"                                             # System name
    startframe           = 1                                              # First frame to analyze
    endframe             = 9999999                                        # Last frame to analyze
    interval             = 1                                              # Number of frames between adjacent frames analyzed
    forcefields          =  "oldff/leaprc.ff03"                           # Define the force field to build the Amber topology
    ions_parameters      = 1                                              # Define ions parameters to build the Amber topology
    PBRadii              = 3                                              # Define PBRadii to build amber topology from GROMACS files
    temperature          = {temperature}                                        # Temperature
    qh_entropy           = 0                                              # Do quasi-harmonic calculation
    interaction_entropy  = 0                                              # Do Interaction Entropy calculation
    ie_segment           = 25                                             # Trajectory segment to calculate interaction entropy
    c2_entropy           = 0                                              # Do C2 Entropy calculation
    assign_chainID       = 0                                              # Assign chains ID
    exp_ki               = 0.0                                            # Experimental Ki in nM
    full_traj            = 0                                              # Print a full traj. AND the thread trajectories
    gmx_path             = ""                                             # Force to use this path to get GROMACS executable
    keep_files           = 2                                              # How many files to keep after successful completion
    netcdf               = 0                                              # Use NetCDF intermediate trajectories
    solvated_trajectory  = 1                                              # Define if it is necessary to cleanup the trajectories
    verbose              = 1                                              # How many energy terms to print in the final output
    /

    # (AMBER) Generalized-Born namelist variables
    &gb
    igb                  = 5                                              # GB model to use
    intdiel              = {intdiel}                                           # Internal dielectric constant for sander
    extdiel              = 78.5                                           # External dielectric constant for sander
    saltcon              = 0.0                                            # Salt concentration (M)
    surften              = 0.0072                                         # Surface tension
    surfoff              = 0.0                                            # Surface tension offset
    molsurf              = 0                                              # Use Connelly surface ('molsurf' program)
    msoffset             = 0.0                                            # Offset for molsurf calculation
    probe                = 1.4                                            # Solvent probe radius for surface area calc
    ifqnt                = 0                                              # Use QM on part of the system
    qm_theory            = ""                                             # Semi-empirical QM theory to use
    qm_residues          = ""                                             # Residues to treat with QM
    com_qmmask           = ""                                             # Mask specifying the quantum atoms in complex
    rec_qmmask           = ""                                             # Mask specifying the quantum atoms in receptor
    lig_qmmask           = ""                                             # Mask specifying the quantum atoms in ligand
    qmcharge_com         = 0                                              # Charge of QM region in complex
    qmcharge_lig         = 0                                              # Charge of QM region in ligand
    qmcharge_rec         = 0                                              # Charge of QM region in receptor
    qmcut                = 9999.0                                         # Cutoff in the QM region
    scfconv              = 1e-08                                          # Convergence criteria for the SCF calculation, in kcal/mol
    peptide_corr         = 0                                              # Apply MM correction to peptide linkages
    writepdb             = 1                                              # Write a PDB file of the selected QM region
    verbosity            = 0                                              # Controls the verbosity of QM/MM related output
    alpb                 = 0                                              # Use Analytical Linearized Poisson-Boltzmann (ALPB)
    arad_method          = 1                                              # Selected method to estimate the effective electrostatic size
    /

    """


    # 写入文件
    with open('mmpbsa.in', 'w') as f:
        f.write(mmpbsa_content.strip())  # strip() 移除开头和结尾的多余空行

    print("mmpbsa 文件已成功生成！")



    command = "mpirun   -np 16     gmx_MMPBSA -O -i mmpbsa.in -cs   md_500ps.tpr  -ct md_500ps_nojump_fit.xtc  -ci index.ndx -cg 17 18 -cp topol.top -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv -nogui  "


   # 执行命令（捕获输出和错误）
    result = subprocess.run(command, shell=True, check=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    print("命令执行成功！")
    print("输出内容：", result.stdout)
    if result.stderr:
        print("错误信息：", result.stderr)

    print("md  已经完成！")





    # 使用示例
    file_path = "FINAL_RESULTS_MMPBSA.dat"
    delta_total, sd = extract_gbsa_energy(file_path)

    if delta_total is not None and sd is not None:
        print(f"ΔTOTAL 能量: {delta_total:.2f} ± {sd:.2f} kcal/mol")
    else:
        print("未找到 ΔTOTAL 能量信息")

    outline="%s %s %s\n"%(base_name,delta_total,sd)
    fw.write(outline)

    # 切换回上级目录（继续处理下一个文件）
    os.chdir("..")



fw.close()













