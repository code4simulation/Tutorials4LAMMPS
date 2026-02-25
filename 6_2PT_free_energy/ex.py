import os
import subprocess
import numpy as np
from scipy.fft import fft, fftfreq
from scipy.optimize import brentq
from scipy import constants  # Add this import
import matplotlib.pyplot as plt

# =====================================================================
# Physical Constants (using scipy.constants)
# =====================================================================
kB = constants.k  # Boltzmann constant (J/K)
h_planck = constants.h  # Planck constant (J·s)
NA = constants.N_A  # Avogadro's number (1/mol)
eV_to_J = constants.e  # Elementary charge (converts eV to Joules)

# =====================================================================
# 1. NPT 평형화 및 부피(Volume) 이동 평균 추출
# =====================================================================
def run_npt_and_get_volume(phase, temp, npt_steps=50000):
    if phase == "solid":
        init_cmd = f"velocity all create {temp} 12345 dist gaussian"
    else:
        init_cmd = """
        velocity all create 3000 12345 dist gaussian
        fix melt all nvt temp 3000 3000 0.1
        run 20000
        unfix melt
        """
        
    lammps_npt = f"""
    units metal
    boundary p p p
    atom_style atomic
    lattice diamond 5.431
    region box block 0 4 0 4 0 4
    create_box 1 box
    create_atoms 1 box
    mass 1 28.0855
    pair_style sw
    pair_coeff * * Si.sw Si

    {init_cmd}
    
    # NPT 실행 및 부피 로깅
    fix 1 all npt temp {temp} {temp} 0.1 iso 0.0 0.0 1.0
    variable V equal vol
    fix vol_log all ave/time 1 100 100 v_V file vol_npt_{phase}_{temp}.txt
    
    run {npt_steps}
    write_data data.npt_{phase}_{temp}
    """
    
    script_name = f"in.npt_{phase}_{temp}"
    with open(script_name, "w") as f:
        f.write(lammps_npt)
    subprocess.run(["lmp", "-in", script_name], capture_output=True)
    
    # Python에서 이동 평균(Moving Average) 계산 (후반부 50% 데이터 사용)
    vol_data = np.loadtxt(f"vol_npt_{phase}_{temp}.txt", skiprows=2)[:, 1]
    half_idx = len(vol_data) // 2
    V_avg = np.mean(vol_data[half_idx:]) # 수렴 구간의 평균 부피
    print(f"[{phase.upper()} @ {temp}K] Converged Volume: {V_avg:.4f} A^3")
    return V_avg

# =====================================================================
# 2. VACF 수렴 검사기 (Convergence Checker)
# =====================================================================
def check_vacf_convergence(cv_file, mean_threshold=0.005, std_diff_threshold=0.005):
    """
    고체상의 영구적인 진동 노이즈를 고려하여 VACF의 수렴성을 평가합니다.
    - mean_threshold: 0으로의 거시적 수렴성 (Drift 확인)
    - std_diff_threshold: 감쇄(Decay)가 멈추고 노이즈 플로어에 도달했는지 확인
    """
    data = np.loadtxt(cv_file, skiprows=2)
    cv = data[:, 1]
    cv_norm = cv / cv[0]
    
    # 마지막 40% 데이터를 두 개의 블록(20%씩)으로 분할
    total_len = len(cv_norm)
    block1 = cv_norm[-int(total_len * 0.4) : -int(total_len * 0.2)] # 끝에서 40% ~ 20% 구간
    block2 = cv_norm[-int(total_len * 0.2) :]                       # 끝에서 20% ~ 0% 구간
    
    # 마지막 구간의 평균 (0 근처인지 확인)
    mean_tail = np.abs(np.mean(block2))
    
    # 두 구간의 진폭(Std) 변화량 비교
    std1 = np.std(block1)
    std2 = np.std(block2)
    std_diff = np.abs(std1 - std2)
    
    # 평균이 충분히 0에 가깝고, 진폭이 더 이상 감소하지 않으면 수렴으로 판정
    is_converged = (mean_tail < mean_threshold) and (std_diff < std_diff_threshold)
    
    return is_converged, mean_tail, std2, std_diff

# =====================================================================
# 3. NVT 실행 및 동적 시간 연장 (Auto-extension)
# =====================================================================
def run_nvt_vacf_loop(phase, temp, V_avg, init_steps=50000, max_steps=200000):
    # 평균 부피(V_avg)에 맞춰 Box 길이(L) 계산 (정육면체 가정)
    L = V_avg**(1/3)
    current_steps = init_steps
    
    while current_steps <= max_steps:
        print(f"  -> Running NVT for {current_steps} steps ({current_steps/1000} ps)...")
        lammps_nvt = f"""
        units metal
        boundary p p p
        atom_style atomic
        
        read_data data.npt_{phase}_{temp}
        pair_style sw
        pair_coeff * * Si.sw Si
        
        # 앞서 구한 V_avg에 맞게 박스 크기 리스케일링
        change_box all x final 0 {L} y final 0 {L} z final 0 {L} remap units box
        
        # NVT 및 VACF 추출
        fix 1 all nvt temp {temp} {temp} 0.1
        compute myVACF all vacf
        fix vacf_log all ave/time 1 1 1 c_myVACF[1] file vacf_{phase}_{temp}.txt
        
        variable E equal pe
        fix energy_log all ave/time 1 1000 1000 v_E file energy_{phase}_{temp}.txt
        
        run {current_steps}
        """
        
        script_name = f"in.nvt_{phase}_{temp}"
        with open(script_name, "w") as f:
            f.write(lammps_nvt)
        subprocess.run(["lmp", "-in", script_name], capture_output=True)
        
        # 수렴 여부 확인
        vacf_file = f"vacf_{phase}_{temp}.txt"
        converged, mean_val, std_val, std_diff = check_vacf_convergence(vacf_file)

        if converged:
            print(f"  -> [PASS] VACF converged. (Mean: {mean_val:.4f}, Std Diff: {std_diff:.4f})")
            break
        else:
            print(f"  -> [FAIL] Not converged (Mean: {mean_val:.4f}, Std Diff: {std_diff:.4f}). Extending...")
            current_steps += 50000 # 50ps 연장
            
    if current_steps > max_steps:
        print("  -> [WARNING] Reached maximum steps without perfect convergence. Proceeding anyway.")

# =====================================================================
# 4. 2PT 자유 에너지 계산 (NumPy 1.x / 2.x 호환 및 표기법 업데이트)
# =====================================================================
def calculate_free_energy_with_v(cv_file, energy_file, T, V_avg_A3, N=512, m_avg=28.0855, dt=0.001):
    """
    속도 자기상관 함수(VACF, Cv)와 NPT 평균 부피를 기반으로 2PT 자유 에너지를 계산합니다.
    
    Parameters:
    - cv_file: VACF correlation file
    - energy_file: Potential energy file
    - T: Temperature (K)
    - V_avg_A3: Average volume (Ångströms³)
    - N: Number of atoms
    - m_avg: Average atomic mass (g/mol)
    - dt: Time step (picoseconds)
    """
    # [호환성 처리] NumPy 2.0 이상이면 trapezoid, 미만이면 trapz 사용
    if hasattr(np, 'trapezoid'):
        integrate_func = np.trapezoid
    else:
        integrate_func = np.trapz

    # 1. Cv(t) 데이터 로드 및 Mirroring 처리
    data = np.loadtxt(cv_file, skiprows=2)
    cv_norm = data[:, 1] / data[0, 1]
    cv_mirrored = np.concatenate((cv_norm[::-1], cv_norm[1:]))
    
    # 2. FFT를 통한 진동 상태 밀도 g(w) 추출
    n_total = len(cv_mirrored)
    dos_fft = fft(cv_mirrored) * dt
    freqs = fftfreq(n_total, d=dt)
    pos_mask = freqs > 0
    f_pos, g_omega = freqs[pos_mask], 2 * np.real(dos_fft[pos_mask])
    
    # g(w) 정규화 (호환성 적분 함수 적용)
    g_omega *= (3 * N) / integrate_func(g_omega, f_pos)
    
    # 3. 유체화 인자 (f) 계산을 위한 파라미터 설정
    g0 = g_omega[0]
    m_kg = m_avg / (NA * 1000)  # Convert g/mol to kg
    V_m3 = V_avg_A3 * 1e-30  # Convert Ångströms³ to m³
    rho = N / V_m3  # Density (atoms/m³)
    
    # 무차원 확산 파라미터 (Delta)
    delta = (g0 / (12 * N)) * np.sqrt(np.pi * kB * T / m_kg) * (rho**(1/3)) * (6/np.pi)**(1/3)
    
    # 방정식: 2 * Delta^-1 * f^2.5 + f^0.5 - 1 = 0 풀이
    try:
        f_sol = brentq(lambda f: 2 * (f**2.5) / delta + f**0.5 - 1, 1e-10, 1.0)
    except ValueError:
        f_sol = 0.0 # 해가 없거나 수렴하지 않으면 고체(f=0)로 간주
        
    # 4. g(w) 분리 (Solid vs Gas)
    if f_sol > 0:
        g_gas = g0 / (1 + (np.pi * g0 * f_pos / (6 * max(f_sol, 1e-10) * N))**2)
    else:
        g_gas = np.zeros_like(g_omega)
        
    g_solid = np.clip(g_omega - g_gas, 0, None)
    
    # 5. 엔트로피(S) 및 역학적 가중치(W) 계산
    beta = 1.0 / (kB * T)
    x = h_planck * (f_pos * 1e12) * beta # 주파수를 Hz 단위로 변환 후 적용
    
    # 고체 성분 가중치 (Quantum Harmonic Oscillator)
    W_S_solid = x / (np.exp(x) - 1) - np.log(1 - np.exp(-x))
    
    # 엔트로피 적분 (호환성 적분 함수 적용)
    Entropy_solid = kB * integrate_func(g_solid * W_S_solid, f_pos)
    Entropy_gas = kB * integrate_func(g_gas, f_pos) * f_sol # 약식 기체 가중치 적용
    Total_Entropy = Entropy_solid + Entropy_gas
    
    # 6. 내부 에너지(E) 및 헬름홀츠 자유 에너지(A) 계산
    pe_data = np.loadtxt(energy_file, skiprows=2)
    E_pot = np.mean(pe_data[:, 1]) * eV_to_J  # Convert eV to Joules
    E_kin = 1.5 * N * kB * T
    Total_Energy = E_pot + E_kin
    
    A = Total_Energy - T * Total_Entropy
    return (A / eV_to_J) / N # eV/atom 단위로 반환

# =====================================================================
# 5. 메인 실행 루프
# =====================================================================
temperatures = [1500, 1600, 1700, 1800] # 시간 단축을 위해 간격을 넓힘 (필요시 조정)
G_solid, G_liquid = [], []

for T in temperatures:
    print(f"\n--- [Solid Phase @ {T}K] ---")
    V_avg_s = run_npt_and_get_volume("solid", T)
    run_nvt_vacf_loop("solid", T, V_avg_s)
    Gs = calculate_free_energy_with_v(f"vacf_solid_{T}.txt", f"energy_solid_{T}.txt", T, V_avg_s)
    G_solid.append(Gs)
    
    print(f"\n--- [Liquid Phase @ {T}K] ---")
    V_avg_l = run_npt_and_get_volume("liquid", T)
    run_nvt_vacf_loop("liquid", T, V_avg_l)
    Gl = calculate_free_energy_with_v(f"vacf_liquid_{T}.txt", f"energy_liquid_{T}.txt", T, V_avg_l)
    G_liquid.append(Gl)

print("\n계산 완료! 자유 에너지(eV/atom) 결과를 확인하세요.")

# =====================================================================
# 6. 녹는점(Tm) 계산 및 자유 에너지 시각화
# =====================================================================

# 데이터 배열화
temps = np.array(temperatures)
G_s = np.array(G_solid)
G_l = np.array(G_liquid)

# 온도가 1개뿐이면 피팅을 할 수 없으므로 예외 처리
if len(temps) < 2:
    print("\n[알림] 녹는점을 계산하려면 최소 2개 이상의 온도 데이터가 필요합니다.")
    print(f"현재 1500K에서의 ΔG (Liquid - Solid): {G_l[0] - G_s[0]:.4f} eV/atom")

else:
    # 1. ΔG 계산 (Liquid - Solid)
    delta_G = G_l - G_s
    
    # 2. 다항식 피팅 (온도 구간이 좁으면 1차 선형 피팅, 넓으면 2차 피팅 권장)
    # 데이터가 3개 이상이면 2차(Quadratic), 2개면 1차(Linear) 피팅 적용
    fit_order = 2 if len(temps) >= 3 else 1
    poly_coeffs = np.polyfit(temps, delta_G, fit_order)
    poly_func = np.poly1d(poly_coeffs)
    
    # 3. 방정식 해(Roots) 찾기 (ΔG = 0이 되는 온도 Tm)
    roots = poly_func.roots
    # 물리적으로 타당한 범위(예: 1000K ~ 2500K) 내의 실수 해만 선택
    valid_roots = [r.real for r in roots if np.isreal(r) and 1000 < r.real < 2500]
    
    if not valid_roots:
        print("\n[경고] 지정된 온도 범위 내에서 교차점(녹는점)을 찾지 못했습니다.")
        Tm = None
    else:
        Tm = valid_roots[0] # 가장 가까운 해 선택
        print(f"\n✅ 계산된 예측 녹는점 (Tm): {Tm:.2f} K")

    # 4. 시각화 (Plotting)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # [Subplot 1] 절대 자유 에너지 (G vs T)
    ax1.plot(temps, G_s, 'bo-', label='Solid (a-Si / c-Si)')
    ax1.plot(temps, G_l, 'ro-', label='Liquid')
    if Tm:
        ax1.axvline(Tm, color='k', linestyle='--', alpha=0.5, label=f'$T_m$ = {Tm:.1f} K')
    ax1.set_xlabel('Temperature (K)')
    ax1.set_ylabel('Free Energy (eV/atom)')
    ax1.set_title('Absolute Free Energy')
    ax1.legend()
    ax1.grid(True, linestyle=':', alpha=0.7)

    # [Subplot 2] 자유 에너지 차이 (ΔG vs T)
    ax2.plot(temps, delta_G, 'ko', label='MD Data')
    
    # 피팅 커브 생성을 위한 촘촘한 온도 배열
    T_dense = np.linspace(temps[0] - 50, temps[-1] + 50, 100)
    ax2.plot(T_dense, poly_func(T_dense), 'b-', label=f'{fit_order}nd Order Fit')
    
    ax2.axhline(0, color='r', linestyle='-', linewidth=1)
    if Tm:
        ax2.axvline(Tm, color='k', linestyle='--', alpha=0.5)
        ax2.plot(Tm, 0, 'r*', markersize=12, label='Melting Point')
        
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel('$\Delta G$ (Liquid - Solid) [eV/atom]')
    ax2.set_title('Free Energy Difference')
    ax2.legend()
    ax2.grid(True, linestyle=':', alpha=0.7)

    plt.tight_layout()
    plt.show()
