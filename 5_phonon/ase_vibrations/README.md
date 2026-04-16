# Theoretical Foundations: Partial Hessian Vibrational Analysis (PHVA) & Statistical Thermodynamics

이 문서는 표면 화학 및 전산 재료 과학에서 분자와 고체 표면 간의 상호작용을 열역학적으로 모델링하기 위해 작성된 `analyze_phva.py` 및 `calc_free_energy.py` 모듈의 배경 이론을 고체물리 및 통계역학 교과서 수준으로 해설합니다.

---

## 1. 진동 분석 (Vibrational Analysis)의 수학적 정식화

### 1.1 퍼텐셜 에너지 표면 (Potential Energy Surface, PES) 곡률
Born-Oppenheimer 근사 하에서 원자핵의 거동은 전자의 바닥 상태 에너지가 만들어내는 퍼텐셜 에너지 표면(PES), $V(\mathbf{R})$ 위에서 움직이는 고전적/양자적 입자로 취급됩니다. 

시스템이 국소적 최소점(Local Minimum) $\mathbf{R}_0$ 에 도달했을 때, 에너지를 미세 변위 $\Delta\mathbf{R}$ 에 대해 테일러 전개(Taylor Expansion)하면 다음과 같습니다.

$$ V(\mathbf{R}) = V(\mathbf{R}_0) + \sum_{i=1}^{3N} \left( \frac{\partial V}{\partial R_i} \right)_0 \Delta R_i + \frac{1}{2} \sum_{i=1}^{3N} \sum_{j=1}^{3N} \left( \frac{\partial^2 V}{\partial R_i \partial R_j} \right)_0 \Delta R_i \Delta R_j + \cdots $$

최소점에서는 1차 미분(힘, Force)이 0이므로, 에너지를 결정하는 핵심 항은 **Hessian 행렬(H)**이라 불리는 2차 미분항이 됩니다.
$$ H_{ij} = \frac{\partial^2 V}{\partial R_i \partial R_j} $$

### 1.2 질량-가중 헤시안 (Mass-Weighted Hessian)
원자마다 질량이 다르므로, Newton의 운동 방정식($F = ma$)을 풀기 위해 질량으로 규격화한 질량-가중 헤시안 행렬 $\tilde{H}$를 정의합니다.

$$ \tilde{H}_{ij} = \frac{1}{\sqrt{m_i m_j}} H_{ij} $$

이 행렬의 고유값(Eigenvalue) $\lambda_k$ 를 구하면 그것이 곧 고유 진동수의 제곱이 됩니다.
$$ \tilde{H} \mathbf{v}_k = \lambda_k \mathbf{v}_k \quad \implies \quad \lambda_k = \omega_k^2 = (2\pi \nu_k)^2 $$

---

## 2. Partial Hessian Vibrational Analysis (PHVA)의 근사 원리

총 원자 수가 $N$개인 슬랩(Slab) 구조 시스템에서 전체 자유도는 $3N$입니다. 시스템이 크더라도 표면 반응(Chemidesorption 등)에 참여하는 원소는 국소적입니다.

### 2.1 블록 행렬화 및 무한 질량 근사
헤시안 행렬을 진동에 참여하는 **활성 원자 그룹(Active, $A$)**과 완전히 멈춰있다고 가정하는 **동결 원자 그룹(Frozen, $F$)**으로 분할합니다.

동결 원자는 수학적으로 **무한대의 질량($m_F \to \infty$)**을 가진다고 가정하는 것과 동치입니다. 따라서 질량-가중 헤시안의 동결 원자 관여 항($\tilde{H}_{AF}, \tilde{H}_{FA}, \tilde{H}_{FF}$)은 모두 0으로 수렴합니다.

$$ \tilde{H} = \begin{pmatrix} \tilde{H}_{AA} & 0 \\ 0 & 0 \end{pmatrix} $$

결과적으로 $3N \times 3N$ 행렬을 물리적 손실을 최소화하면서 $3N_A \times 3N_A$ 행렬($N_A$: 활성 원자 수)로 축소할 수 있으며, 이로써 연산량(Force calculation 횟수)을 극적으로 절감합니다. 구현된 파이썬 코드에서 `ase.vibrations.Vibrations` 객체의 `indices` 인자가 바로 활성 부분공간($A$)을 정의하는 역할을 합니다.

---

## 3. 전이 상태 표상과 허수 진동수 (Imaginary Frequencies)

물리적으로 안정된 바닥 상태(Ground State)는 모든 방향에 대해 곡률이 양수(Convex)이므로, 모든 $\lambda_k > 0$ 이고 따라서 $\nu_k$ 는 실수(Real)로 주어집니다. 

하지만 반응 경로(Reaction coordinate)를 따라 형성되는 1차 안장점(1st-Order Saddle Point), 즉 **전이 상태(Transition State, TS)**에서는 딱 하나의 방향에 대해 위치 에너지가 극대값을 가집니다.

$$ \frac{\partial^2 V}{\partial R_{TS}^2} < 0 \quad \implies \quad \lambda_{TS} < 0 $$

$\omega^2$ 이 음수이므로, 진동수 $\nu$ 는 복소수의 허수(Imaginary, $i$) 형태를 띱니다. 
> [!WARNING]
> 허수 진동수는 그 방향성이 "복원력(Restoring Force)이 아닌 파괴력(Accelerating Force)"으로 작용함을 뜻합니다. `calc_free_energy.py` 모듈에서 설계한 `clean_frequencies` 로직은 수치 속 허수의 개수(`n_imag`)를 파악하여, 이 구조가 안정 상태인지(0개), 전이 상태인지(1개), 아니면 물리적으로 붕괴하는 상태인지(>1개) 판단하는 잣대로 기능합니다.

---

## 4. 통계역학과 열역학적 자유 에너지 모델

구조의 전자 에너지($E_{DFT}$)와 진동 주파수를 모두 구했다면, 거시적 물성인 깁스 자유 에너지($G$)를 통계역학적 분배 함수(Partition Function, $Z$)로부터 도출합니다.
$$ G(T, P) = -k_B T \ln Z_{tot} + pV $$

계산 코드는 대상의 상(Phase)에 따라 두 가지 철저히 구분된 열역학 모델을 차용합니다.

### 4.1 고체 기판 및 흡착물 (Harmonic Approximation)
기판을 이루는 고체나 단단히 흡착된 종(Adsorbate)은 병진 운동(Trans)과 자유 회전(Rot)이 억제된 채, 특정 격자점에 구속되어 진동(Frustrated Vibration)만 수행합니다. 따라서 분배 함수는 진동항 $Z_{vib}$ 만을 포함합니다.

1. **영점 에너지 (Zero-Point Energy, ZPE)**: 양자역학적 불확정성 원리에 의해 절대 영도($0$ K)에서도 존재하는 기저 진동 에너지.
   $$ E_{ZPE} = \sum_i \frac{1}{2} h \nu_i $$
2. **진동 내부에너지 ($U_{vib}$)** 및 **진동 엔트로피 ($S_{vib}$)**:
   $$ U_{vib}(T) = \sum_i \frac{h\nu_i}{e^{h\nu_i / k_BT} - 1} $$
   $$ S_{vib}(T) = \sum_i k_B \left[ \frac{h\nu_i / k_BT}{e^{h\nu_i / k_BT} - 1} - \ln(1 - e^{-h\nu_i / k_BT}) \right] $$

고체는 부피 변화로 인한 일($pV$)이 무시할 수 있을 만큼 작습니다. 즉, 압력(P)에 상관 없이 다음과 같은 **고체 조화 진동자 방정식**이 도출되며, 코드는 이를 따릅니다.
$$ G^{solid}(T) \approx F(T) = E_{DFT} + E_{ZPE} + U_{vib}(T) - T \cdot S_{vib}(T) $$

### 4.2 기체 분자 (Ideal Gas Approximation)
기체 분자는 공간을 자유롭게 날아다니고 회전하므로 모든 분배 함수가 발현됩니다 ($Z_{tot} = Z_{trans} Z_{rot} Z_{vib} Z_{elec}$).
코드는 `ase.thermochemistry.IdealGasThermo` 규약을 엄격히 따릅니다.

1. **병진 엔트로피 (Sackur-Tetrode Equation)**: 여기서 시스템의 **분압(Partial Pressure, $P$)** 의존성이 강하게 발생합니다. 압력이 커질수록 부피($V$)가 수축하여 가용 미시상태(Microstates)가 줄어들기 때문입니다.
   $$ S_{trans} = n R \left[ \ln \left( \frac{V}{n} \left( \frac{2\pi m k_B T}{h^2} \right)^{3/2} \right) + \frac{5}{2} \right] \quad \left(V = \frac{nRT}{P} \right) $$
2. **회전 엔트로피 ($S_{rot}$)**: 분자의 관성 모멘트 기하($Linear$ vs $Nonlinear$)와 **대칭수(Symmetry Number, $\sigma$)**에 의존합니다.

기체의 깁스 자유 에너지는 온도와 더불어 사용자 입력 **압력($P$)**에 따라 즉각적으로 보정됩니다.
$$ G^{gas}(T, P) = E_{DFT} + E_{ZPE} + \Delta H_{trans+rot+vib}(T) - T \cdot \left[ S_{trans}(T,P) + S_{rot}(T) + S_{vib}(T) \right] $$

---

## 5. 결론 및 코드 연동 정리

1. [analyze_phva.py](file:///c:/Users/user/.gemini/antigravity/playground/exo-photosphere/analyze_phva.py) 스크립트는 Born-Oppenheimer 및 무한 질량 근사를 통해 $3N_A \times 3N_A$ 의 질량-가중 헤시안 행렬을 풀고 $(\nu_i)$를 획득합니다.
2. [calc_free_energy.py](file:///c:/Users/user/.gemini/antigravity/playground/exo-photosphere/calc_free_energy.py)의 `clean_frequencies` 함수는 $\nu_i$의 부호를 검사하여 음수/허수를 분류하고, $\nu_i$를 반응 좌표(Reaction coordinate)로 쓸지, 제외할지, 에러를 던질지 판별합니다.
3. 남은 온전한 $\nu_i$ 스펙트럼은 상(Phase)의 물리적 속성에 따라 **체적 의존성**이 있는 기체(Ideal Gas) 공식과, 체적 의존성이 탈각된 고체 조화 진동(Harmonic) 공식으로 흘러가 최종 열역학적 거시량 $G$를 출력합니다.
