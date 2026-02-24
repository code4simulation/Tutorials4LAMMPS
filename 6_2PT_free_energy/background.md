# 2PT(Two-Phase Thermodynamics) 기반 자유 에너지 계산 이론 및 시뮬레이션 워크플로우

## 1. 2PT 모델의 이론적 배경

일반적인 결정질 고체의 열역학적 특성은 Phonon에 기반한 조화 진동자(Harmonic Oscillator) 모델을 통해 정확하게 계산할 수 있습니다. 하지만 액체나 비정질 상태의 경우, 원자들이 고정된 위치에서 진동만 하는 것이 아니라 시스템 내를 확산하기 때문에 단순 고체 물리 이론만으로는 엔트로피를 과소평가하게 됩니다.

2PT 모델은 이러한 시스템의 진동 상태 밀도(Vibrational Density of States, VDOS)를 고체적 성분(Solid-like)과 기체적 성분(Gas-like)의 두 가지 상으로 분리하여 정확한 열역학적 물성치를 도출하는 방법론입니다.

### 1.1 속도 자기상관 함수(VACF)와 진동 상태 밀도 ($g(\omega)$)

원자들의 움직임은 분자 동역학(MD) 시뮬레이션에서 속도 자기상관 함수(Velocity Autocorrelation Function, $C_v(t)$)를 통해 추적할 수 있습니다.


$$C_v(t) = \frac{\langle \sum_{i=1}^{N} \mathbf{v}_i(t) \cdot \mathbf{v}_i(0) \rangle}{\langle \sum_{i=1}^{N} \mathbf{v}_i(0) \cdot \mathbf{v}_i(0) \rangle}$$

이 $C_v(t)$를 푸리에 변환(Fourier Transform)하여 주파수 영역으로 넘기면, 시스템의 포논 상태 밀도를 포함하는 전체 진동 상태 밀도 $g(\omega)$를 얻을 수 있습니다. (여기서 $\omega$는 각진동수, Angular frequency 입니다.)


$$g(\omega) = \frac{1}{2\pi} \int_{-\infty}^{\infty} C_v(t) e^{-i\omega t} dt$$

*참고: $C_v(t)$는 짝함수(Even function)이므로, 실제 수치 해석에서는 코사인 변환인 $\int_{0}^{\infty} C_v(t) \cos(\omega t) dt$ 형태로 단순화되어 계산되기도 합니다.*

### 1.2 상태 밀도의 2상 분리 (Two-Phase Partitioning)

2PT 모델은 전체 상태 밀도 $g(\omega)$를 확산 모드와 진동 모드로 나눕니다.


$$g(\omega) = g_{solid}(\omega) + g_{gas}(\omega)$$

* **Gas-like Phase ($g_{gas}$):** 원자의 확산(Diffusion)에 기여하는 성분입니다. 주파수가 0일 때의 값 $g(0)$은 자가 확산 계수(Self-diffusion coefficient)와 직접적으로 연관됩니다.
* **Solid-like Phase ($g_{solid}$):** 원자의 국부적인 진동(Phonon)에 기여하는 성분입니다. 전체 $g(\omega)$에서 기체 성분을 뺀 나머지로 정의됩니다.

### 1.3 유체화 인자 (Fluidicity Factor, $f$)의 도출

$g(\omega)$를 분리하기 위해서는 전체 자유도 중 기체적 성분이 차지하는 비율을 나타내는 유체화 인자 $f$ (where $0 \le f \le 1$)를 구해야 합니다.

2PT 모델은 기체적 성분이 **강체구 유체(Hard-sphere fluid)**의 거동을 따른다고 가정합니다. Hard-sphere 모델의 상태 방정식(Enskog theory)과 분자 동역학적 확산 관계식을 연립하면, 시스템의 무차원 확산 파라미터 $\Delta$를 다음과 같이 정의할 수 있습니다.


$$\Delta = \frac{g(0)}{12N} \left( \frac{\pi k_B T}{m} \right)^{1/2} \rho^{1/3} \left( \frac{6}{\pi} \right)^{1/3}$$


(여기서 $m$은 원자 질량, $\rho$는 수 밀도 $N/V$ 입니다.)

이 $\Delta$ 값을 바탕으로, 강체구 모델의 충돌 빈도와 관련된 보편적인 비선형 방정식을 풀어 유체화 인자 $f$를 도출합니다.


$$2\Delta^{-1} f^{5/2} + f^{1/2} - 1 = 0$$


이 방정식의 해 $f$를 구하면, 이를 이용해 $g_{gas}(\omega)$의 형태를 특정 짓고 나머지 $g_{solid}(\omega)$를 확정하게 됩니다.

### 1.4 역학적 가중치 함수 (Weighting Functions, $W$)

분리된 상태 밀도를 열역학적 에너지와 엔트로피로 변환하기 위해, 각 주파수 $\omega$에 해당하는 양자 역학 및 통계 역학적 가중치 함수 $W(\omega)$를 곱하여 적분합니다.

**1) 고체 성분 가중치 ($W_{solid}$)**
고체 성분은 **양자 조화 진동자(Quantum Harmonic Oscillator)** 모델을 따릅니다. 따라서 엔트로피를 구하기 위한 가중치 $W_S^{solid}(\omega)$는 포논 통계 역학의 정확한 해를 사용합니다.


$$W_S^{solid}(\omega) = \frac{\hbar \omega / k_B T}{\exp(\hbar \omega / k_B T) - 1} - \ln(1 - \exp(-\hbar \omega / k_B T))$$


이 가중치를 $g_{solid}(\omega)$에 곱해 적분하면 양자 효과가 반영된 고체의 진동 엔트로피가 산출됩니다.

**2) 기체 성분 가중치 ($W_{gas}$)**
기체 성분은 **강체구 기체(Hard-sphere gas)** 모델을 따릅니다. 기체 성분의 엔트로피 가중치 $W_S^{gas}$는 이상 기체의 병진 엔트로피(Sackur-Tetrode equation)에 입자들의 부피에 의한 패킹 효과(Packing fraction)를 보정한 Carnahan-Starling 상태 방정식을 결합하여 유도됩니다. (수식이 매우 복잡하여 Python 코드 내부에서는 해석적 해의 근사치 테이블이나 모듈화된 함수를 통해 계산됩니다.)

최종적으로 총 엔트로피($S$)와 내부 에너지($E$)를 구한 뒤, 헬름홀츠 자유 에너지($A$)를 산출합니다.


$$S_{total} = k_B \int_{0}^{\infty} g_{solid}(\omega) W_S^{solid}(\omega) d\omega + k_B \int_{0}^{\infty} g_{gas}(\omega) W_S^{gas}(\omega) d\omega$$

$$A = E - T S_{total}$$
