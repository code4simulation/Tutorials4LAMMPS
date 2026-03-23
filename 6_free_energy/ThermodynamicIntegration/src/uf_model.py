import numpy as np
import scipy.constants as const
from scipy.interpolate import CubicSpline
from scipy.integrate import trapezoid
import subprocess
import re
import os

class UhlenbeckFord:
    """
    Unified Uhlenbeck-Ford (UF) Model Handler.
    Supports both legacy Calphy splines and the extended Master Table (0-10.0+).
    """
    KB = const.physical_constants["Boltzmann constant in eV/K"][0]
    H = const.physical_constants["Planck constant in eV/Hz"][0]
    NA = const.physical_constants["Avogadro constant"][0]

    # Universal Master Table Data (x in [0, 10.0])
    _X_MASTER = [
        0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 
        0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 
        0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 
        0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 
        0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 
        1.0, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 
        1.2, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.3, 1.31, 1.32, 1.33, 1.34, 1.35, 1.36, 1.37, 1.38, 1.39, 
        1.4, 1.41, 1.42, 1.43, 1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.5, 1.51, 1.52, 1.53, 1.54, 1.55, 1.56, 1.57, 1.58, 1.59, 
        1.6, 1.61, 1.62, 1.63, 1.64, 1.65, 1.66, 1.67, 1.68, 1.69, 1.7, 1.71, 1.72, 1.73, 1.74, 1.75, 1.76, 1.77, 1.78, 1.79, 
        1.8, 1.81, 1.82, 1.83, 1.84, 1.85, 1.86, 1.87, 1.88, 1.89, 1.9, 1.91, 1.92, 1.93, 1.94, 1.95, 1.96, 1.97, 1.98, 1.99, 
        2.0, 2.01, 2.02, 2.03, 2.04, 2.05, 2.06, 2.07, 2.08, 2.09, 2.1, 2.11, 2.12, 2.13, 2.14, 2.15, 2.16, 2.17, 2.18, 2.19, 
        2.2, 2.21, 2.22, 2.23, 2.24, 2.25, 2.26, 2.27, 2.28, 2.29, 2.3, 2.31, 2.32, 2.33, 2.34, 2.35, 2.36, 2.37, 2.38, 2.39, 
        2.4, 2.41, 2.42, 2.43, 2.44, 2.45, 2.46, 2.47, 2.48, 2.49, 2.5, 2.51, 2.52, 2.53, 2.54, 2.55, 2.56, 2.57, 2.58, 2.59, 
        2.6, 2.61, 2.62, 2.63, 2.64, 2.65, 2.66, 2.67, 2.68, 2.69, 2.7, 2.71, 2.72, 2.73, 2.74, 2.75, 2.76, 2.77, 2.78, 2.79, 
        2.8, 2.81, 2.82, 2.83, 2.84, 2.85, 2.86, 2.87, 2.88, 2.89, 2.9, 2.91, 2.92, 2.93, 2.94, 2.95, 2.96, 2.97, 2.98, 2.99, 
        3.0, 3.01, 3.02, 3.03, 3.04, 3.05, 3.06, 3.07, 3.08, 3.09, 3.1, 3.11, 3.12, 3.13, 3.14, 3.15, 3.16, 3.17, 3.18, 3.19, 
        3.2, 3.21, 3.22, 3.23, 3.24, 3.25, 3.26, 3.27, 3.28, 3.29, 3.3, 3.31, 3.32, 3.33, 3.34, 3.35, 3.36, 3.37, 3.38, 3.39, 
        3.4, 3.41, 3.42, 3.43, 3.44, 3.45, 3.46, 3.47, 3.48, 3.49, 3.5, 3.51, 3.52, 3.53, 3.54, 3.55, 3.56, 3.57, 3.58, 3.59, 
        3.6, 3.61, 3.62, 3.63, 3.64, 3.65, 3.66, 3.67, 3.68, 3.69, 3.7, 3.71, 3.72, 3.73, 3.74, 3.75, 3.76, 3.77, 3.78, 3.79, 
        3.8, 3.81, 3.82, 3.83, 3.84, 3.85, 3.86, 3.87, 3.88, 3.89, 3.9, 3.91, 3.92, 3.93, 3.94, 3.95, 3.96, 3.97, 3.98, 3.99, 
        4.0, 5.0, 7.5, 10.0
    ]
    _FEX_MASTER = [
        0.0, 0.0099878, 0.0199818, 0.0300256, 0.0401135, 0.0502185, 0.0603514, 0.0705192, 0.0807039, 0.0909044, 
        0.1011256, 0.1113718, 0.1216435, 0.1319369, 0.1422539, 0.1525951, 0.1629586, 0.1733428, 0.1837483, 0.1941762, 
        0.2046265, 0.2150975, 0.2255877, 0.2360966, 0.2466244, 0.2571715, 0.267738, 0.2783239, 0.2889294, 0.2995542, 
        0.3101984, 0.3208619, 0.3315445, 0.3422456, 0.3529645, 0.363701, 0.3744557, 0.3852288, 0.3960199, 0.4068283, 
        0.4176538, 0.4284963, 0.4393558, 0.4502317, 0.4611234, 0.4720306, 0.4829538, 0.4938935, 0.5048499, 0.5158229, 
        0.5268117, 0.5378154, 0.5488333, 0.5598657, 0.5709135, 0.5819767, 0.5930547, 0.6041466, 0.6152525, 0.6263728, 
        0.6375072, 0.648655, 0.6598153, 0.6709884, 0.6821751, 0.6933757, 0.7045893, 0.7158152, 0.7270535, 0.7383046, 
        0.7495687, 0.7608457, 0.7721351, 0.7834363, 0.7947488, 0.8060724, 0.8174074, 0.8287542, 0.8401129, 0.8514838, 
        0.8628666, 0.8742608, 0.8856657, 0.8970811, 0.9085071, 0.9199439, 0.9313917, 0.9428501, 0.9543183, 0.9657954, 
        0.9772817, 0.9887787, 1.0002878, 1.0118085, 1.0233395, 1.0348796, 1.0464284, 1.0579862, 1.0695534, 1.0811306, 
        1.092718, 1.1043154, 1.1159226, 1.1275391, 1.1391649, 1.1507995, 1.1624428, 1.1740947, 1.1857549, 1.1974234, 
        1.2091001, 1.2207848, 1.2324775, 1.2441782, 1.2558867, 1.2676032, 1.2793275, 1.2910597, 1.3027997, 1.3145475, 
        1.3263031, 1.3380666, 1.3498379, 1.3616169, 1.3734038, 1.3851984, 1.3970005, 1.4088102, 1.4206274, 1.4324518, 
        1.4442833, 1.4561218, 1.4679672, 1.4798194, 1.4916781, 1.5035435, 1.5154153, 1.5272937, 1.5391785, 1.5510698, 
        1.5629676, 1.574872, 1.5867829, 1.5987005, 1.6106246, 1.6225554, 1.6344926, 1.6464363, 1.6583864, 1.6703428, 
        1.6823054, 1.694274, 1.7062486, 1.718229, 1.7302151, 1.7422069, 1.7542043, 1.7662072, 1.7782156, 1.7902294, 
        1.8022487, 1.8142733, 1.8263034, 1.8383389, 1.8503798, 1.8624261, 1.8744776, 1.8865346, 1.8985968, 1.9106644, 
        1.9227372, 1.9348152, 1.9468985, 1.9589869, 1.9710805, 1.9831792, 1.9952828, 2.0073915, 2.019505, 2.0316234, 
        2.0437466, 2.0558744, 2.0680068, 2.0801437, 2.0922852, 2.1044311, 2.1165814, 2.1287362, 2.1408954, 2.1530591, 
        2.1652272, 2.1773999, 2.1895771, 2.2017589, 2.2139453, 2.2261363, 2.2383318, 2.2505319, 2.2627364, 2.2749454, 
        2.2871587, 2.2993764, 2.3115982, 2.3238242, 2.3360543, 2.3482884, 2.3605265, 2.3727685, 2.3850145, 2.3972643, 
        2.4095179, 2.4217754, 2.4340367, 2.4463018, 2.4585706, 2.4708433, 2.4831196, 2.4953997, 2.5076835, 2.519971, 
        2.5322622, 2.544557, 2.5568555, 2.5691577, 2.5814634, 2.5937728, 2.6060857, 2.6184022, 2.6307223, 2.6430458, 
        2.6553729, 2.6677035, 2.6800375, 2.692375, 2.7047159, 2.7170602, 2.7294078, 2.7417588, 2.7541131, 2.7664706, 
        2.7788314, 2.7911955, 2.8035626, 2.815933, 2.8283064, 2.840683, 2.8530627, 2.8654455, 2.8778313, 2.8902202, 
        2.9026122, 2.9150073, 2.9274055, 2.9398068, 2.9522112, 2.9646186, 2.9770291, 2.9894427, 3.0018593, 3.0142789, 
        3.0267014, 3.0391269, 3.0515552, 3.0639864, 3.0764205, 3.0888573, 3.1012968, 3.1137391, 3.1261841, 3.1386317, 
        3.151082, 3.1635349, 3.1759904, 3.1884484, 3.2009091, 3.2133723, 3.2258381, 3.2383064, 3.2507774, 3.2632509, 
        3.275727, 3.2882058, 3.3006872, 3.3131712, 3.325658, 3.3381473, 3.3506393, 3.3631339, 3.3756311, 3.3881309, 
        3.4006331, 3.4131378, 3.4256449, 3.4381543, 3.4506661, 3.4631802, 3.4756965, 3.4882151, 3.500736, 3.5132591, 
        3.5257845, 3.5383123, 3.5508423, 3.5633747, 3.5759094, 3.5884464, 3.6009858, 3.6135275, 3.6260715, 3.6386179, 
        3.6511664, 3.6637173, 3.6762703, 3.6888256, 3.701383, 3.7139425, 3.7265041, 3.7390679, 3.7516337, 3.7642017, 
        3.7767716, 3.7893437, 3.8019177, 3.8144938, 3.827072, 3.8396522, 3.8522343, 3.8648186, 3.8774048, 3.889993, 
        3.9025833, 3.9151756, 3.92777, 3.9403663, 3.9529647, 3.9655651, 3.9781674, 3.9907718, 4.0033781, 4.0159864, 
        4.0285966, 4.0412087, 4.0538226, 4.0664385, 4.0790561, 4.0916756, 4.1042968, 4.1169199, 4.1295448, 4.1421714, 
        4.1547998, 4.1674301, 4.1800621, 4.1926959, 4.2053315, 4.217969, 4.2306082, 4.2432492, 4.255892, 4.2685366, 
        4.281183, 4.2938312, 4.3064811, 4.3191328, 4.3317863, 4.3444414, 4.3570983, 4.369757, 4.3824173, 4.3950793, 
        4.407743, 4.4204083, 4.4330752, 4.4457438, 4.458414, 4.4710857, 4.4837591, 4.4964341, 4.5091106, 4.5217888, 
        4.5344685, 4.5471498, 4.5598327, 4.5725173, 4.5852034, 4.5978912, 4.6105805, 4.6232715, 4.6359641, 4.6486582, 
        4.6613539, 4.6740512, 4.68675, 4.6994503, 4.7121521, 4.7248555, 4.7375603, 4.7502665, 4.7629742, 4.7756833, 
        4.7883937, 6.0632909, 9.281996, 12.5355347
    ]

    def __init__(self, mode='master'):
        self.mode = mode
        self.update_spline()

    def update_spline(self):
        self._cs_fex = CubicSpline(self._X_MASTER, self._FEX_MASTER)

    @staticmethod
    def get_x(rho, sigma):
        """Dimensionless density x = B2 * rho"""
        return 0.5 * (np.pi * sigma**2)**1.5 * rho

    @staticmethod
    def get_rho(x, sigma):
        """Number density rho = x / B2"""
        return x / (0.5 * (np.pi * sigma**2)**1.5)

    def get_excess_fe_dimless(self, x):
        """Returns F_ex / (N * kB * T)"""
        if x < 0: return 0.0
        if x > max(self._X_MASTER):
            print(f"Warning: x={x:.2f} out of range ({max(self._X_MASTER)}). Consider calling expand_to().")
        return self._cs_fex(x)

    def get_ideal_gas_fe(self, temp, rho, mass_g_mol):
        m = (mass_g_mol / self.NA) * 1e-3
        h_Js = const.h
        kb_JK = const.k
        Lambda_m = np.sqrt(h_Js**2 / (2 * np.pi * m * kb_JK * temp))
        Lambda_A = Lambda_m * 1e10
        return self.KB * temp * (np.log(rho * (Lambda_A**3)) - 1)

    def get_total_fe(self, temp, rho, sigma, mass_g_mol=None):
        f_ex = self.get_excess_fe_dimless(self.get_x(rho, sigma)) * self.KB * temp
        if mass_g_mol is not None:
            f_ideal = self.get_ideal_gas_fe(temp, rho, mass_g_mol)
            return f_ideal + f_ex
        return f_ex

    def expand_to(self, x_target, lammps_exe, input_file, temp=1000.0, sigma=2.5):
        """
        Dynamically expands the master table by running simulations up to x_target.
        """
        current_max_x = max(self._X_MASTER)
        if x_target <= current_max_x:
            print(f"x_target={x_target} is already within range.")
            return

        print(f"Expanding UF Table: {current_max_x} -> {x_target}")
        
        # Determine points to simulate
        # 1.0 step or smaller if needed
        new_x_points = np.linspace(current_max_x, x_target, int((x_target - current_max_x) / 1.0) + 2)[1:]
        
        ev_ang3_to_bar = 1.6021766208e6
        kBT_bar = self.KB * temp * ev_ang3_to_bar
        s_factor = 0.5 * (np.pi * sigma**2)**1.5

        new_data = [] # (x, z)
        for x in new_x_points:
            log_file = f"log.expand_x_{x}.lammps"
            cmd = [lammps_exe, "-in", input_file, "-var", "x_in", str(x), "-var", "T_in", str(temp), "-log", log_file]
            print(f"Running simulation for x={x}...")
            subprocess.run(cmd, capture_output=True)
            
            # Extract P
            with open(log_file, "r") as f:
                match = re.search(r"DONE: x=[\d\.]+ P_avg=([\d\.\-]+)", f.read())
                if match:
                    p_avg = float(match.group(1))
                    rho = x / s_factor
                    z = p_avg / (rho * kBT_bar)
                    new_data.append((x, z))
                else:
                    print(f"Failed to extract P for x={x}")
        
        if not new_data:
            print("No new data points collected.")
            return

        # Integrate and update internal lists
        # We need a unified list of (x, z) to integrate (Z-1)/x
        # For simplicity, we'll append fex points one by one
        for x_new, z_new in new_data:
            # Trapezoidal integration for the new segment
            x_prev = self._X_MASTER[-1]
            z_prev = self.get_z(x_prev) # Need a way to get Z
            
            # Integrand at prev and new
            integrand_prev = (z_prev - 1.0) / x_prev
            integrand_new = (z_new - 1.0) / x_new
            
            # delta_fex
            delta_fex = 0.5 * (integrand_prev + integrand_new) * (x_new - x_prev)
            fex_new = self._FEX_MASTER[-1] + delta_fex
            
            self._X_MASTER.append(x_new)
            self._FEX_MASTER.append(fex_new)
        
        self.update_spline()
        print(f"Table successfully expanded to x={max(self._X_MASTER)}")

    def get_z(self, x):
        """Returns compressibility factor Z at x (using current splines)"""
        # d(Fex/NkBT)/dx = (Z-1)/x  => Z = 1 + x * d(Fex/NkBT)/dx
        return 1.0 + x * self._cs_fex(x, 1) # First derivative

# Singleton instance
_UF_INSTANCE = UhlenbeckFord(mode='master')

def get_uhlenbeck_ford_fe(temp, rho, sigma):
    return _UF_INSTANCE.get_total_fe(temp, rho, sigma)

def get_ideal_gas_fe(temp, rho, mass_g_mol):
    return _UF_INSTANCE.get_ideal_gas_fe(temp, rho, mass_g_mol)

def expand_uf_table(x_target, lammps_exe, input_file, temp=1000.0, sigma=2.5):
    """Utility to expand the singleton table."""
    _UF_INSTANCE.expand_to(x_target, lammps_exe, input_file, temp, sigma)

def rho_to_x(rho, sigma):
    return UhlenbeckFord.get_x(rho, sigma)

def x_to_rho(x, sigma):
    return UhlenbeckFord.get_rho(x, sigma)

if __name__ == "__main__":
    # Test
    x_test = 4.0
    rho_test = x_to_rho(x_test, 2.5)
    fex = get_uhlenbeck_ford_fe(1000, rho_test, 2.5)
    print(f"Default range test: F_ex(x=4.0) = {fex / (_UF_INSTANCE.KB * 1000):.6f} kBT")
