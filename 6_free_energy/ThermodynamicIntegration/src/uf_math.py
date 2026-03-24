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
    Master Table is managed externally in 'uf_master_table.npz'.
    """
    KB = const.physical_constants["Boltzmann constant in eV/K"][0]
    H = const.physical_constants["Planck constant in eV/Hz"][0]
    NA = const.physical_constants["Avogadro constant"][0]

    def __init__(self, table_path="uf_master_table.npz"):
        self.table_path = table_path
        self._X_MASTER = []
        self._FEX_MASTER = []
        self._cs_fex = None
        self.load_table()

    def load_table(self):
        """Loads the master table from npz file. Fallbacks to initialization if missing."""
        if os.path.exists(self.table_path):
            data = np.load(self.table_path)
            self._X_MASTER = list(data['x_data'])
            self._FEX_MASTER = list(data['fex_data'])
            print(f"Loaded UF Master Table from '{self.table_path}' (x_max={max(self._X_MASTER)})")
        else:
            # If missing, we could throw error or initialize with hardcoded defaults.
            # For robustness, I'll allow an empty start if needed, but usually npz should be there.
            raise FileNotFoundError(f"UF Master Table file '{self.table_path}' not found. Please ensure it exists.")
        
        self.update_spline()

    def save_table(self):
        """Saves current master table to npz file."""
        np.savez(self.table_path, x_data=np.array(self._X_MASTER), fex_data=np.array(self._FEX_MASTER))
        print(f"UF Master Table saved to '{self.table_path}'")

    def update_spline(self):
        if len(self._X_MASTER) > 1:
            self._cs_fex = CubicSpline(self._X_MASTER, self._FEX_MASTER)

    @staticmethod
    def get_x(rho, sigma):
        return 0.5 * (np.pi * sigma**2)**1.5 * rho

    @staticmethod
    def get_rho(x, sigma):
        return x / (0.5 * (np.pi * sigma**2)**1.5)

    def get_excess_fe_dimless(self, x):
        if x < 0: return 0.0
        if x > max(self._X_MASTER):
            print(f"Warning: x={x:.2f} out of range ({max(self._X_MASTER)}). Call expand_to().")
        return self._cs_fex(x)

    def get_total_fe(self, temp, rho, sigma, mass_g_mol=None):
        f_ex = self.get_excess_fe_dimless(self.get_x(rho, sigma)) * self.KB * temp
        if mass_g_mol is not None:
            # Ideal gas calculation
            m = (mass_g_mol / self.NA) * 1e-3
            h_Js = const.h
            kb_JK = const.k
            Lambda_m = np.sqrt(h_Js**2 / (2 * np.pi * m * kb_JK * temp))
            Lambda_A = Lambda_m * 1e10
            f_ideal = self.KB * temp * (np.log(rho * (Lambda_A**3)) - 1)
            return f_ideal + f_ex
        return f_ex

    def expand_to(self, x_target, lammps_exe, input_file, temp=1000.0, sigma=2.5, save=True):
        """Runs simulations to expand the table up to x_target."""
        current_max_x = max(self._X_MASTER)
        if x_target <= current_max_x:
            print(f"x_target={x_target} already in range.")
            return

        print(f"Expanding UF Table: {current_max_x} -> {x_target}")
        new_x_points = np.linspace(current_max_x, x_target, int((x_target - current_max_x) / 1.0) + 2)[1:]
        
        ev_ang3_to_bar = 1.6021766208e6
        kBT_bar = self.KB * temp * ev_ang3_to_bar
        s_factor = 0.5 * (np.pi * sigma**2)**1.5

        for x in new_x_points:
            log_file = f"log.expand_x_{x}.lammps"
            cmd = [lammps_exe, "-in", input_file, "-var", "x_in", str(x), "-var", "T_in", str(temp), "-log", log_file]
            subprocess.run(cmd, capture_output=True)
            
            with open(log_file, "r") as f:
                match = re.search(r"DONE: x=[\d\.]+ P_avg=([\d\.\-]+)", f.read())
                if match:
                    p_avg = float(match.group(1))
                    rho = x / s_factor
                    z_new = p_avg / (rho * kBT_bar)
                    
                    # Integrate
                    x_prev = self._X_MASTER[-1]
                    z_prev = self.get_z(x_prev)
                    delta_fex = 0.5 * ((z_prev-1)/x_prev + (z_new-1)/x) * (x - x_prev)
                    
                    self._X_MASTER.append(x)
                    self._FEX_MASTER.append(self._FEX_MASTER[-1] + delta_fex)
                else:
                    print(f"Failed to expand for x={x}")
        
        self.update_spline()
        if save: self.save_table()

    def get_z(self, x):
        return 1.0 + x * self._cs_fex(x, 1)

# Instance
_UF_INSTANCE = UhlenbeckFord()

def get_uhlenbeck_ford_fe(temp, rho, sigma):
    return _UF_INSTANCE.get_total_fe(temp, rho, sigma)

def rho_to_x(rho, sigma):
    return UhlenbeckFord.get_x(rho, sigma)

def x_to_rho(x, sigma):
    return UhlenbeckFord.get_rho(x, sigma)

if __name__ == "__main__":
    print(f"Current UF Data Range: x [0, {max(_UF_INSTANCE._X_MASTER)}]")