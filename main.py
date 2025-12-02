import sys

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from geopy.geocoders import Nominatim

from utils import (
    J2,
    JD2alphaG,
    clr,
    date2JD,
    dOmegadt_orbit,
    earth_angular_speed_wrt_sun,
    find_la0lo0,
    groundtrack,
    mu,
    np,
    omega_earth,
    os,
    period,
    print_animated_gt,
    print_gt,
    rE,
)

current_file_path = os.path.abspath(__file__)
current_directory = os.path.dirname(current_file_path)
parent_directory = os.path.dirname(current_directory)
sys.path.append(parent_directory)

clr()
JD = date2JD(2025, 11, 25, 21, 30, 0)
alphag = np.deg2rad(JD2alphaG(JD))

# Define altitude and inclination ranges
h = np.linspace(250, 1200, 200)  # Altitude range for the main data
a = [x + rE for x in h]
i = np.linspace(96, 101, 200)  # Inclination range in degrees
ecc = 1e-6

# Calculate values for meshgrid
A, I = np.meshgrid(a, np.radians(i))
P = A * (1 - ecc**2)
T = period(A)
dOmega_dt = dOmegadt_orbit(ecc, A, I)
delta_dOmega_dt = dOmega_dt - earth_angular_speed_wrt_sun()
delta_min, delta_max = np.min(delta_dOmega_dt), np.max(delta_dOmega_dt)
delta_norm = mcolors.Normalize(vmin=delta_min, vmax=delta_max)

h = np.linspace(250, 1200, 200)  # Altitude range for the main data
a = [x + rE for x in h]  # Semi-major axis corresponding to altitude
i = np.linspace(0, 101, 200)  # Inclination range in degrees

city = input("Enter the name of a city (or type custom): ")
if city == "custom":
    La0 = np.deg2rad(float(input("Enter the latitude in degrees: ")))
    Lo0 = np.deg2rad(float(input("Enter the longitude in degrees: ")))
else:
    geolocator = Nominatim(user_agent="geoapi")
    location = geolocator.geocode(city)

    if location:
        print(f"City: {city}")
        print(f"Info: {location}")
        print(f"Latitude: {location.latitude}, Longitude: {location.longitude}")
        Lo0 = np.deg2rad(location.longitude)
        La0 = np.deg2rad(location.latitude)
    else:
        print(f"City {city} not found.")

print("Press to continue")
input()
clr()

print("Let's test our guessing abilities!")
print("You can provide three different orbital altitudes:")
print("Each time I'll plot the groundtrack, and you should improve your guess")
print("so that the groundtrack passes again over the city you provided")
print("First of all, pick an orbit:")
print("[1] SSO")
print("[2] Critical Inclination [prograde, 63.43 deg]")
print("[3] Polar ['prograde', 90 deg]")
print("[4] Custom Inclination")

orbit_type = input(">>> ")
sso = False
if orbit_type == "1":
    sso = True
elif orbit_type == "2":
    i0 = np.arccos(np.sqrt(1 / 5))
elif orbit_type == "3":
    i0 = np.deg2rad(90)
elif orbit_type == "4":
    i0 = np.deg2rad(float(input("Enter the inclination in degrees: ")))
clr()

e0 = 1e-8
rand_ix = []
for tentative in range(3):
    if tentative > 0:
        print("Let's try again!")

    guess_h = float(input(f"Tentative {tentative + 1}: Enter the altitude in km "))
    if guess_h < 0:
        break
    a0 = rE + guess_h
    if sso:
        i0 = dOmegadt_orbit(e0, a=a0, dOmega_dt=earth_angular_speed_wrt_sun())
    Omega0 = alphag + Lo0
    La, Lo, nu0, Omega0 = find_la0lo0(La0, Lo0, a0, e0, i0, Omega0, alphag, 1e-5)
    KE = [a0, e0, i0, Omega0, 0, nu0]
    nTf = np.ceil(86400 / period(a0))
    Lat, Lon, DLo = groundtrack(
        KE, alphag, nTf, 60, dOdt=earth_angular_speed_wrt_sun(), dodt=True
    )
    print_gt(Lat, Lon, period(a0), nTf)

print("Let's work in a different way now.")
print("Press to continue")
input()
clr()

plt.rcParams.update(
    {
        "font.size": 18,
        "axes.titlesize": 20,
        "axes.labelsize": 18,
        "xtick.labelsize": 16,
        "ytick.labelsize": 16,
        "legend.fontsize": 16
    }
)


fig = plt.figure(figsize=(16, 16))
plt.ion()
# Subplot 1: Altitude vs. N_R_revisit_period_day
ax1 = fig.add_subplot(2, 2, 1)
(line1,) = ax1.plot([], [], "o")
ax1.set_xlabel("Altitude (km)")
ax1.set_ylabel(r"$R$, Revisit Period (days)")
ax1.set_title(r"Altitude vs. revisit period $R$")
ax1.grid(which="both")

# Subplot 2: Altitude vs. N_T_orbits_day
ax2 = fig.add_subplot(2, 2, 2)
(line2,) = ax2.plot([], [], "o")
ax2.set_xlabel("Altitude (km)")
ax2.set_ylabel(r"$N_T$, Revolution number")
ax2.set_title(r"Altitude vs. revolution number $N_T$")
ax2.grid(which="both")

# Subplot 3: Altitude vs. tau
ax3 = fig.add_subplot(2, 2, 3)
(line4,) = ax3.plot([], [], "o")
ax3.set_xlabel("Altitude (km)")
ax3.set_ylabel(r"$\tau = N_T / R$")
ax3.set_title(r"Altitude vs. $\tau$")
ax3.grid(which="both")

# 3D plot
ax_3d = fig.add_subplot(2, 2, 4, projection="3d")
scatter_3d = ax_3d.scatter([], [], [])
ax_3d.set_xlabel("Altitude (km)")
ax_3d.set_ylabel(r"$R$, Revisit Period (days)")
ax_3d.set_zlabel(r"$N_T$, Revolution number")
ax_3d.set_title(r"3D Plot: altitude vs. $R$ vs. $N_T$")

h_values, tau_values, Revisit_days, Revolutions, i_points = [], [], [], [], []

saved_points = []
unique_ratios = set()
N_REV_MAX = 20
N_ORBS_DAY_MAX = 300
rev_ix = []
index = 0
total_cases = N_REV_MAX * (N_ORBS_DAY_MAX - 13)
ten_integer_checkpoints = np.linspace(0, total_cases, 11, dtype=int)
for N_R_revisit_period_day in range(1, N_REV_MAX + 1):
    for N_T_orbits_day in range(13, N_ORBS_DAY_MAX + 1):
        index += 1
        tau = N_T_orbits_day / N_R_revisit_period_day
        if tau in unique_ratios:
            continue  # Skip if ratio was already found
        unique_ratios.add(tau)
        h0 = 200
        hf = 1200
        tolerance = 1e-12
        max_iterations = 2000
        for iteration in range(max_iterations):
            h = (h0 + hf) / 2
            a = rE + h
            i = i0
            Omegadot = dOmegadt_orbit(ecc, a=a, i=i)
            if sso:
                i = dOmegadt_orbit(ecc, a=a, dOmega_dt=earth_angular_speed_wrt_sun())
                Omegadot = earth_angular_speed_wrt_sun()
            n = np.sqrt(mu / a**3)
            delta_n = (
                3
                * (rE**2)
                * J2
                * np.sqrt(mu)
                / (4 * a ** (7 / 2))
                * (2 - 3 * np.sin(i) ** 2)
            )
            omegadot = (
                3
                / 4
                * J2
                * (rE**2)
                * np.sqrt(mu)
                / (a ** (7 / 2))
                * (4 - 5 * np.sin(i) ** 2)
            )
            RR = rE**2 * J2 * np.sqrt(mu) / (a ** (7 / 2))

            omega_relative = omega_earth() - Omegadot
            nodal_sc_vel = n + delta_n + omegadot
            target_vel = omega_relative * tau

            tot = nodal_sc_vel - target_vel

            if abs(tot) < tolerance:
                h_values.append(h)
                tau_values.append(tau)
                Revisit_days.append(N_R_revisit_period_day)
                Revolutions.append(N_T_orbits_day)
                i_points.append(i)

                line1.set_data(h_values, Revisit_days)
                line2.set_data(h_values, Revolutions)
                line4.set_data(h_values, tau_values)

                ax_3d.cla()
                ax_3d.scatter(h_values, Revisit_days, Revolutions, c="blue")
                ax_3d.set_xlabel("Altitude (km)")
                ax_3d.set_ylabel(r"$R$, Revisit Period (days)")
                ax_3d.set_zlabel(r"$N_T$, Revolution number")
                ax_3d.set_title(r"3D Plot: altitude vs. $R$ vs. $N_T$")

                for ax in [ax1, ax2, ax3]:
                    ax.relim()
                    ax.autoscale_view()

                if h_values and Revisit_days and Revolutions:
                    ax_3d.set_xlim(min(h_values) - 50, max(h_values) + 50)
                    ax_3d.set_ylim(0, max(Revisit_days))
                    ax_3d.set_zlim(0, max(Revolutions))

                plt.pause(0.01)
                saved_points.append((h, N_T_orbits_day, N_R_revisit_period_day, tau, i))
                break

            hf = h
            if tot > 0:
                h0 = h

# filter h_values for N_R_revisit_period_day only
clr()
N_R_revisit_period_day = int(
    input(f"Filter by desired revisit period (days) (1 - {N_REV_MAX}): ")
)
h_remaining = [point for point in saved_points if point[2] == N_R_revisit_period_day]
h_excluded = [point for point in saved_points if point[2] != N_R_revisit_period_day]

print(f"Remaining cases with R = {N_R_revisit_period_day}: {len(h_remaining)}")

h_values_remaining = [point[0] for point in h_remaining]
nd_values_remaining = [point[2] for point in h_remaining]
ns_values_remaining = [point[1] for point in h_remaining]
tau_values_remaining = [point[3] for point in h_remaining]

h_values_excluded = [point[0] for point in h_excluded]
nd_values_excluded = [point[2] for point in h_excluded]
ns_values_excluded = [point[1] for point in h_excluded]
tau_values_excluded = [point[3] for point in h_excluded]

# Update the existing plots
# Subplot 1: Altitude vs. N_R_revisit_period_day
ax1.cla()
ax1.scatter(
    h_values_remaining,
    nd_values_remaining,
    c="blue",
    label="Included Points",
    alpha=0.8,
    s=20,
)
ax1.scatter(
    h_values_excluded,
    nd_values_excluded,
    c="red",
    label="Excluded Points",
    alpha=0.5,
    s=10,
)
ax1.set_xlabel("Altitude (km)")
ax1.set_ylabel(r"$R$, Revisit Period (days)")
ax1.set_title(r"Altitude vs. revisit period $R$")
ax1.legend(loc="best")
ax1.grid(which="both")

# Subplot 2: Altitude vs. Revolution Number (N_T)
ax2.cla()
ax2.scatter(
    h_values_remaining,
    ns_values_remaining,
    c="blue",
    label="Included Points",
    alpha=0.8,
    s=20,
)
ax2.scatter(
    h_values_excluded,
    ns_values_excluded,
    c="red",
    label="Excluded Points",
    alpha=0.5,
    s=10,
)
ax2.set_xlabel("Altitude (km)")
ax2.set_ylabel(r"$N_T$, Revolution number")
ax2.set_title(r"Altitude vs. revolution number $N_T$")
ax2.legend(loc="best")
ax2.grid(which="both")

# Subplot 3: Altitude vs. Tau
ax3.cla()
ax3.scatter(
    h_values_remaining,
    tau_values_remaining,
    c="blue",
    label="Included Points",
    alpha=0.8,
    s=20,
)
ax3.scatter(
    h_values_excluded,
    tau_values_excluded,
    c="red",
    label="Excluded Points",
    alpha=0.5,
    s=10,
)
ax3.set_xlabel("Altitude (km)")
ax3.set_ylabel(r"$\tau = N_T / R$")
ax3.set_title(r"Altitude vs. $\tau$")
ax3.legend(loc="best")
ax3.grid(which="both")


# Subplot 4: 3D Plot
ax_3d.cla()
ax_3d.scatter(
    h_values_remaining,
    nd_values_remaining,
    ns_values_remaining,
    c="blue",
    label="Included Points",
)
ax_3d.scatter(
    h_values_excluded,
    nd_values_excluded,
    ns_values_excluded,
    c="red",
    label="Excluded Points",
    alpha=0.5,
    s=10,
)
ax_3d.set_xlabel("Altitude (km)")
ax_3d.set_ylabel(r"$R$, Revisit Period (days)")
ax_3d.set_zlabel(r"$N_T$, Revolution number")
ax_3d.set_title(r"3D Plot: altitude vs. $R$ vs. $N_T$")
ax_3d.legend()

plt.draw()
plt.pause(0.01)

print("Filter by desired number of orbits to revisit:")
print("Available options:")
for point in ns_values_remaining:
    print(point, end=", ")
mp = ns_values_remaining[0]
nn = 20
if mp > nn:
    nn = mp
N_T_orbits_day = int(input(": "))
h_remaining = [
    point
    for point in saved_points
    if point[1] == N_T_orbits_day and point[2] == N_R_revisit_period_day
]
h_excluded = [
    point
    for point in saved_points
    if point[1] != N_T_orbits_day or point[2] != N_R_revisit_period_day
]
print(
    f"Remaining altitudes for T_N =  {N_T_orbits_day}, R = {N_R_revisit_period_day}: {len(h_remaining)}"
)
h_values_remaining = [point[0] for point in h_remaining]
nd_values_remaining = [point[2] for point in h_remaining]
ns_values_remaining = [point[1] for point in h_remaining]
tau_values_remaining = [point[3] for point in h_remaining]

h_values_excluded = [point[0] for point in h_excluded]
nd_values_excluded = [point[2] for point in h_excluded]
ns_values_excluded = [point[1] for point in h_excluded]
tau_values_excluded = [point[3] for point in h_excluded]

# Update the existing plots
# Subplot 1: Altitude vs. N_R_revisit_period_day
ax1.cla()
ax1.scatter(
    h_values_remaining,
    nd_values_remaining,
    c="blue",
    label="Included Points",
    alpha=0.8,
    s=20,
)
ax1.scatter(
    h_values_excluded,
    nd_values_excluded,
    c="red",
    label="Excluded Points",
    alpha=0.5,
    s=10,
)
ax1.set_xlabel("Altitude (km)")
ax1.set_ylabel(r"$R$, Revisit Period (days)")
ax1.set_title(r"Altitude vs. revisit period $R$")
ax1.legend(loc="best")
ax1.grid(which="both")

# Subplot 2: Altitude vs. Revolution Number (N_T)
ax2.cla()
ax2.scatter(
    h_values_remaining,
    ns_values_remaining,
    c="blue",
    label="Included Points",
    alpha=0.8,
    s=20,
)
ax2.scatter(
    h_values_excluded,
    ns_values_excluded,
    c="red",
    label="Excluded Points",
    alpha=0.5,
    s=10,
)
ax2.set_xlabel("Altitude (km)")
ax2.set_ylabel(r"$N_T$, Revolution number")
ax2.set_title(r"Altitude vs. revolution number $N_T$")
ax2.legend(loc="best")
ax2.grid(which="both")

# Subplot 3: Altitude vs. Tau
ax3.cla()
ax3.scatter(
    h_values_remaining,
    tau_values_remaining,
    c="blue",
    label="Included Points",
    alpha=0.8,
    s=20,
)
ax3.scatter(
    h_values_excluded,
    tau_values_excluded,
    c="red",
    label="Excluded Points",
    alpha=0.5,
    s=10,
)
ax3.set_xlabel("Altitude (km)")
ax3.set_ylabel(r"$\tau = N_T / R$")
ax3.set_title(r"Altitude vs. $\tau$")
ax3.legend(loc="best")
ax3.grid(which="both")

# Subplot 4: 3D Plot
ax_3d.cla()
ax_3d.scatter(
    h_values_remaining,
    nd_values_remaining,
    ns_values_remaining,
    c="blue",
    label="Included Points",
)
ax_3d.scatter(
    h_values_excluded,
    nd_values_excluded,
    ns_values_excluded,
    c="red",
    label="Excluded Points",
)
ax_3d.set_xlabel("Altitude (km)")
ax_3d.set_ylabel(r"$R$, Revisit Period (days)")
ax_3d.set_zlabel(r"$N_T$, Revolution number")
ax_3d.set_title(r"3D Plot: altitude vs. $R$ vs. $N_T$")
ax_3d.legend()

# Redraw the figure
plt.draw()
plt.pause(0.01)

print("Press to continue")
input()
clr()

h = h_remaining[0][0]
tau = h_remaining[0][3]
T = h_remaining[0][1]
RV = h_remaining[0][2]
nTf = T * RV

print(f"Selected altitude: {h:.3f} km")
print(f"Selected inclination: {np.degrees(i):.3f} deg")
print(f"Selected revisit period: every {RV} day(s)")
print(f"Selected complete orbits to revisit: {nTf}")
print(rf"Selected ratio tau: {tau:.3f}")

a0 = rE + h

dOdt = earth_angular_speed_wrt_sun()
if sso:
    i0 = dOmegadt_orbit(e0, a=a0, dOmega_dt=dOdt)
Omega0 = alphag + Lo0
omega0 = 0
nu0 = 0

err_lat = np.inf
err_lon = np.inf

La0, Lo0, nu0, Omega0 = find_la0lo0(La0, Lo0, a0, e0, i0, Omega0, alphag, 1e-5)
KE = [a0, e0, i0, Omega0, omega0, nu0]
Lat, Lon, DLo = groundtrack(
    KE, alphag, nTf, 60, dOdt=earth_angular_speed_wrt_sun(), dodt=True
)
print(KE, alphag, np.degrees(Lat[0]), np.degrees(Lon[0]))
print_gt(Lat, Lon, period(a0), nTf)
print_animated_gt(Lat, Lon, period(a0), nTf, La0=La0, Lo0=Lo0)
