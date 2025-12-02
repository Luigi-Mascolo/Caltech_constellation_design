import os

import numpy as np
import plotly.graph_objects as go

current_file_path = os.path.abspath(__file__)
current_directory = os.path.dirname(current_file_path)
parent_directory = os.path.dirname(current_directory)

g0 = 9.80665
sb = 5.67e-8
ENABLE_LOGGING = False


def clr():
    os.system("cls" if os.name == "nt" else "clear")


J2 = 1.08263e-3
J4 = -1.58e-6

mu = 398600
mup = {"Earth": mu, "Mars": 4.282837e4, "Sun": 1.32712440018e11, "Moon": 4902.801076}
raR = 149597870.707
rE = 6378.137  # in km
rpl = {"Earth": rE, "Mars": 3389.5, "Venus": 6051.8, "Moon": 1737.4}
rplhe = {"Earth": raR, "Mars": 2.279e8, "Venus": 1.0821e8}


def mean_motion(a, mu=mu):
    return np.sqrt(mu / (a**3))


def mean_anomaly_from_eccentric_anomaly(e, E):
    return E - e * np.sin(E)


def mean_to_true_anomaly(M, e, tol=1e-12, max_iter=100):
    M = wrapToY(M)

    def f(E):
        return E - e * np.sin(E)

    def f_prime(E):
        return 1 - e * np.cos(E)

    E = M  # Approx. initialization
    for i in range(max_iter):
        Eold = E
        E = E + (M - f(E)) / f_prime(E)
        err = E - Eold
        if abs(err) < tol:
            break

    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))

    if nu < 0:
        nu += 2 * np.pi
    return nu


def eccentric_anomaly_from_enu(e, nu):
    E = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(nu / 2))
    return E


def earth_angular_speed_wrt_sun():
    return 2 * np.pi / (86400 * 365.2425)


def omega_earth(flag="solar"):
    omega = 2 * np.pi / 86400
    if flag == "sideral":
        omega = 2 * np.pi / (86164.09)
    return omega


def date2JD(y, m, d, ore, mmin, sec):
    month_correction = np.floor((14 - m) / 12)
    correct_month = m + 12 * month_correction - 3
    correct_year = y + 4800 - month_correction
    correction_leap_year = np.floor(correct_year / 4) + np.floor(
        correct_year / 400
    )
    correction_skip_year = np.floor(correct_year / 100)
    decimal_hour = ore + mmin / 60 + sec / 3600
    day_and_hour = d + (decimal_hour - 12) / 24

    JD = (
        365 * correct_year
        + correction_leap_year
        - correction_skip_year
        + np.floor((153 * correct_month + 2) / 5)
        + day_and_hour
        - 32045
    )
    return JD


def JD2alphaG(JD, JD2000=2451545.0):
    ElapsedDays = JD - JD2000
    T = ElapsedDays / 36525

    ag1 = 280.46061837 + 360.98564736629 * ElapsedDays + 0.0003875 * T**2

    return ag1 % 360


def dOmegadt_orbit(e, a=None, i=None, dOmega_dt=None):
    if e is None:
        raise ValueError("Parameter 'e' must be provided (cannot be None).")
    e = np.asarray(e, dtype=float)
    if a is not None:
        a = np.asarray(a, dtype=float)
    if i is not None:
        i = np.asarray(i, dtype=float)
    if dOmega_dt is not None:
        dOmega_dt = np.asarray(dOmega_dt, dtype=float)

    none_count = int(a is None) + int(i is None) + int(dOmega_dt is None)
    if none_count != 1:
        raise ValueError("Provide exactly 3 of (a, e, i, dOmega_dt).")

    if a is None:   # compute a from dOmega/dt formula
        a = (
            -3 / 2 * J2 * (rE / (1 - e**2)) ** 2 * np.cos(i) * np.sqrt(mu) / (dOmega_dt)
        ) ** (2 / 7)
        return a

    elif i is None:     # compute i from dOmega/dt formula
        i = np.arccos(
            -2
            / 3
            * (dOmega_dt / J2)
            * ((1 - e**2) / (rE)) ** 2
            * np.sqrt((a) ** 7 / (mu))
        )
        return i

    elif dOmega_dt is None:     # compute dOmega/dt
        dOmega_dt = (
            -3
            / 2
            * J2
            * (rE / (1 - e**2)) ** 2
            * np.cos(i)
            * np.sqrt(mu)
            / (a ** (7 / 2))
        )
        return dOmega_dt

def euler_rotation_matrix(angle, axis):
    if axis == "x":
        return np.array(
            [
                [1, 0, 0],
                [0, np.cos(angle), np.sin(angle)],
                [0, -np.sin(angle), np.cos(angle)],
            ]
        )
    elif axis == "y":
        return np.array(
            [
                [np.cos(angle), 0, np.sin(angle)],
                [0, 1, 0],
                [-np.sin(angle), 0, np.cos(angle)],
            ]
        )
    elif axis == "z":
        return np.array(
            [
                [np.cos(angle), np.sin(angle), 0],
                [-np.sin(angle), np.cos(angle), 0],
                [0, 0, 1],
            ]
        )
    else:
        raise ValueError("The axis must be 'x', 'y' or 'z'.")


def eci_to_perifocal(Omega, inc, omega):
    Rz1 = euler_rotation_matrix(Omega, "z")
    Rx = euler_rotation_matrix(inc, "x")
    Rz2 = euler_rotation_matrix(omega, "z")
    rot_matrix = np.dot(np.dot(Rz2, Rx), Rz1)
    return rot_matrix


def perifocal_to_eci(Omega, inc, omega):
    return np.transpose(eci_to_perifocal(Omega, inc, omega))


def eci_to_ecef(GST):
    return euler_rotation_matrix(GST, "z")


def ecef_to_eci(GST):
    return np.transpose(eci_to_ecef(GST))


def eci_to_ecef_state(r_eci, v_eci, GST):
    R = eci_to_ecef(GST)
    r_ecef = np.dot(R, r_eci)
    vel_transpose_earth = np.cross([0, 0, omega_earth()], r_eci)
    v_ecef = np.dot(R, v_eci - vel_transpose_earth)
    return r_ecef, v_ecef


def ecef_to_eci_state(r_ecef, v_ecef, GST):
    R = ecef_to_eci(GST)
    r_eci = np.dot(R, r_ecef)
    vel_transpose_earth = np.cross([0, 0, omega_earth()], r_eci)
    v_eci = np.dot(R, v_ecef) + vel_transpose_earth
    return r_eci, v_eci


def cart2sph(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arctan2(y, x)
    phi = np.arccos(z / r) if r != 0 else 0
    return r, theta, phi


def sph2cart(r, theta, phi):
    x = r * np.cos(phi) * np.cos(theta)
    y = r * np.cos(phi) * np.sin(theta)
    z = r * np.sin(phi)
    return x, y, z


def period(a, mu=mu):
    return 2 * np.pi * np.sqrt(a**3 / (mu))


def deltalon(a, mu=mu):
    T = period(a, mu=mu)
    return omega_earth() * T


def wrapToPi(xin):
    flaglist = False
    if not isinstance(xin, list):
        xin = [xin]
        flaglist = True
    xout = []
    for x in xin:
        xwrap = np.remainder(x, 2 * np.pi)
        mask = np.abs(xwrap) > np.pi
        xwrap = np.where(mask, xwrap - 2 * np.pi * np.sign(xwrap), xwrap)
        mask1 = x < 0
        mask2 = np.remainder(x, np.pi) == 0
        mask3 = np.remainder(x, 2 * np.pi) != 0
        xwrap = np.where(mask1 & mask2 & mask3, xwrap - 2 * np.pi, xwrap)
        xout.append(xwrap)
    if flaglist:
        xout = xout[0]
    return xout


def wrapToY(xin, Y0=0, Y1=2 * np.pi):
    flaglist = False
    if not isinstance(xin, list):
        xin = [xin]
        flaglist = True
    range_length = Y1 - Y0

    xout = []
    for x in xin:
        xwrap = Y0 + np.remainder(x - Y0, range_length)
        xwrap = np.where(xwrap > Y1, xwrap - range_length, xwrap)
        xout.append(xwrap)
    if flaglist:
        xout = xout[0]
    return xout


def groundtrack(Kepler, alphag0, nTf, plot_interval=60, dOdt=0, dodt=False):
    a, e, i, Omega, omega, nu = Kepler
    n = mean_motion(a)
    E = eccentric_anomaly_from_enu(e, nu)

    MA = mean_anomaly_from_eccentric_anomaly(e, E)
    T = period(a)
    time = np.arange(0, nTf * T, plot_interval)  # Time vector
    time = np.append(time, np.array([nTf * T]))
    Lat_rad = np.zeros(len(time))
    Long_rad = np.zeros(len(time))

    for j in range(len(time)):
        alphag = alphag0 + omega_earth() * time[j]
        M = MA + n * time[j]
        nu = mean_to_true_anomaly(M, e, 1e-12)
        Omega_true = Omega + dOdt * time[j]
        omega_true = omega + 3 / 2 * J2 / ((rE / (a * (1 - e**2))) ** 2) * (
            1 - 3 / 2 * np.sin(i) ** 2
        ) * n * np.sqrt(1 - e**2)
        r_ECI, v_ECI = state_vectors_from_orbital_elements(
            a, e, i, Omega_true, omega_true, nu
        )
        r_ECEF, v_ECEF = eci_to_ecef_state(r_ECI, v_ECI, alphag)

        Lat_rad[j], Long_rad[j] = ecef_to_lat_lon(r_ECEF)

    # Convert to degrees
    Lat = np.rad2deg(Lat_rad)
    Long = np.rad2deg(Long_rad)

    return Lat, Long, deltalon(a)


global_fig = None


def print_gt(Lat, Long, T, nTf=1, update=False, save_to_file=False):
    global global_fig

    # Create a new figure if it doesn't exist or update is False
    if global_fig is None or not update:
        global_fig = go.Figure()

    # Add the new data as a trace to the global figure
    global_fig.add_trace(
        go.Scattergeo(
            lon=Long,
            lat=Lat,
            mode="lines",
            line=dict(color="black"),
            name=f"Groundtrack ({nTf} period(s))",
        )
    )
    global_fig.add_trace(
        go.Scattergeo(
            lon=[Long[0]],
            lat=[Lat[0]],
            mode="markers",
            marker=dict(size=10, color="green"),
            name="Start",
        )
    )
    global_fig.add_trace(
        go.Scattergeo(
            lon=[Long[-1]],
            lat=[Lat[-1]],
            mode="markers",
            marker=dict(size=10, color="red"),
            name="End",
        )
    )

    # Update layout only for the first call or when update is False
    if len(global_fig.data) == 3 or not update:
        global_fig.update_layout(
            title=f"Groundtrack - Observation over {nTf:.1f} period(s) (T = {T:.2f} seconds per orbit)",
            geo=dict(
                showland=True,
                subunitcolor="rgb(255, 255, 255)",
                countrycolor="rgb(255, 255, 255)",
                showlakes=True,
                lakecolor="rgb(150, 255, 255)",
                showsubunits=True,
                showcountries=True,
                countrywidth=0.5,
                subunitwidth=0.5,
                resolution=110,
                # Add grid properties
                lonaxis=dict(
                    showgrid=True, gridwidth=0.5, gridcolor="rgb(200, 200, 200)"
                ),
                lataxis=dict(
                    showgrid=True, gridwidth=0.5, gridcolor="rgb(200, 200, 200)"
                ),
            ),
            showlegend=True,
        )

        # add "label" trace
        x = list(range(-180, 180 + 30, 30))
        y = list(range(-90, 90 + 10, 10))
        xpos = -175
        ypos = -85
        global_fig.add_trace(
            go.Scattergeo(
                {
                    "lon": x[1:-1] + [xpos] * (len(y) - 2),
                    "lat": [ypos] * (len(x) - 2) + y[1:-1],
                    "showlegend": False,
                    "text": x[1:-1] + y[1:-1],
                    "mode": "text",
                }
            )
        )

    # SaTe or show the figure
    if save_to_file:
        global_fig.write_html(os.path.join(current_directory, "gt.html"))
        global_fig.write_image(os.path.join(current_directory, "gt.png"), scale=5)
    else:
        global_fig.show()


def print_animated_gt(
    Lat, Long, T, nTf=1, update=False, save_to_file=False, La0=False, Lo0=False
):
    global global_fig

    # Create a new figure for the animation
    fig = go.Figure()

    # Add the ground track animation frames
    frames = []
    for i in range(len(Lat)):
        frames.append(
            go.Frame(
                data=[
                    go.Scattergeo(
                        lon=Long[: i + 1],
                        lat=Lat[: i + 1],
                        mode="lines",
                        line=dict(color="black"),
                        name=f"Groundtrack {nTf}",
                    ),
                    go.Scattergeo(
                        lon=[Long[i]],
                        lat=[Lat[i]],
                        mode="markers",
                        marker=dict(size=10, color="red"),
                        name="Current Position",
                    ),
                ]
            )
        )

    # Initial frame of the animation
    fig.add_trace(
        go.Scattergeo(
            lon=[Long[0]],
            lat=[Lat[0]],
            mode="markers",
            marker=dict(size=10, color="green"),
            name=f"Start {nTf}",
        )
    )

    # Define layout for animation controls and geospatial properties
    fig.update_layout(
        title=f"Groundtrack - Observation over {nTf:.1f} period(s) (T = {T:.2f} seconds)",
        geo=dict(
            showland=True,
            subunitcolor="rgb(255, 255, 255)",
            countrycolor="rgb(255, 255, 255)",
            showlakes=True,
            lakecolor="rgb(150, 255, 255)",
            showsubunits=True,
            showcountries=True,
            countrywidth=0.5,
            subunitwidth=0.5,
            resolution=110,
            center=dict(lon=0),  # Center the map at Greenwich (Longitude 0)
            projection_scale=1,
            # Add grid properties
            lonaxis=dict(showgrid=True, gridwidth=0.5, gridcolor="rgb(200, 200, 200)"),
            lataxis=dict(showgrid=True, gridwidth=0.5, gridcolor="rgb(200, 200, 200)"),
        ),
        updatemenus=[
            {
                "buttons": [
                    {
                        "args": [
                            None,
                            {
                                "frame": {"duration": 50, "redraw": True},
                                "fromcurrent": True,
                            },
                        ],
                        "label": "Play 1x",
                        "method": "animate",
                    },
                    {
                        "args": [
                            None,
                            {
                                "frame": {"duration": 17, "redraw": True},
                                "fromcurrent": True,
                            },
                        ],
                        "label": "3x",
                        "method": "animate",
                    },
                    {
                        "args": [
                            [None],
                            {
                                "frame": {"duration": 0, "redraw": True},
                                "mode": "immediate",
                                "transition": {"duration": 0},
                            },
                        ],
                        "label": "Pause",
                        "method": "animate",
                    },
                ],
                "direction": "left",
                "pad": {"r": 10, "t": 87},
                "showactive": False,
                "type": "buttons",
                "x": 0.1,
                "xanchor": "right",
                "y": 0,
                "yanchor": "top",
            }
        ],
    )

    # add "label" trace
    x = list(range(-180, 180 + 30, 30))
    y = list(range(-90, 90 + 10, 10))
    xpos = -175
    ypos = -85
    fig.add_trace(
        go.Scattergeo(
            {
                "lon": x[1:-1] + [xpos] * (len(y) - 2),
                "lat": [ypos] * (len(x) - 2) + y[1:-1],
                "showlegend": False,
                "text": x[1:-1] + y[1:-1],
                "mode": "text",
            }
        )
    )
    # Add frames to figure
    fig.frames = frames

    if save_to_file:
        current_directory = os.getcwd()  # Adjust if you want a different path
        fig.write_html(os.path.join(current_directory, "animated_gt.html"))
    else:
        fig.show()


def lat_lon_error(a0, e0, i0, Omega0, nu0, alphag):
    # Compute ECI state vectors
    r_ECI, v_ECI = state_vectors_from_orbital_elements(a0, e0, i0, Omega0, 0, nu0)
    # Convert to ECEF
    r_ECEF, _ = eci_to_ecef_state(r_ECI, v_ECI, alphag)
    # Convert ECEF to latitude and longitude
    La0, Lo0 = ecef_to_lat_lon(r_ECEF)
    # Return absolute errors
    return La0, Lo0


def find_la0lo0(La, Lo, a0, e0, i0, Omega0, alphag, force_tol=False):
    retrograde = False
    if i0 > np.pi / 2:
        retrograde = True
    Omega_min, Omega_max = [Lo - np.pi / 2, Lo]
    if retrograde:
        Omega_min, Omega_max = [Lo, Lo + np.pi / 2]
    Omega0 = (Omega_min + Omega_max) / 2

    err = np.inf
    tol = 1e-6
    if retrograde:
        tol = 5e-2
        if force_tol:
            tol = force_tol
    counter = 0
    maxc = 100
    while abs(err) > tol:
        counter += 1
        # Step 1: Find the best nu0 for the current Omega0
        nu_candidates = np.linspace(0, np.pi / 2, 10000)  # Test nu0 in [0, pi/2]
        lats, lons = zip(
            *[lat_lon_error(a0, e0, i0, Omega0, nu, alphag) for nu in nu_candidates]
        )
        for ix, lt in enumerate(lats):
            if lt > La:
                break
        nu0_best = nu_candidates[ix - 1]
        La0 = lats[ix - 1]
        Lo0 = lons[ix - 1]
        min_err = abs(lons[ix - 1] - Lo)

        if Lo0 < Lo:  # If latitude is too high, adjust bounds
            Omega_min = Omega0
        else:
            Omega_max = Omega0

        Omega0 = (Omega_min + Omega_max) / 2  # Update midpoint

        # Update overall error
        err = min_err
        if counter >= maxc:
            break

    return La0, Lo0, nu0_best, Omega0


def ecef_to_lat_lon(r_ECEF):
    x, y, z = r_ECEF
    r_delta = np.linalg.norm(r_ECEF[1:])
    sinA = y / r_delta
    cosA = x / r_delta

    # Calcola latitudine
    lat = np.arcsin(z / np.linalg.norm(r_ECEF))

    # Calcola longitudine
    lon = np.arctan2(sinA, cosA)

    return lat, lon


def state_vectors_from_orbital_elements(a, ecc, inc, Omega, omega, nu, mu=mu):
    """
    Calcola i vettori di posizione e velocità da elementi orbitali di Kepler.

    Args:
    a (float): Semiasse maggiore in metri.
    ecc (float): Eccentricità dell'orbita.
    inc (float): Inclinazione dell'orbita in radianti.
    Omega (float): Longitudine del nodo ascendente in radianti.
    omega (float): Argomento del perigeo in radianti.
    nu (float): Anomalia vera in radianti.

    Returns:
    np.array, np.array: Vettori di posizione e velocità.
    """
    p = a * (1 - ecc**2)
    r = r_conicp(p, ecc, nu)

    r_perifocal = np.array([(r * np.cos(nu)), (r * np.sin(nu)), 0])

    v_perifocal = np.array(
        [-np.sqrt(mu / p) * np.sin(nu), np.sqrt(mu / p) * (ecc + np.cos(nu)), 0]
    )
    R = perifocal_to_eci(Omega, inc, omega)

    r_vec = np.dot(R, r_perifocal)
    v_vec = np.dot(R, v_perifocal)

    return r_vec, v_vec

def r_conicp(p, e, nu):
    return p / (1 + e * np.cos(nu))
