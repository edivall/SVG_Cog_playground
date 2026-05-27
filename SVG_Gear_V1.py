import streamlit as st
import svgwrite
import math
import numpy as np

#time only for checking speed
import time 

#for stl creation ---> pip install trimesh shapely mapbox-earcut
from shapely.geometry import Polygon
import trimesh
import io

import json

# ----------------------------
# Save / Load parameters - as json
# ----------------------------

SAVE_KEYS = [
    "n_teeth",
    "angle_span",
    "tooth_type",
    "r_base",
    "tooth_height",
    "power",
    "width",

    "tooth_reflect",
    "radial_shear",
    "radial_wiggles",
    "wiggle_strength",

    "n_points",

    "n2",
    "slop",
    "steps_per_rev",
    "smooth_iters",
    "rotg1",
    "show_mating",

    "stl_thickness",

    "center_hole",
    "bolt_holes",
    "bolt_radius",
    "bolt_size",

    "spokes",
    "spoke_width",
    "spoke_inner",
    "spoke_outer",
    "dt"
]


def get_save_data():
    """
    Collect all user parameters from session_state.
    """

    data = {}

    for k in SAVE_KEYS:

        if k in st.session_state:
            data[k] = st.session_state[k]

    return data


def load_save_data(data):
    """
    Restore parameters into Streamlit session_state.
    """

    for k, v in data.items():

        if k in SAVE_KEYS:
            st.session_state[k] = v



#Auto change toothspan (can be overridden later manually)

def update_tooth_span():
    st.session_state.angle_span = 360 / st.session_state.n_teeth
# ----------------------------
# x-y graphs
# ----------------------------


def build_debug_data_from_points(polar_points):
    n = len(polar_points)

    angles = []
    radius_vals = []
    theta_vals = []
    tweak_vals = []

    for i, (r, theta) in enumerate(polar_points):
        t = i / (n - 1) if n > 1 else 0

        angles.append(t * 360)  # normalized display angle
        radius_vals.append(r)
        theta_vals.append(theta)

    # compute "tweak" as deviation from linear angle span
    base_theta = theta_vals[0]
    for theta in theta_vals:
        tweak_vals.append(theta - base_theta)

    return angles, radius_vals, tweak_vals

# ----------------------------
# Geometry helpers
# ----------------------------

def rotate_point(x, y, angle_deg):
    a = math.radians(angle_deg)
    return (
        x * math.cos(a) - y * math.sin(a),
        x * math.sin(a) + y * math.cos(a),
    )


def polar_to_cartesian(r, theta_deg):
    t = math.radians(theta_deg)
    return (r * math.cos(t), r * math.sin(t))


# ----------------------------
# Tooth profiles
# ----------------------------
def radius_profile(t, r_base, tooth_height, power, width, tooth_type,nteeth): #t=0 to 1
    if t>1:
        return r_base    
    
    if tooth_type == "Gaussian":
        x = (t - 0.5) / width
        f = math.exp(-(x ** 2)) ** power
        return r_base + tooth_height * f

    elif tooth_type == "Sinusoidal":
        s = math.sin(2 * math.pi * t - math.pi / 2)
        f = ((s + 1) / 2) ** power
        return r_base + tooth_height * f

    elif tooth_type == "Spike/Square":
        # -----------------------------
        # Pure linear trapezoid tooth
        # -----------------------------

        phase = t % 1.0
        center = 0.5

        half_width = width / 2.0

        # sharpness controls slope region size
        slope_width = max(1e-6, half_width / max(power, 1e-6))

        d = abs(phase - center)

        r0 = r_base
        r1 = r_base + tooth_height

        # outside tooth
        if d >= half_width:
            return r0

        # flat top region
        elif d <= (half_width - slope_width):
            return r1

        # linear slope region
        else:
            x = (d - (half_width - slope_width)) / slope_width
            x = min(max(x, 0.0), 1.0)

            # straight-line interpolation
            return r1 - x * (r1 - r0)
            
    elif tooth_type == "Rounded Square": #use power as radius
        # 0                  
        # 0___radpos1----¦cen-hw---radpos2--¦cen
        cirdist=2 * 3.1415 * r_base / nteeth  #approx width of tooth 
        center=0.5
        rad1=power/10 * (1 -width)
        rad2=power/10 * width 
        half_width = width / 2.0 
        radpos1=center-half_width-rad1
        radpos2=center-half_width+rad2
        if radpos2>center:
            radpos2=center
        th=tooth_height  
        th1=th *2.0
        th2=th *2.0   # * (1-width) / (width)
        if t>0.5:
            t=1-t
        safety = 1.02     # stop squareroot -1
        if t<=(radpos1):
            return r_base 
        
        elif  t<(center-half_width): 
            return r_base + th1*(rad1-(rad1 * rad1 *safety - (radpos1 - t) ** 2) ** 0.5)
            
        elif  t<(radpos2): 
            return r_base + tooth_height -th2*(rad2-(rad2 * rad2 *safety- (radpos2 - t) ** 2) ** 0.5 )     

        else:
            return r_base+tooth_height
     
    elif tooth_type == "Ratchet":
        return r_base + t * tooth_height
     
    else:
        return r_base


def angle_profile(t, angle_span, rfrac, rad_wig, wig_strength, reflect, radial_shear):
    theta = t * angle_span 
    if (reflect == "On") and (t> 0.5):
        wig=-1.0
        t=1-t
    else:
        wig=1.0
     
    if rfrac>1: rfrac=1
    if rfrac<=0: rfrac=0.0001    
    cycles=8
    dist  =3
    s =  math.sin(cycles    * math.pi * t) #t is 0 to 1
    
    #mag =s
    #if (t < (dist/cycles)) or (t > ((cycles-dist)/cycles)):
    #    mag = (abs(s)) ** angle_power
    #    mag =theta + angle_wave * mag * math.copysign(1, s)
    #else:
    mag=theta
  
    mag=mag + math.sin(rad_wig    * math.pi * rfrac) * wig_strength*wig + rfrac * radial_shear
    
    #print(mag,angle_power)
    #mag = mag - 0.5
    return mag


# ----------------------------
# Tooth generator
# ----------------------------

def generate_tooth(params):

    pts = []
    fullang=360/params["n_teeth"]
    
    for i in range(params["n_points"] + 1):
        
        t =  i /params["n_points"] 
        t= t * fullang / params["angle_span"] 
        #widthFract=params["width"]/params["angle_span"]
        
        r = radius_profile(
            t,
            params["r_base"],
            params["tooth_height"],
            params["power"],
            params["width"]/100,
            params["tooth_type"],
            params["n_teeth"]
        )
        
        if params["tooth_height"]==0:
            rfrac=1000
        else:    
            rfrac= (r-params["r_base"])/(params["tooth_height"])
        
        theta = angle_profile(
            t,
            params["angle_span"],
            #params["angle_wave"],
            #params["angle_power"],
            rfrac,
            params["radial_wiggles"],
            params["wiggle_strength"],
            params["tooth_reflect"],
            params["radial_shear"]
        )

        pts.append((r, theta))

    return pts


# ----------------------------
# Outer gear
# ----------------------------

def build_outer_profile(polar_points_deg, n_copies, center):
    cx, cy = center

    base = [
        (r * math.cos(math.radians(t)), r * math.sin(math.radians(t)))
        for r, t in polar_points_deg
    ]

    pts = []

    for i in range(n_copies):
        ang = 360 * i / n_copies
        rot = [rotate_point(x, y, ang) for x, y in base]
        pts.extend([(x + cx, y + cy) for x, y in rot])

    pts.append(pts[0])
    return pts


# ----------------------------
# Spokes
# ----------------------------

def spoke_wedge(cx, cy, r_inner, r_outer, angle_deg, width_deg, dt):

    a1 = angle_deg - width_deg / 2
    a2 = angle_deg + width_deg / 2

    p1 = polar_to_cartesian(r_inner, a1)
    p2 = polar_to_cartesian(r_outer, a1 - dt)
    p3 = polar_to_cartesian(r_outer, a2 + dt)
    p4 = polar_to_cartesian(r_inner, a2)

    p1 = (p1[0] + cx, p1[1] + cy)
    p2 = (p2[0] + cx, p2[1] + cy)
    p3 = (p3[0] + cx, p3[1] + cy)
    p4 = (p4[0] + cx, p4[1] + cy)

    return f"""
    M {p1[0]},{p1[1]}
    L {p2[0]},{p2[1]}
    A {r_outer},{r_outer} 0 0,1 {p3[0]},{p3[1]}
    L {p4[0]},{p4[1]}
    A {r_inner},{r_inner} 0 0,0 {p1[0]},{p1[1]}
    Z
    """


# ----------------------------
# SVG builder
# ----------------------------

def build_svg(polar_points, p):

    size = svg_size
    cx, cy = size / 2, size / 2
    #svg for display
    dwg = svgwrite.Drawing(size=(size, size), viewBox=f"0 0 {size} {size}")
    #svg to save
    r_outer = max(p["r_base"], p["r_base"] + p["tooth_height"])
    size = 2 * r_outer
    sdwg = svgwrite.Drawing(
        size=(f"{size}mm", f"{size}mm"),
        viewBox=f"0 0 {size} {size}"
    )
    outer = build_outer_profile(polar_points, p["n_teeth"], (cx, cy))
    path_d = "M " + " L ".join(f"{x},{y}" for x, y in outer) + " Z"

    # center hole
    r = p["center_hole"]
    path_d += f"""
    M {cx+r},{cy}
    A {r},{r} 0 1,0 {cx-r},{cy}
    A {r},{r} 0 1,0 {cx+r},{cy}
    """

    # bolt holes
    for i in range(p["bolt_holes"]):
        ang = 360 * i / p["bolt_holes"]
        x, y = polar_to_cartesian(p["bolt_radius"], ang)
        x += cx
        y += cy

        br = p["bolt_size"]

        path_d += f"""
        M {x+br},{y}
        A {br},{br} 0 1,0 {x-br},{y}
        A {br},{br} 0 1,0 {x+br},{y}
        """

    # spokes
    r_inner = p["center_hole"] + p["spoke_inner"]
    r_outer = p["bolt_radius"] - p["bolt_size"] - p["spoke_outer"]

    for i in range(p["spokes"]):
        ang = 360 * i / p["spokes"]
        path_d += spoke_wedge(
            cx, cy,
            r_inner, r_outer,
            ang,
            p["spoke_width"],
            p["dt"]
        )

    dwg.add(dwg.path(
        d=path_d,
        fill="lightsteelblue",
        stroke="black",
        stroke_width=1,
        fill_rule="evenodd"
    ))
    
    sdwg.add(sdwg.path(
        d=path_d,
        fill="lightsteelblue",
        stroke="black",
        stroke_width=1,
        fill_rule="evenodd"
    ))

    return dwg.tostring(),sdwg.tostring()


# ----------------------------
# Streamlit UI
# ----------------------------

# Save / Load streamlit UI parameters


save_data = json.dumps(
    get_save_data(),
    indent=2
)



st.title("⚙️ SVG Cog / Gear shape Generator!")

with st.sidebar:

    st.header("Settings")
    
    st.download_button(
        "💾 Save Current Configuration",
        data=save_data,
        file_name="gear_preset.json",
        mime="application/json"
    )

    uploaded_files = st.file_uploader(
        "📂 Reload a Gear Configuration (multiple allowed)",
        type=["json"],
        #label_visibility="collapsed",
        accept_multiple_files=True
    )

    if uploaded_files:
        st.write("### Loaded saved files")

        for i, f in enumerate(uploaded_files):
            try:
                # IMPORTANT: read file safely
                data = json.loads(f.getvalue().decode("utf-8"))
            except Exception as e:
                st.error(f"Failed to load {f.name}: {e}")
                continue

            col1, col2 = st.columns([4, 1])

            with col1:
                st.write(f"📄 {f.name}")

            with col2:
                if st.button("Use file", key=f"use_preset_{i}_{f.name}"):
                    load_save_data(data)
                    st.success(f"Loaded: {f.name}")
                    st.rerun()

        
    
    
    
    
    svg_size = st.slider("SVG image size", 50, 800, 300, key="svg_size")

    st.header("Gear")
    n_teeth = st.slider("Teeth", 3, 100, 20, on_change=update_tooth_span, key="n_teeth")
    angle_span = st.slider("Tooth span (deg)", 2.0, 130.0, 18.0, key="angle_span")
    tooth_type = st.selectbox("Profile", ["Sinusoidal", "Spike/Square", "Gaussian", "Rounded Square", "Ratchet"], key="tooth_type")
    r_base = st.slider("Base radius", 1.0, 120.0, 50.0, key="r_base")
    tooth_height = st.slider("Tooth height", -50.0, 50.0, 10.0, key="tooth_height")
    power = st.slider("Sharpness (1 for pure gauss)", 0.1, 7.0, 1.0, key="power")
    width = st.slider("Tooth width (not sinusoidal)", 1.0, 99.0, 50.0, key="width")

    st.header("Radial wave")
    tooth_reflect = st.selectbox("Mid point reflection", ["Off", "On"], key="tooth_reflect")
    radial_shear = st.slider("Radial shear (spiral)", -15.0, 15.0, 0.0, key="radial_shear")
    radial_wiggles = st.slider("Radial wiggles (half oscillations)", -5.0, 5.0, 0.0, key="radial_wiggles")
    wiggle_strength = st.slider("Radial Wiggle amplitude (0=Off)", -6.0, 6.0, 2.0, key="wiggle_strength")

    st.header("Samples")
    n_points = st.slider("Resolution (Samples per tooth on gear 1)", 1, 250, 100, key="n_points")

    st.header("Mating Gear")
    n2 = st.slider("Gear 2 teeth", 2, 150, 10, key="n2")
    slop = st.slider("Gap between", 0.0, 20.0, 0.5, key="slop")
    steps_per_rev = st.slider("Samples per complete gear (360 deg)", 180, 2880, 720, key="steps_per_rev")
    smooth_iters = st.slider("Smoothing iterations for gear 2", 0, 10, 0, key="smooth_iters")
    rotg1 = st.number_input("rotate gear1", 0.0, 90.0, 0.0, 1.0, key="rotg1")
    show_mating = st.checkbox("Show mating gear", True, key="show_mating")

    st.header("3D Export")
    stl_thickness = st.slider("STL thickness", 1.0, 50.0, 5.0, key="stl_thickness")

    st.header("Holes / Spokes - SVG save only")
    center_hole = st.slider("Center hole (0=Off)", 0, 50, 5, key="center_hole")
    bolt_holes = st.slider("Outer holes (0=Off)", 0, 120, 6, key="bolt_holes")
    bolt_radius = st.slider("Outer hole radius", 20, 150, 70, key="bolt_radius")
    bolt_size = st.slider("Outer hole size", 1, 15, 3, key="bolt_size")

    spokes = st.slider("Spokes (0=Off)", 0, 12, 4, key="spokes")
    spoke_width = st.slider("Spoke hole width ( < 360 / num spokes)", 5, 120, 60, key="spoke_width")
    spoke_inner = st.slider("Spoke inner margin", 0, 50, 10, key="spoke_inner")
    spoke_outer = st.slider("Spoke outer margin", 0, 50, 10, key="spoke_outer")
    dt = st.slider("Spoke width tweak (dt)", 0.0, 20.0, 5.0, key="dt")
    

params = {
    "svg_size": svg_size,
    "r_base": r_base,
    "tooth_height": tooth_height,
    "power": power,
    "width": width,
    "tooth_type": tooth_type,
    "tooth_reflect": tooth_reflect,
    #"angle_wave": angle_wave,
    #"angle_power": angle_power,
    "radial_shear": radial_shear,
    "radial_wiggles": radial_wiggles,
    "wiggle_strength": wiggle_strength,
    "angle_span": angle_span,
    "n_points": n_points,

    "n_teeth": n_teeth,
    "center_hole": center_hole,
    "bolt_holes": bolt_holes,
    "bolt_radius": bolt_radius,
    "bolt_size": bolt_size,
    "spokes": spokes,
    "spoke_width": spoke_width,
    "spoke_inner": spoke_inner,
    "spoke_outer": spoke_outer,
    "dt": dt,
    "slop": slop,
    "rotg1": rotg1
}


# Gear 2
def build_gear1_points(polar_points_deg, n_teeth, rotg1):
    """
    Returns full list of (x, y) points for gear 1 in Cartesian coords,
    centered at (0,0), NOT shifted to SVG center.
    """

    # convert one tooth from polar → cartesian
    base = [
        (
            r * math.cos(math.radians(theta+rotg1)),
            r * math.sin(math.radians(theta+rotg1))
        )
        for r, theta in polar_points_deg
    ]

    pts = []

    # replicate around circle
    for i in range(n_teeth):
        ang = 360 * i / n_teeth
        rot = [rotate_point(x, y, ang) for x, y in base]
        pts.extend(rot)

    return pts







def generate_mating_gear_level1_works_but_could_be_faster(gear1_pts, n1, n2, d, slop=0.0, steps=120):

    ratio = n1 / n2
    rots =1 
    if n2>n1:
        rots=n2/n1

    # --- envelope storage: one best point per angle bin
    bins = {}

    def ang(x, y):
        return math.atan2(y, x)

    gear1_pts = np.array(gear1_pts)
    
    for phi in np.linspace(0, rots*2*np.pi, steps):

        phi2 = -phi * ratio

        c1, s1 = math.cos(phi), math.sin(phi)
        c2, s2 = math.cos(phi2), math.sin(phi2)

        for x, y in gear1_pts:

            # --- rotate gear1
            x1 = x * c1 - y * s1
            y1 = x * s1 + y * c1

            # --- translate to gear2 center
            x1 = d -x1

            # --- rotate into gear2 frame
            x2 = x1 * c2 - y1 * s2
            y2 = x1 * s2 + y1 * c2

            # --- bin by angle
            a = ang(x2, y2)

            # 720 bins gives smooth but stable outline
            k = int((a + math.pi) / (2 * math.pi) * 720)
            

            r2 = x2*x2 + y2*y2
            #r2 = x2 * math.cos(a) + y2 * math.sin(a)

            # keep only outer envelope (max radius per direction)

            if k not in bins or r2 < bins[k][0]:
                bins[k] = (r2, (x2, y2))

    # --- extract final contour
    gear2 = [v[1] for v in bins.values()]

    # --- sort into continuous loop
    gear2.sort(key=lambda p: math.atan2(p[1], p[0]))

    if slop != 0.0:
        new_gear2 = []

        for x, y in gear2:
            r = math.hypot(x, y)

            if r > 1e-9:
                r_new = max(r - slop, 1e-6)  # prevent collapse
                scale = r_new / r
                new_gear2.append((x * scale, y * scale))
            else:
                new_gear2.append((x, y))

        gear2 = new_gear2



    return gear2

def generate_mating_gear_level1(gear1_pts, n1, n2, d, slop=0.0, steps=120, nbins=720):

    """
    Faster vectorized mating gear generator.
    gear1_pts : [(x,y), ...]        Cartesian points of gear1 centered at origin
    n1, n2 : int       Tooth counts
    d : float          Center distance between gears
    slop : float       Radial clearance
    steps : int        Rotation samples
    nbins : int        Angular envelope resolution
    """

    # ---------------------------------
    # gear ratio
    # ---------------------------------
    ratio = n1 / n2

    # ---------------------------------
    # preserve original behavior
    # ---------------------------------
    rots = 1    #number of rotations to make

    if n2 > n1:
        rots = n2 / n1

    # ---------------------------------
    # convert once to numpy
    # ---------------------------------
    pts = np.asarray(gear1_pts, dtype=np.float64)

    x = pts[:, 0]
    y = pts[:, 1]

    # ---------------------------------
    # envelope storage
    # ---------------------------------
    best_r2 = np.full(nbins, np.inf)

    best_x = np.zeros(nbins)
    best_y = np.zeros(nbins)

    # ---------------------------------
    # precompute rotations
    # ---------------------------------
    phis = np.linspace(
        0,
        rots * 2*np.pi,
        steps,
        endpoint=False
    )

    cos1 = np.cos(phis)
    sin1 = np.sin(phis)

    phi2s = -phis * ratio

    cos2 = np.cos(phi2s)
    sin2 = np.sin(phi2s)

    # ---------------------------------
    # main loop
    # ---------------------------------
    for i in range(steps):

        c1 = cos1[i]
        s1 = sin1[i]

        c2 = cos2[i]
        s2 = sin2[i]

        # -----------------------------
        # rotate gear1
        # -----------------------------
        x1 = x * c1 - y * s1
        y1 = x * s1 + y * c1

        # -----------------------------
        # original mirrored transform
        # -----------------------------
        x1 = d - x1

        # -----------------------------
        # rotate into gear2 frame
        # -----------------------------
        x2 = x1 * c2 - y1 * s2
        y2 = x1 * s2 + y1 * c2

        # -----------------------------
        # angular bins
        # -----------------------------
        ang = np.arctan2(y2, x2)

        k = (
            ((ang + np.pi) / (2*np.pi) * nbins)
            .astype(np.int32)
        ) % nbins

        # -----------------------------
        # radius²
        # -----------------------------
        r2 = x2*x2 + y2*y2

        # -----------------------------
        # envelope extraction
        # keep INNER envelope
        # -----------------------------
        for j in range(len(r2)):

            kk = k[j]

            if r2[j] < best_r2[kk]:

                best_r2[kk] = r2[j]

                best_x[kk] = x2[j]
                best_y[kk] = y2[j]

    # ---------------------------------
    # build contour
    # ---------------------------------
    gear2 = []

    for i in range(nbins):

        if best_r2[i] != np.inf:

            xx = best_x[i]
            yy = best_y[i]

            # -------------------------
            # apply radial slop
            # -------------------------
            if slop != 0.0:

                r = math.hypot(xx, yy)

                if r > 1e-9:

                    r_new = max(r - slop, 1e-6)

                    s = r_new / r

                    xx *= s
                    yy *= s

            gear2.append((xx, yy))

    # ---------------------------------
    # continuous polygon ordering
    # ---------------------------------
    gear2.sort(
        key=lambda p: math.atan2(p[1], p[0])
    )

    return gear2

# option to smooth gear 2 teeth
def chaikin_smooth(points, iterations=2):       #Test 1  - not that effective
    pts = np.array(points)

    for _ in range(iterations):
        new_pts = []

        for i in range(len(pts)):
            p0 = pts[i]
            p1 = pts[(i + 1) % len(pts)]

            Q = 0.75 * p0 + 0.25 * p1
            R = 0.25 * p0 + 0.75 * p1

            new_pts.extend([Q, R])

        pts = np.array(new_pts)

    return pts.tolist()

def smooth_radius(points, alpha=0.2, iterations=2):
    pts = np.array(points)

    x = pts[:, 0]
    y = pts[:, 1]

    r = np.sqrt(x*x + y*y)
    theta = np.arctan2(y, x)

    # ensure sorted by angle
    order = np.argsort(theta)
    r = r[order]
    theta = theta[order]

    # circular smoothing
    for _ in range(iterations):
        r_new = r.copy()

        for i in range(len(r)):
            r_prev = r[i - 1]
            r_next = r[(i + 1) % len(r)]

            r_new[i] = (1 - alpha) * r[i] + alpha * 0.5 * (r_prev + r_next)

        r = r_new

    # rebuild cartesian
    x_new = r * np.cos(theta)
    y_new = r * np.sin(theta)

    return list(zip(x_new, y_new))


def build_svg_two_gears(gear1_pts, gear2_pts, d, r1, svg_size, hole_r ):

    cx, cy = svg_size / 2 - r1, svg_size / 2    

    dwg = svgwrite.Drawing(size=(svg_size, svg_size),
                           viewBox=f"0 0 {svg_size} {svg_size}")

    # --- Gear 1 (centered)
    pts1 = [(x + cx, y + cy) for (x, y) in gear1_pts]
    if pts1:
        path1 = "M " + " L ".join(f"{x},{y}" for x, y in pts1) + " Z"

        r = hole_r
        path1 += f"""
        M {cx+r},{cy}
        A {r},{r} 0 1,0 {cx-r},{cy}
        A {r},{r} 0 1,0 {cx+r},{cy}
        """
        
        dwg.add(dwg.path(
            d=path1,
            fill="lightsteelblue",
            stroke="black",
            stroke_width=1,
            fill_rule="evenodd"
        ))

    # --- Gear 2 (shifted by center distance)
    pts2 = [(cx + d - x, cy + y) for (x, y) in gear2_pts]
    if pts2:
        path2 = "M " + " L ".join(f"{x},{y}" for x, y in pts2) + " Z"

        r = hole_r
        path2 += f"""
        M {cx+r+d},{cy}
        A {r},{r} 0 1,0 {cx-r+d},{cy}
        A {r},{r} 0 1,0 {cx+r+d},{cy}
        """

        dwg.add(dwg.path(
            d=path2,
            fill="lightcoral",
            stroke="black",
            stroke_width=1,
            fill_rule="evenodd"
        ))

    return dwg.tostring()




# ----------------------------
# STL Export
# ----------------------------
def make_circle_points(radius, segments=64, center=(0, 0), reverse=True):

    cx, cy = center

    pts = []

    for i in range(segments):
        a = 2 * math.pi * i / segments

        x = cx + radius * math.cos(a)
        y = cy + radius * math.sin(a)

        pts.append((x, y))

    if reverse:
        pts.reverse()

    return pts


def points_to_stl_bytes(gear1_pts, gear2_pts=None, thickness=5.0, mirror_gear=0, spacing_shift=(0, 0), center_hole_radius=0):

    meshes = []

    # ------------------------
    # Gear 1 with center hole
    # ------------------------

    holes = []

    if center_hole_radius > 0:
        hole_pts = make_circle_points(center_hole_radius)
        holes.append(hole_pts)
    
    if mirror_gear==0:    
        poly1 = Polygon(gear1_pts, holes=holes)
    else:
       mirrored = [(-x, y ) for x, y in gear1_pts]
       poly1 = Polygon(mirrored, holes=holes)
       

    if not poly1.is_valid:
        poly1 = poly1.buffer(0)

    mesh1 = trimesh.creation.extrude_polygon(
        poly1,
        height=thickness
    )

    meshes.append(mesh1)

    # ------------------------
    # Gear 2
    # ------------------------
      
    if gear2_pts and len(gear2_pts) > 2:

        dx, dy = spacing_shift

        # IMPORTANT:
        # mirror X exactly like SVG rendering
        shifted = [(-x + dx, y + dy) for x, y in gear2_pts]

        holes2 = []

        if center_hole_radius > 0:
            hole2 = make_circle_points(
                center_hole_radius,
                center=(dx, dy)
            )
            holes2.append(hole2)

        poly2 = Polygon(shifted, holes=holes2)

        if not poly2.is_valid:
            poly2 = poly2.buffer(0)

        mesh2 = trimesh.creation.extrude_polygon(
            poly2,
            height=thickness
        )    
              

        meshes.append(mesh2)

    # ------------------------
    # Combine meshes
    # ------------------------

    combined = trimesh.util.concatenate(meshes)

    return combined.export(file_type='stl')








def normalize_to_pitch(gear_pts, target_radius):
    out = []
    for x, y in gear_pts:
        r = math.hypot(x, y) + 1e-9
        s = target_radius / r
        out.append((x * s, y * s))
    return out
# ----------------------------
# Generate + render
# ----------------------------

polar_points = generate_tooth(params)
gear1_pts = build_gear1_points(polar_points, params["n_teeth"], rotg1)  #make full gear in cartesian coords

n1 = params["n_teeth"]
n2 = n2  # your slider / choice

# --- estimate outer radius of gear1
r_vals = [r for r, _ in polar_points]

R1_pitch = (    max(r_vals) + min(r_vals)) / 2

# --- enforce correct gear ratio constraint
R2_pitch = R1_pitch * (n2 / n1)

# --- center distance (true gear law)
R1_outer = max(r_vals)
R1_inner = min(r_vals)
R1_pitch = (R1_outer + R1_inner) / 2

R2_pitch = R1_pitch * (n2 / n1)

d = R1_pitch + R2_pitch


#gear1_pts = normalize_to_pitch(gear1_pts, R1_pitch)

slop=params["slop"]

#t0 = time.perf_counter()
if show_mating:
   
    gear2_pts = generate_mating_gear_level1(gear1_pts, n1, n2, d,slop, steps=120, nbins=steps_per_rev ) # over 4 times faster
    
    
    if smooth_iters!=0 and gear2_pts:
        #gear2_pts = chaikin_smooth(gear2_pts, iterations=smooth_iters)
        gear2_pts = smooth_radius(gear2_pts, alpha=0.3, iterations=smooth_iters)
    #gear2_pts = generate_mating_gear_level1_works_but_could_be_faster(gear1_pts, n1, n2, d,slop, steps=120  )
else:
    gear2_pts = []


#t1 = time.perf_counter()
#print(f"generate_mating_gear_level1: {(t1 - t0):.4f} sec")

#angles, radius_vals, tweak_vals = build_debug_data(params)  #xy graphs
angles, radius_vals, tweak_vals = build_debug_data_from_points(polar_points)

#svg_str,ssvg_str = build_svg(polar_points, params)
if show_mating and gear2_pts:
    svg_str = build_svg_two_gears(gear1_pts, gear2_pts, d, R1_outer, svg_size, params["center_hole" ])
    ssvg_str = svg_str
else:
    svg_str, ssvg_str = build_svg(polar_points, params)
    
    

with st.container():
    st.image(svg_str, use_container_width=True)
    
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("📏 Radius vs Angle")
        st.line_chart({"radius": radius_vals}, height=150)

    with col2:
        st.subheader("🔄 Angle tweak vs Angle")
        st.line_chart({"tweak": tweak_vals}, height=150)








st.download_button(
    "💾 Download SVG",
    ssvg_str,
    file_name="gear.svg",
    mime="image/svg+xml"
)

st.download_button(
    "💾 Download SVG - Screenshot",
    svg_str,
    file_name="gear.svg",
    mime="image/svg+xml"
)

# ----------------------------
# STL generation - in one file
# ----------------------------

if st.button("🧊 Generate single STL file"):

    with st.spinner("Generating STL..."):

        if show_mating and gear2_pts:

            gear2_shift_x = d

            stl_data = points_to_stl_bytes(
                gear1_pts,
                gear2_pts,
                thickness=stl_thickness,
                mirror_gear=0,
                spacing_shift=(gear2_shift_x, 0),
                center_hole_radius=params["center_hole"]
            )
            

        else:

            stl_data = points_to_stl_bytes(
                gear1_pts,
                None,
                thickness=stl_thickness,
                mirror_gear=0,
                center_hole_radius=params["center_hole"]
            )

    st.download_button(
        "💾 Download STL",
        data=stl_data,
        file_name="gears.stl",
        mime="model/stl"
    )

# ----------------------------
# STL generation - separate files
# ----------------------------
def clear_stl_downloads():  #Clear save buttons after gear 2 saved
    st.session_state.gear1_stl = None
    st.session_state.gear2_stl = None

if "gear1_stl" not in st.session_state:
    st.session_state.gear1_stl = None

if "gear2_stl" not in st.session_state:
    st.session_state.gear2_stl = None


if st.button("🧊 Generate 2 STL files"):

    with st.spinner("Generating STL..."):

        # ------------------------
        # Gear 1 STL
        # ------------------------

        st.session_state.gear1_stl = points_to_stl_bytes(
            gear1_pts,
            thickness=stl_thickness,
            mirror_gear=0,
            center_hole_radius=params["center_hole"]
        )

        # ------------------------
        # Gear 2 STL
        # ------------------------

        st.session_state.gear2_stl = None

        if show_mating and gear2_pts:

            st.session_state.gear2_stl = points_to_stl_bytes(
                gear2_pts,
                thickness=stl_thickness,
                mirror_gear=1,
                center_hole_radius=params["center_hole"]
            )

# ------------------------
# Persistent download buttons
# ------------------------

if st.session_state.gear1_stl is not None:

    col1, col2 = st.columns(2)

    with col1:

        st.download_button(
            "💾 Download Gear 1 STL",
            data=st.session_state.gear1_stl,
            file_name="gear1.stl",
            mime="model/stl"
        )

    with col2:

        if st.session_state.gear2_stl is not None:

            st.download_button(
                "💾 Download Gear 2 STL",
                data=st.session_state.gear2_stl,
                file_name="gear2.stl",
                mime="model/stl",
                on_click=clear_stl_downloads
            )







#impath="Gears\"
impath="https://raw.githubusercontent.com/edivall/SVG_Cog_playground/main/images/"

# Local data array for user info
items = [
    {
        "title": "Zeroth Item",
        "image": impath+"GearMesh4.png",
        "caption": "Primary and meshing gears",
        "header": "Meshing gear calc",
        "description": "Given a primary gear (blue/grey) a second (meshing) gear (red) is calculated - as best as possible! If it has a noisy profile then try increasing samples on the first tooth.",
        "points": ["Number of teeth on meshing gear can be set.", "Additional gap between the 2 gears can be added.", "Primary gear can be rotated to see meshing in action!"]
    },
    {
        "title": "First Item",
        "image": impath+"GearSin1.png",
        "caption": "Sinusoidal",
        "header": "Sinusoidal profile",
        "description": "Tooth shape based on simple sinusoid.",
        "points": ["Set sharpness to 0 to give pure sinusoidal teeth.", "Shape can be tailored using sharpness which raises sinusoid to that power.", "Tooth width inactive"]
    },
    {
        "title": "Second Item",
        "image": impath+"GearSpike1.png",
        "caption": "Spike/Square",
        "header": "Spike/Square",
        "description": "Tooth shape based on straight lines. With adjustable slope. ",
        "points": ["Tooth width (%) sets width of tooth.", "Sharpness adjusts the tooth slope.", "Point Z"]
    },
    {
        "title": "3rd Item",
        "image": impath+"GearRatchet1.png",
        "caption": "Ratchet",
        "header": "Ratchet",
        "description": "Tooth shape based on spiral shape.",
        "points": ["To give straight (not curved) shape then reduce samples.","Use radial shear to tweak angles"]
    },
    {
        "title": "4th Item",
        "image": impath+"GearRounded1.png",
        "caption": "Rounded square",
        "header": "Rounded Square",
        "description": "Tooth shape based on basic square profile with the corners rounded.",
        "points": ["Tooth width (%) sets width of tooth.","Sharpness adjusts the amount of rounding."]
    },
    {
        "title": "5th Item",
        "image": impath+"GearGauss1.png",
        "caption": "Gaussian",
        "header": "Gaussian profile",
        "description": "Tooth shape based on simple gaussian.",
        "points": ["Set sharpness to 0 to give pure gaussian teeth.", "Shape can be tailored using sharpness which raises sinusoid to that power.", "Tooth width changes gaussian width."]
    },
    {
        "title": "6th Item",
        "image": impath+"GearRadShear1.png",
        "caption": "Radial shear",
        "header": "Radial shear - Spiral",
        "description": "As tooth moves away from the center it is twisted left or right.",
        "points": ["Tweak ratchets.", "Spiral effects", "Make saws..."]
    },
    {
        "title": "7th Item",
        "image": impath+"GearRadWig1.png",
        "caption": "Radial wiggle 1",
        "header": "Radial wiggle 1",
        "description": "As tooth moves away from the center it is twisted wider and narrower following a sine wave. This example uses 2 half oscillations of a sinewave to first make the tooth narrower and then wider as it moves out from the center. Midpoint reflection is on to ensure the tooth is symetrical.",
        "points": [""]
    },
    {
        "title": "8th Item",
        "image": impath+"GearRadWig2.png",
        "caption": "Radial wiggle 2",
        "header": "Radial wiggle 2",
        "description": "This example on a gaussian gear uses 4 half oscillations of a sinewave to first make the tooth wave to and fro as it moves out from the center. Midpoint reflection is OFF to ensure the teeth both sides of tooth move together.",
        "points": [""]
    }
]

# Render each item in items to give help list
for item in items:
    with st.container(border=True):
        col1, col2 = st.columns([1, 2])

        with col1:
            #st.title(item["title"])
            st.image(item["image"], caption=item["caption"])

        with col2:
            st.header(item["header"])
            st.write(item["description"])
            st.markdown("\n".join([f"- {p}" for p in item["points"]]))
