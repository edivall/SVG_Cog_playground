import streamlit as st
import svgwrite
import math
import numpy as np







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


def angle_profile(t, angle_span, angle_wave, angle_power, rfrac, rad_wig, wig_strength, reflect, radial_shear):
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
    if (t < (dist/cycles)) or (t > ((cycles-dist)/cycles)):
        mag = (abs(s)) ** angle_power
        mag =theta + angle_wave * mag * math.copysign(1, s)
    else:
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

    for i in range(params["n_points"] + 1):
        t = i / params["n_points"]
        widthFract=params["width"]/params["angle_span"]
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
            params["angle_wave"],
            params["angle_power"],
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

st.title("⚙️ SVG Cog shape Generator!")

with st.sidebar:

    st.header("Settings")
    svg_size = st.slider("SVG size", 50, 800, 300)
    st.header("Gear")
    n_teeth = st.slider("Teeth", 2, 100, 20, on_change = update_tooth_span,key="n_teeth")
    angle_span = st.slider("Tooth span (deg)", 2.0, 130.0, 18.0,  key="angle_span")
    tooth_type = st.selectbox("Profile", ["Sinusoidal", "Spike/Square","Gaussian","Rounded Square","Ratchet"])
    r_base = st.slider("Base radius", 1, 120, 72)
    tooth_height = st.slider("Tooth height", -100, 100, 10)
    power = st.slider("Sharpness (1 for pure gauss)", 0.1, 7.0, 1.0)
    width = st.slider("Tooth width (not sinusoidal)", 1.0, 99.0, 50.0)

    st.header("Radial wave")
    tooth_reflect = st.selectbox("Mid point reflection", ["Off", "On"])
    radial_shear = st.slider("Radial shear (spiral)", -90.0, 90.0, 0.0)
    radial_wiggles = st.slider("Radial wiggles (half oscillations)", -5.0, 5.0, 0.0)
    wiggle_strength = st.slider("Radial Wiggle amplitude (0=Off)", -6.0, 6.0, 2.0)
    angle_wave = st.slider("Wave strength (0=Off)", -5.0, 5.0, 0.0)
    angle_power = st.slider("Wave power", 0.5, 6.0, 2.0)


    st.header("Samples")
    n_points = st.slider("Resolution (Samples per tooth)", 1, 200, 60)

    st.header("Holes / Spokes")
    center_hole = st.slider("Center hole (0=Off)", 0, 50, 5)
    bolt_holes = st.slider("Outer holes (0=Off)", 0, 120, 6)
    bolt_radius = st.slider("Outer hole radius", 20, 150, 70)
    bolt_size = st.slider("Outer hole size", 1, 15, 3)

    spokes = st.slider("Spokes (0=Off)", 0, 12, 4)
    spoke_width = st.slider("Spoke hole width ( < 360 / num spokes)", 5, 120, 60)
    spoke_inner = st.slider("Spoke inner margin", 0, 50, 10)
    spoke_outer = st.slider("Spoke outer margin", 0, 50, 10)
    dt = st.slider("Spoke width tweak (dt)", 0.0, 20.0, 5.0)


params = {
    "svg_size": svg_size,
    "r_base": r_base,
    "tooth_height": tooth_height,
    "power": power,
    "width": width,
    "tooth_type": tooth_type,
    "tooth_reflect": tooth_reflect,
    "angle_wave": angle_wave,
    "angle_power": angle_power,
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
    "dt": dt
}


# ----------------------------
# Generate + render
# ----------------------------

polar_points = generate_tooth(params)

#angles, radius_vals, tweak_vals = build_debug_data(params)  #xy graphs
angles, radius_vals, tweak_vals = build_debug_data_from_points(polar_points)

svg_str,ssvg_str = build_svg(polar_points, params)


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
    "💾 Download Scaled SVG",
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










# Local data array for user info
items = [
    {
        "title": "First Item",
        "image": "images\GearSin1.png",
        "caption": "Sinusoidal",
        "header": "Sinusoidal profile",
        "description": "Tooth shape based on simple sinusoid.",
        "points": ["Set sharpness to 0 to give pure sinusoidal teeth.", "Shape can be tailored using sharpness which raises sinusoid to that power.", "Tooth width inactive"]
    },
    {
        "title": "Second Item",
        "image": "images\GearSpike1.png",
        "caption": "Spike/Square",
        "header": "Spike/Square",
        "description": "Tooth shape based on straight lines. With adjustable slope. ",
        "points": ["Tooth width (%) sets width of tooth.", "Sharpness adjusts the tooth slope.", "Point Z"]
    },
    {
        "title": "3rd Item",
        "image": "imagesGearRatchet1.png",
        "caption": "Ratchet",
        "header": "Ratchet",
        "description": "Tooth shape based on spiral shape.",
        "points": ["To give straight (not curved) shape then reduce samples.","Use radial shear to tweak angles"]
    },
    {
        "title": "4th Item",
        "image": "images\GearRounded1.png",
        "caption": "Rounded square",
        "header": "Rounded Square",
        "description": "Tooth shape based on basic square profile with the corners rounded.",
        "points": ["Tooth width (%) sets width of tooth.","Sharpness adjusts the amount of rounding."]
    },
    {
        "title": "5th Item",
        "image": "images\GearGauss1.png",
        "caption": "Gaussian",
        "header": "Gaussian profile",
        "description": "Tooth shape based on simple gaussian.",
        "points": ["Set sharpness to 0 to give pure gaussian teeth.", "Shape can be tailored using sharpness which raises sinusoid to that power.", "Tooth width changes gaussian width."]
    },
    {
        "title": "6th Item",
        "image": "images\GearRadShear1.png",
        "caption": "Radial shear",
        "header": "Radial shear - Spiral",
        "description": "As tooth moves away from the center it is twisted left or right.",
        "points": ["Tweak ratchets.", "Spiral effects", "Make saws..."]
    },
    {
        "title": "7th Item",
        "image": "images\GearRadWig1.png",
        "caption": "Radial wiggle 1",
        "header": "Radial wiggle 1",
        "description": "As tooth moves away from the center it is twisted wider and narrower following a sine wave. This example uses 2 half oscillations of a sinewave to first make the tooth narrower and then wider as it moves out from the center. Midpoint reflection is on to ensure the tooth is symetrical.",
        "points": [""]
    },
    {
        "title": "8th Item",
        "image": "images\GearRadWig2.png",
        "caption": "Radial wiggle 2",
        "header": "Radial wiggle 2",
        "description": "This example on a gaussian gear uses 4 half oscillations of a sinewave to first make the tooth wave to and fro as it moves out from the center. Midpoint reflection is OFF to ensure the teeth both sides of tooth move together.",
        "points": [""]
    }
]

# Render each item
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
