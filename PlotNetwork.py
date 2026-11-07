
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Circle, Ellipse, Polygon

SCALE = 10000
BUILDING = [(1300, 980), (1500, 980), (1500, 820),
            (1630, 820), (1630, 740), (1300, 740)]

name = input("Enter the network name (e.g. iter8, preliminary): ").strip()

stations = pd.read_csv(f"{name}_stations.csv")
lines = pd.read_csv(f"{name}_sightlines.csv")
ellipses = pd.read_csv(f"{name}_ellipses.csv")

for i, (design, points) in enumerate(stations.groupby("design", sort=False)):
    fig, ax = plt.subplots(figsize=(8, 8))
    coords = points.set_index("name")

    # building and tank 
    ax.add_patch(Polygon(BUILDING, fill=False, color="black", label="Building"))
    ax.add_patch(Circle((1050, 750), 120, fill=False, color="red", label="Tank")) # Circle() uses radius

    # lines of sight (pair of coordinates pairs, make a line between the two)
    for line in lines[lines.design == design].itertuples():
        # coordinate of station a
        a = coords.loc[line.station_a]
        # coordinate of station b
        b = coords.loc[line.station_b]
        ax.plot([a.x, b.x], [a.y, b.y], color="gray", lw=0.8)

    # error ellipses
    for e in ellipses[ellipses.design == design].itertuples():
        relative = e.name == "P1-P2"

        # creates an "Ellipse" object which is used to create a empty circle for each ellipse
        # note: Ellipse() wants width = 2 * desired width (major/minor diameter sort of), so half the axis values
        ax.add_patch(Ellipse(
            (e.center_x, e.center_y),
            e.major_mm * SCALE / 500,
            e.minor_mm * SCALE / 500,
            angle=e.orientation_deg,
            fill=False,
            color="red" if relative else "blue"
        ))

    # making control points black triangles
    # making non-control points blue o's
    for fixed, marker, color in [("yes", "^", "black"), ("no", "o", "blue")]:
        p = points[points.fixed == fixed]
        ax.scatter(p.x, p.y, marker=marker, color=color)

    # point names
    for p in points.itertuples():
        ax.annotate(p.name, (p.x, p.y), xytext=(5, 5),
                    textcoords="offset points")

    # plot boilerplate, attempted to match the lab graphic
    ax.set(xlim=(750, 1800), ylim=(300, 1200),
           xlabel="X (m)", ylabel="Y (m)",
           aspect="equal",
           title=f"{design} - 95% Error Ellipses (Scale: {SCALE:,}x)")

    ax.grid(True)
    ax.legend()
    fig.savefig(f"{name}_design_{i}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
