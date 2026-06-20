"""
Generate a white image with evenly spaced black circles.

Usage:
    python circle_grid.py
(edit the parameters below, or import generate_circle_grid() elsewhere)
"""

import math

from PIL import Image, ImageDraw


def generate_circle_grid(Nx, Ny, Nd, spacing=None, xmargin=None, ymargin=None,
                          pattern="square", orientation="horizontal"):
    """
    Create a white image of size (Nx, Ny) with a grid of evenly spaced
    black circles of diameter Nd.

    Args:
        Nx, Ny: image width and height in pixels
        Nd: circle diameter in pixels
        spacing: center-to-center distance between circles along a line
                 (defaults to 2*Nd, i.e. circles separated by one diameter
                 of empty space)
        xmargin, ymargin: empty border before the first column/row of circles
                (default to spacing / 2, so circles are centered in the image)
        pattern: "square" for a regular grid, or "triangular" for a
                 hexagonal/triangular packing (alternating lines offset by
                 half a spacing, lines packed closer together)
        orientation: "horizontal" (default) packs rows tightly together
                 vertically, with circles offset along x from row to row.
                 "vertical" rotates this 90 degrees: columns are packed
                 tightly together horizontally, with circles offset along y
                 from column to column. Only affects "triangular" pattern.

    Returns:
        PIL.Image
    """
    if spacing is None:
        spacing = 2 * Nd
    if xmargin is None:
        xmargin = spacing / 2
    if ymargin is None:
        ymargin = spacing / 2

    # "Tight" spacing is the distance between successive lines (rows for
    # horizontal orientation, columns for vertical orientation). It's
    # reduced for a triangular pattern so the packing is equilateral.
    # "Loose" spacing is the distance between circles along a single line.
    tight_spacing = spacing * (math.sqrt(3) / 2) if pattern == "triangular" else spacing
    loose_spacing = spacing

    img = Image.new("L", (Nx, Ny), color=255)  # white background
    draw = ImageDraw.Draw(img)
    r = Nd / 2

    if orientation == "vertical":
        # Lines run vertically (columns); tight axis is x, loose axis is y.
        u_max, u_margin = Nx, xmargin
        v_max, v_margin = Ny, ymargin
    else:
        # Lines run horizontally (rows); tight axis is y, loose axis is x.
        u_max, u_margin = Ny, ymargin
        v_max, v_margin = Nx, xmargin

    u = u_margin
    line = 0
    while u <= (u_max - u_margin):
        v_offset = (loose_spacing / 2) if (pattern == "triangular" and line % 2 == 1) else 0
        v = v_margin + v_offset
        while v <= (v_max - v_margin):
            if orientation == "vertical":
                x, y = u, v
            else:
                x, y = v, u
            draw.ellipse([x - r, y - r, x + r, y + r], fill=0)  # black circle
            v += loose_spacing
        u += tight_spacing
        line += 1

    return img


if __name__ == "__main__":
    Nx, Ny = 200, 90   # image size
    Nd = 20              # circle diameter

    img = generate_circle_grid(Nx, Ny, Nd, pattern="triangular", orientation="vertical",
                                xmargin=30, ymargin=5, spacing=Nd * 2)
    img.save("circle_grid.png")
    print("Saved circle_grid.png")