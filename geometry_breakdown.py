import math
import tkinter as tk
from tkinter import messagebox, filedialog
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import screeninfo
from screeninfo import get_monitors
import numpy as np
from matplotlib.patches import Polygon
import os
import datetime
from PIL import Image, ImageDraw, ImageFont
from io import BytesIO
from scipy.integrate import quad

# Golden ratio constant
φ = (1 + math.sqrt(5)) / 2

# Color palette for shapes
shape_colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:gray"]

# Global Button Styling
BUTTON_NORMAL = "#f4d25d"
BUTTON_HOVER = "#d1b139"
BUTTON_FONT = ("Segoe UI", 12, "bold")

def styled_button(master, text, command, **kwargs):
    btn = tk.Button(master, text=text, command=command,
                    bg=BUTTON_NORMAL,
                    activebackground=BUTTON_HOVER,
                    font=BUTTON_FONT,
                    relief="raised", bd=2,
                    **kwargs)  # allow width, height, etc.
    btn.bind("<Enter>", lambda e: btn.config(bg=BUTTON_HOVER))
    btn.bind("<Leave>", lambda e: btn.config(bg=BUTTON_NORMAL))
    return btn

# Global Font Styling
def apply_global_styles(root):
    default_font = ("Arial", 14)
    default_fg = "#0a2e4f"       # Deep blue from your palette
    default_bg = "#cddbe8"       # Light blue from your palette

    # Global font and text color
    root.option_add("*Font", default_font)
    root.option_add("*Foreground", default_fg)

    # Optional background settings
    root.option_add("*Label.Background", default_bg)
    root.option_add("*Frame.Background", default_bg)
    root.option_add("*Checkbutton.Background", default_bg)
    root.option_add("*Menubutton.Background", default_bg)
    root.option_add("*Menu.Background", default_bg)

# Base class for all shapes
class Shape:
    def summary(self):
        raise NotImplementedError("Each shape must implement its own summary method.")

class Circle:
    def __init__(self, radius):
        self.radius = radius

    def summary(self):
        return {
            "Shape": "Circle",
            "Radius": self.radius,
            "Area": math.pi * self.radius ** 2,
            "Circumference": 2 * math.pi * self.radius,
            "Diameter": 2 * self.radius,
            "Is Golden Ratio Approx": math.isclose(self.radius * 2, φ, rel_tol=0.1)
        }

    def visualize(self, ax, color="black", label=None, annotate=False):
        circle = plt.Circle((0, 0), self.radius, fill=False, color=color, label=label)
        ax.add_patch(circle)

        # Annotations
        if annotate:
            ax.plot(0, 0, marker="o", color=color)
            ax.text(0, 0, "Center", ha="right", va="bottom", fontsize=8)

            ax.plot([0, self.radius], [0, 0], color=color, linestyle="--")
            ax.text(self.radius / 2, 0, f"r={self.radius:.1f}", ha="center", va="top", fontsize=8)

        ax.set_xlim(-self.radius * 1.5, self.radius * 1.5)
        ax.set_ylim(-self.radius * 1.5, self.radius * 1.5)
        ax.set_aspect("equal")
        ax.grid(True)

# Oval (Ellipse)
class Ellipse(Shape):
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def area(self):
        return math.pi * self.a * self.b

    def circumference(self):
        h = ((self.a - self.b) ** 2) / ((self.a + self.b) ** 2)
        return math.pi * (self.a + self.b) * (1 + (3 * h) / (10 + math.sqrt(4 - 3 * h)))

    def summary(self):
        return {
            "Shape": "Ellipse",
            "Major Axis (a)": self.a,
            "Minor Axis (b)": self.b,
            "Area": self.area(),
            "Circumference (Approx)": self.circumference(),
            "Eccentricity": math.sqrt(1 - (self.b ** 2) / (self.a ** 2)),
            "Is Golden Ratio Approx": math.isclose(self.a / self.b, φ, rel_tol=0.05)
        }

    def visualize(self, ax, color="black", label=None, annotate=False):
        t = np.linspace(0, 2 * math.pi, 500)
        x = self.a * np.cos(t)
        y = self.b * np.sin(t)
        ax.plot(x, y, color=color, label=label)

        if annotate:
            ax.plot(0, 0, marker="o", color=color)
            ax.text(0, 0, "Center", ha="right", va="bottom", fontsize=8)

            ax.plot([0, self.a], [0, 0], color=color, linestyle="--")
            ax.text(self.a / 2, 0, f"a={self.a:.1f}", ha="center", va="top", fontsize=8)

            ax.plot([0, 0], [0, self.b], color=color, linestyle=":")
            ax.text(0, self.b / 2, f"b={self.b:.1f}", ha="left", va="center", fontsize=8)

        ax.set_xlim(-self.a * 1.5, self.a * 1.5)
        ax.set_ylim(-self.b * 1.5, self.b * 1.5)
        ax.set_aspect("equal")
        ax.grid(True)

class Triangle:
    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c

    def summary(self):
        s = (self.a + self.b + self.c) / 2
        area = math.sqrt(s * (s - self.a) * (s - self.b) * (s - self.c))
        angles = self._angles_deg()
        return {
            "Shape": "Triangle",
            "Sides": (self.a, self.b, self.c),
            "Area": area,
            "Perimeter": self.a + self.b + self.c,
            "Angles (deg)": angles,
            "Is Equilateral": self._is_equilateral(),
            "Is Isosceles": self._is_isosceles(),
            "Is Scalene": self._is_scalene(),
            "Is Right Triangle": self._is_right(),
            "Is Golden Ratio Approx": math.isclose(max(self.a, self.b, self.c) / min(self.a, self.b, self.c), φ, rel_tol=0.1)
        }

    def _angles_deg(self):
        a, b, c = self.a, self.b, self.c
        A = math.degrees(math.acos((b**2 + c**2 - a**2) / (2*b*c)))
        B = math.degrees(math.acos((a**2 + c**2 - b**2) / (2*a*c)))
        C = 180 - A - B
        return (round(A, 1), round(B, 1), round(C, 1))

    def _is_equilateral(self):
        return math.isclose(self.a, self.b) and math.isclose(self.b, self.c)

    def _is_isosceles(self):
        return self.a == self.b or self.b == self.c or self.a == self.c

    def _is_scalene(self):
        return self.a != self.b and self.b != self.c and self.a != self.c

    def _is_right(self):
        sides = sorted([self.a, self.b, self.c])
        return math.isclose(sides[0]**2 + sides[1]**2, sides[2]**2)

    def visualize(self, ax, color="black", label=None, annotate=False):
        A = (0, 0)
        B = (self.c, 0)
        cos_angle = (self.a**2 + self.c**2 - self.b**2) / (2 * self.a * self.c)
        sin_angle = math.sqrt(1 - cos_angle**2)
        C = (self.a * cos_angle, self.a * sin_angle)

        triangle = Polygon([A, B, C], closed=True, fill=False, edgecolor=color, label=label)
        ax.add_patch(triangle)

        if annotate:
            def midpoint(p1, p2):
                return ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)

            def dist(p1, p2):
                return math.hypot(p1[0] - p2[0], p1[1] - p2[1])

            ab = dist(A, B)
            bc = dist(B, C)
            ca = dist(C, A)

            # Match actual side lengths to input values
            side_map = {
                ab: ("c", midpoint(A, B)),
                bc: ("a", midpoint(B, C)),
                ca: ("b", midpoint(C, A))
            }

            input_sides = {
                round(self.a, 5): "a",
                round(self.b, 5): "b",
                round(self.c, 5): "c"
            }

            for d, (label, mid) in side_map.items():
                d_rounded = round(d, 5)
                for input_val, input_label in input_sides.items():
                    if math.isclose(d_rounded, input_val, rel_tol=1e-3):
                        ax.text(*mid, f"{input_label}={input_val:.1f}", fontsize=8, ha="center", va="center")
                        break

        xs = [A[0], B[0], C[0]]
        ys = [A[1], B[1], C[1]]
        ax.set_xlim(min(xs)-1, max(xs)+1)
        ax.set_ylim(min(ys)-1, max(ys)+1)
        ax.set_aspect("equal")
        ax.grid(True)

class Rectangle:
    def __init__(self, width, height):
        self.width = width
        self.height = height

    def summary(self):
        return {
            "Shape": "Rectangle",
            "Width": self.width,
            "Height": self.height,
            "Area": self.width * self.height,
            "Perimeter": 2 * (self.width + self.height),
            "Is Golden Ratio Approx": math.isclose(max(self.width, self.height) / min(self.width, self.height), φ, rel_tol=0.1)
        }

    def visualize(self, ax, color="black", label=None, annotate=False):
        corners = [(0, 0), (self.width, 0), (self.width, self.height), (0, self.height)]
        rect = Polygon(corners, closed=True, fill=False, edgecolor=color, label=label)
        ax.add_patch(rect)

        if annotate:
            def midpoint(p1, p2):
                return ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)

            ax.text(*midpoint(corners[0], corners[1]), f"Width = {self.width:.1f}", fontsize=8, ha="center", va="top")
            ax.text(*midpoint(corners[1], corners[2]), f"Height = {self.height:.1f}", fontsize=8, ha="left", va="center")

        ax.set_xlim(-1, self.width + 1)
        ax.set_ylim(-1, self.height + 1)
        ax.set_aspect("equal")
        ax.grid(True)

class RegularPolygon:
    def __init__(self, n_sides, radius):
        self.n = n_sides
        self.radius = radius

    def summary(self):
        side_length = 2 * self.radius * math.sin(math.pi / self.n)
        perimeter = self.n * side_length
        area = 0.5 * self.n * self.radius ** 2 * math.sin(2 * math.pi / self.n)
        return {
            "Shape": f"Regular {self.n}-gon",
            "Sides": self.n,
            "Radius": self.radius,
            "Side Length": side_length,
            "Perimeter": perimeter,
            "Area": area,
            "Is Golden Ratio Approx": math.isclose(perimeter / self.radius, φ, rel_tol=0.1)
        }

    def visualize(self, ax, color="black", label=None, annotate=False):
        angles = [2 * math.pi * i / self.n for i in range(self.n)]
        points = [(self.radius * math.cos(a), self.radius * math.sin(a)) for a in angles]
        poly = Polygon(points, closed=True, fill=False, edgecolor=color, label=label)
        ax.add_patch(poly)

        if annotate:
            for i, (x, y) in enumerate(points):
                ax.plot(x, y, marker="o", color=color)
                ax.text(x, y, f"V{i+1}", fontsize=8, ha="right", va="bottom")

            for i in range(self.n):
                p1 = points[i]
                p2 = points[(i + 1) % self.n]
                mid = ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)
                side_length = math.hypot(p2[0] - p1[0], p2[1] - p1[1])
                ax.text(*mid, f"s={side_length:.1f}", fontsize=8, ha="center", va="center")

            ax.plot(0, 0, marker="x", color=color)
            ax.text(0, 0, "Center", fontsize=8, ha="left", va="top")

        ax.set_xlim(-self.radius * 1.2, self.radius * 1.2)
        ax.set_ylim(-self.radius * 1.2, self.radius * 1.2)
        ax.set_aspect("equal")
        ax.grid(True)

class IrregularPolygon:
    def __init__(self, points):
        self.points = points

    def summary(self):
        from shapely.geometry import Polygon as ShapelyPolygon
        from shapely.ops import polygonize, unary_union
        from shapely.geometry import LineString

        # Create a shapely polygon from the input points
        ring = LineString(self.points + [self.points[0]])
        merged = unary_union(ring)

        # Break into simple polygons
        polygons = list(polygonize(merged))

        # Compute total area as the sum of simple parts
        total_area = sum(abs(p.area) for p in polygons)

        def perimeter(p):
            return sum(math.hypot(p[i][0] - p[(i+1)%len(p)][0], p[i][1] - p[(i+1)%len(p)][1]) for i in range(len(p)))

        return {
            "Shape": "Irregular Polygon",
            "Vertices": len(self.points),
            "Perimeter": perimeter(self.points),
            "Area": total_area
        }

    def visualize(self, ax, color="black", label=None, annotate=False):
        poly = Polygon(self.points, closed=True, fill=False, edgecolor=color, label=label)
        ax.add_patch(poly)

        if annotate:
            for i, (x, y) in enumerate(self.points):
                ax.plot(x, y, marker="o", color=color)
                ax.text(x, y, f"({x:.1f},{y:.1f})", fontsize=8, ha="right", va="bottom")

            for i in range(len(self.points)):
                p1 = self.points[i]
                p2 = self.points[(i + 1) % len(self.points)]
                mid = ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)
                side_length = math.hypot(p2[0] - p1[0], p2[1] - p1[1])
                ax.text(*mid, f"s={side_length:.1f}", fontsize=8, ha="center", va="center")

        all_x = [x for x, _ in self.points]
        all_y = [y for _, y in self.points]
        ax.set_xlim(min(all_x) - 1, max(all_x) + 1)
        ax.set_ylim(min(all_y) - 1, max(all_y) + 1)
        ax.set_aspect("equal")
        ax.grid(True)

class Spiral:
    def __init__(self, spiral_type, a, b, theta):
        self.spiral_type = spiral_type.lower()
        self.a = a
        self.b = b
        self.theta = theta

    def summary(self):
        try:
            if self.spiral_type == "archimedean":
                r_final = self.a + self.b * self.theta
                arc_length, _ = quad(lambda t: np.sqrt((self.a + self.b * t) ** 2 + self.b ** 2), 0, self.theta)

            elif self.spiral_type == "logarithmic":
                r_final = self.a * math.exp(self.b * self.theta)
                arc_length = (self.a * math.sqrt(1 + self.b ** 2) / self.b) * (math.exp(self.b * self.theta) - 1)

            elif self.spiral_type == "fermat":
                r_final = math.sqrt(self.a**2 + self.b**2 * self.theta**2)
                arc_length, _ = quad(lambda t: math.sqrt((self.b ** 2 * t / math.sqrt(self.a ** 2 + self.b ** 2 * t ** 2)) ** 2 +
                                                         (math.sqrt(self.a ** 2 + self.b ** 2 * t ** 2)) ** 2), 0, self.theta)

            elif self.spiral_type == "hyperbolic":
                r_final = self.a / (1 + self.b * self.theta)
                arc_length, _ = quad(lambda t: math.sqrt((self.a / (1 + self.b * t)) ** 2 +
                                                         (self.a * self.b / (1 + self.b * t) ** 2) ** 2), 0, self.theta)
            else:
                raise ValueError("Unsupported spiral type")

            end_x = r_final * math.cos(self.theta)
            end_y = r_final * math.sin(self.theta)

            return {
                "Shape": "Spiral",
                "Type": self.spiral_type.title(),
                "a": self.a,
                "b": self.b,
                "Theta (rad)": self.theta,
                "Final Radius": round(r_final, 5),
                "End Coordinates": (round(end_x, 5), round(end_y, 5)),
                "Arc Length": round(arc_length, 5)
            }

        except Exception as e:
            return {
                "Shape": "Spiral",
                "Type": self.spiral_type.title(),
                "a": self.a,
                "b": self.b,
                "Theta (rad)": self.theta,
                "Error": str(e)
            }

    def visualize(self, ax, color="black", label=None, annotate=False):
        t = np.linspace(0, self.theta, 1000)

        if self.spiral_type == "archimedean":
            r = self.a + self.b * t
        elif self.spiral_type == "logarithmic":
            r = self.a * np.exp(self.b * t)
        elif self.spiral_type == "fermat":
            r = np.sqrt(self.a**2 + self.b**2 * t**2)
        elif self.spiral_type == "hyperbolic":
            r = self.a / (1 + self.b * t)
        else:
            raise ValueError("Unsupported spiral type")

        x = r * np.cos(t)
        y = r * np.sin(t)

        ax.plot(x, y, color=color, label=label)

        if annotate:
            ax.plot(0, 0, marker="o", color=color)
            ax.text(0, 0, "Origin", fontsize=8, ha="right", va="bottom")

            # Final radius and angle label
            final_x = x[-1]
            final_y = y[-1]
            ax.plot(final_x, final_y, marker="x", color=color)
            ax.text(final_x, final_y, f"θ={self.theta:.2f} rad", fontsize=8, ha="left", va="bottom")

        ax.set_xlim(min(x) - 1, max(x) + 1)
        ax.set_ylim(min(y) - 1, max(y) + 1)
        ax.set_aspect("equal")
        ax.grid(True)

# Display summary
def display_summary(summary):
    output = "\nShape Summary:\n"
    for k, v in summary.items():
        output += f"{k}: {v}\n"
    return output

# Window creation and layout handler
class ShapeAnalyzerUI(tk.Frame):
    def __init__(self, parent, shape_id=None, on_metrics_updated=None):
        super().__init__(parent)
        self.shape_id = shape_id
        self.on_metrics_updated = on_metrics_updated
        self.parameters = {}
        self.metrics = {}
        self.entries = {}
        self.build_ui()

    def build_ui(self):
        main_frame = tk.Frame(self, bg="#f0f4f8")
        main_frame.pack(fill="both", expand=True)

        control_frame = tk.Frame(main_frame, bg="#f0f4f8")
        control_frame.pack(side="left", fill="y", padx=10, pady=10)

        output_frame = tk.Frame(main_frame, bg="#f0f4f8")
        output_frame.pack(side="left", expand=True, fill="both", padx=10, pady=10)

        self.result_label = tk.Label(output_frame, text="Shape Data", font=("Courier", 10), anchor="nw", justify="left", bg="#f0f4f8")
        self.result_label.pack(anchor="nw", fill="both", expand=True)

        canvas_frame = tk.Frame(output_frame, bg="#f0f4f8")
        canvas_frame.pack(fill="both", expand=True)

        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=canvas_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        self.shape_var = tk.StringVar(value="Circle")
        self.param_frame = tk.Frame(control_frame, bg="#f0f4f8")
        self.param_frame.pack(pady=10, fill="x")

        self.options = [
            ("Circle", ["radius"]),
            ("Ellipse", ["a", "b"]),
            ("Triangle", ["a", "b", "c"]),
            ("Rectangle", ["width", "height"]),
            ("Regular Polygon", ["sides", "length"]),
            ("Irregular Polygon", ["vertices"]),
            ("Spiral", ["a", "b", "theta", "type"]),
        ]

        def update_params(*args):
            for name, params in self.options:
                if name == self.shape_var.get():
                    self.prompt_params(params)
                    break

        tk.Label(control_frame, text="Select a shape:", bg="#f0f4f8").pack()
        shape_dropdown = ttk.Combobox(control_frame, textvariable=self.shape_var, values=[name for name, _ in self.options])
        shape_dropdown.pack()
        shape_dropdown.bind("<<ComboboxSelected>>", update_params)
        update_params()

        styled_button(control_frame, text="Go!", command=self.analyze_shape).pack(pady=10)

    def prompt_params(self, params):
        for widget in self.param_frame.winfo_children():
            widget.destroy()
        self.entries.clear()

        format_note = {
            "vertices": "Use: (x1,y1); (x2,y2); (x3,y3); ...\nOR: x1,y1 x2,y2 x3,y3 ...\nOR: x1 y1 x2 y2 x3 y3 ...",
        }

        for param in params:
            label_text = {
                "radius": "Radius:",
                "width": "Width:",
                "height": "Height:",
                "a": "a:",
                "b": "b:",
                "c": "c:",
                "sides": "Number of Sides:",
                "length": "Side Length:",
                "theta": "Max Angle (θ, radians):",
                "type": "Spiral Type:",
                "vertices": "Vertices:"
            }.get(param, param + ":")

            lbl = tk.Label(self.param_frame, text=label_text, bg="#f0f4f8")
            lbl.pack()

            if param in format_note:
                tip = tk.Label(self.param_frame, text=format_note[param], font=("Arial", 12), fg="gray", bg="#f0f4f8")
                tip.pack()

            if self.shape_var.get() == "Spiral" and param == "a":
                note = tk.Label(self.param_frame, text="Initial Radius", font=("Arial", 12), fg="black", bg="#f0f4f8")
                note.pack()

            if self.shape_var.get() == "Spiral" and param == "b":
                note = tk.Label(self.param_frame, text="Spiral Constant", font=("Arial", 12), fg="black", bg="#f0f4f8")
                note.pack()

            if self.shape_var.get() == "Ellipse":
                note = tk.Label(self.param_frame, text="Enter values for semi-minor and semi-major axis. Use a > b for horizontal, b > a for vertical", font=("Arial", 12), fg="gray", bg="#f0f4f8")
                note.pack()

            if self.shape_var.get() == "Triangle":
                note = tk.Label(self.param_frame, text="Enter length of side", font=("Arial", 12), fg="black", bg="#f0f4f8")
                note.pack()

            if param == "type" and self.shape_var.get() == "Spiral":
                combo = ttk.Combobox(self.param_frame, values=["archimedean", "logarithmic", "fermat", "hyperbolic"])
                combo.set("archimedean")
                combo.pack()
                self.entries[param] = combo
            else:
                ent = tk.Entry(self.param_frame)
                ent.pack()
                self.entries[param] = ent

    def get_entry_values(self):
        try:
            result = {}
            for k, v in self.entries.items():
                if k == "type":
                    result[k] = v.get().lower()
                elif k == "vertices":
                    raw = v.get().strip()

                    # Normalize: Replace all parentheses, multiple spaces, and semicolons with commas
                    clean = raw.replace("(", "").replace(")", "").replace(";", " ")
                    tokens = clean.replace(",", " ").split()

                    if len(tokens) % 2 != 0:
                        raise ValueError("Irregular Polygon input must contain x,y pairs.")

                    vertices = []
                    for i in range(0, len(tokens), 2):
                        x = float(tokens[i])
                        y = float(tokens[i+1])
                        vertices.append((x, y))

                    result[k] = vertices

                else:
                    result[k] = float(v.get())
            return result
        except Exception as e:
            messagebox.showerror("Input Error", str(e))
            return None

    def analyze_shape(self):
        try:
            shape_type = self.shape_var.get()
            values = self.get_entry_values()
            if values is None:
                return

            if shape_type == "Circle":
                shape = Circle(values["radius"])
            elif shape_type == "Ellipse":
                shape = Ellipse(values["a"], values["b"])
            elif shape_type == "Triangle":
                shape = Triangle(values["a"], values["b"], values["c"])
            elif shape_type == "Rectangle":
                shape = Rectangle(values["width"], values["height"])
            elif shape_type == "Regular Polygon":
                shape = RegularPolygon(int(values["sides"]), values["length"])
            elif shape_type == "Irregular Polygon":
                shape = IrregularPolygon(values["vertices"])
            elif shape_type == "Spiral":
                shape = Spiral(values["type"], values["a"], values["b"], values["theta"])
            else:
                raise ValueError("Unknown shape")

            self.shape = shape
            color = shape_colors[(self.shape_id - 1) % len(shape_colors)]
            self.ax.clear()
            shape.visualize(self.ax, color=color, label=f"Shape {self.shape_id}", annotate=True)
            self.ax.legend()
            self.canvas.draw()
            self.result_label.config(text=display_summary(shape.summary()))

            if self.on_metrics_updated:
                self.on_metrics_updated(self.parameters, shape.summary())

        except Exception as e:
            messagebox.showerror("Analysis Error", str(e))

# Shape window with tiling and compare support
open_windows = []

class ShapeWindow:
    def __init__(self, parent):
        self.shape_id = len(open_windows) + 1

        self.frame = tk.Frame(parent, bd=2, relief=tk.RIDGE, bg="#f0f4f8")
        self.frame.pack(side=tk.LEFT, fill=tk.Y, padx=2, pady=2)

        self.header = tk.Frame(self.frame, bg="#f0f4f8")
        self.header.pack(fill=tk.X)

        self.label = tk.Label(self.header, text=f"Shape {self.shape_id}", font=("Arial", 12, "bold"), bg="#f0f4f8")
        self.label.pack(side=tk.LEFT, padx=5)

        close_btn = styled_button(self.header, text="✕", command=self.close_window, width=2)
        close_btn.pack(side=tk.RIGHT)

        self.analyzer = ShapeAnalyzerUI(self.frame, shape_id=self.shape_id, on_metrics_updated=self.set_metrics)
        self.analyzer.pack(fill=tk.BOTH, expand=True)

        self.controls = tk.Frame(self.frame, bg="#f0f4f8")
        self.controls.pack(fill=tk.X)

        self.compare_var = tk.IntVar()
        tk.Checkbutton(self.controls, text="Include in Comparison", variable=self.compare_var, command=self.update_comparison).pack(side=tk.LEFT, padx=4)

        save_btn = styled_button(self.controls, text="Save Graph", command=self.save_graph)
        save_btn.pack(side=tk.RIGHT, padx=4)

        summary_btn = styled_button(self.controls, text="Save Graph With Summary", command=self.save_graph_with_summary)
        summary_btn.pack(side=tk.RIGHT, padx=4)

        self.shape_summary = None

        open_windows.append(self)
        self.renumber_all_windows()
        self.resize_windows()

    def spawn_new_window(self):
        app.open_new_shape_window()

    def close_window(self):
        if self in open_windows:
            open_windows.remove(self)
        self.frame.destroy()
        self.renumber_all_windows()
        for win in open_windows:
            win.resize_windows()

    def renumber_all_windows(self):
        for idx, win in enumerate(open_windows, start=1):
            win.shape_id = idx
            win.label.config(text=f"Shape {idx}")

    def set_metrics(self, params, metrics):
        self.parameters = params
        self.metrics = metrics
        if hasattr(self.analyzer, 'shape') and hasattr(self.analyzer.shape, 'summary'):
            self.shape_summary = self.analyzer.shape.summary()

    def update_comparison(self):
        global comparison_window_instance
        selected = [w for w in open_windows if w.compare_var.get() == 1]
        if len(selected) >= 2:
            if comparison_window_instance is None or not comparison_window_instance.window.winfo_exists():
                comparison_window_instance = ComparisonWindow(selected)
            else:
                comparison_window_instance.refresh(selected)

    def resize_windows(self):
        screen = get_monitors()[0]
        num = len(open_windows)
        if num == 0:
            return
        width = screen.width // min(num, 4)
        max_height = int(screen.height * 0.9)
        height = max_height // ((num - 1) // 4 + 1)
        for win in open_windows:
            win.frame.config(width=width, height=height)
    def save_graph(self):
        if not hasattr(self.analyzer, 'shape'):
            messagebox.showerror("Error", "Please analyze a shape first.")
            return
        try:
            file_path = filedialog.asksaveasfilename(
                defaultextension=".png",
                filetypes=[("PNG Image", "*.png")],
                title="Save Shape Graph",
                initialfile=f"shape_{self.shape_id}.png"
            )
            if not file_path:
                return
            self.analyzer.fig.savefig(file_path, dpi=300, bbox_inches="tight")
            messagebox.showinfo("Saved", f"Graph saved as {file_path}")
        except Exception as e:
            messagebox.showerror("Save Error", str(e))

    def save_graph_with_summary(self):
        if not hasattr(self.analyzer, 'shape'):
            messagebox.showerror("Error", "Please analyze a shape first.")
            return
        try:

            # Save graph to memory
            buf = BytesIO()
            self.analyzer.fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
            buf.seek(0)
            graph_img = PIL.Image.open(buf).convert("RGB")

            # Prepare summary text
            summary = display_summary(self.analyzer.shape.summary())
            lines = summary.strip().split("\n")

            # Load font and create dummy image to measure text
            try:
                font = PIL.ImageFont.truetype("arial.ttf", size=32)
            except IOError:
                font = PIL.ImageFont.load_default()
            dummy_img = PIL.Image.new("RGB", (10, 10))
            draw = PIL.ImageDraw.Draw(dummy_img)

            # Measure text size
            padding = 10
            line_height = draw.textbbox((0, 0), "A", font=font)[3] + 2
            text_width = max(draw.textbbox((0, 0), line, font=font)[2] for line in lines) + 2 * padding
            text_height = len(lines) * line_height + 2 * padding

            # Create image for text
            text_img = PIL.Image.new("RGB", (text_width, text_height), "white")
            draw = PIL.ImageDraw.Draw(text_img)
            for i, line in enumerate(lines):
                draw.text((padding, padding + i * line_height), line, fill="black", font=font)

            # Combine graph + text vertically
            total_width = max(graph_img.width, text_img.width)
            total_height = graph_img.height + text_img.height
            combined = PIL.Image.new("RGB", (total_width, total_height), "white")
            combined.paste(graph_img, (0, 0))
            combined.paste(text_img, (0, graph_img.height))

            # Save combined image
            file_path = filedialog.asksaveasfilename(
                defaultextension=".png",
                filetypes=[("PNG Image", "*.png")],
                title="Save Shape Summary",
                initialfile=f"shape_{self.shape_id}_summary.png"
            )
            if not file_path:
                return
            combined.save(file_path)
            messagebox.showinfo("Saved", f"Graph + Summary saved as {file_path}")
        except Exception as e:
            messagebox.showerror("Save Error", str(e))

comparison_window_instance = None  # Global reference to the one comparison window

class ComparisonWindow:
    def __init__(self, shape_windows):
        global comparison_window_instance
        comparison_window_instance = self

        self.shape_windows = shape_windows
        self.window = tk.Toplevel()
        self.window.configure(bg="#f0f4f8")
        self.window.title("Shape Comparison")
        self.window.protocol("WM_DELETE_WINDOW", self.close_window)

        # Frame + Scrollbars
        frame = tk.Frame(self.window, bg="#f0f4f8")
        frame.pack(fill=tk.BOTH, expand=True)

        self.y_scrollbar = tk.Scrollbar(frame, orient=tk.VERTICAL)
        self.y_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        self.x_scrollbar = tk.Scrollbar(frame, orient=tk.HORIZONTAL)
        self.x_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)

        self.table = tk.Text(frame, wrap=tk.NONE,
                             yscrollcommand=self.y_scrollbar.set,
                             xscrollcommand=self.x_scrollbar.set)
        self.table.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.y_scrollbar.config(command=self.table.yview)
        self.x_scrollbar.config(command=self.table.xview)

        save_btn = styled_button(self.window, text="Save Comparison", command=self.save_comparison)
        save_btn.pack(pady=5)

        self.refresh(shape_windows)

    def close_window(self):
        global comparison_window_instance
        comparison_window_instance = None
        self.window.destroy()

    def refresh(self, shape_windows):
        self.shape_windows = shape_windows
        self.table.delete("1.0", tk.END)

        for i in range(len(shape_windows)):
            for j in range(i + 1, len(shape_windows)):
                label1 = shape_windows[i].label.cget("text")
                label2 = shape_windows[j].label.cget("text")
                self.table.insert(tk.END, f"Comparison: {label1} vs {label2}\n")

                summary1 = shape_windows[i].shape_summary
                summary2 = shape_windows[j].shape_summary

                if not summary1 or not summary2:
                    self.table.insert(tk.END, "  ⚠ One or both shapes have not been analyzed yet.\n\n")
                    continue

                keys = set(summary1.keys()).union(summary2.keys())
                for key in sorted(keys):
                    val1 = summary1.get(key, "N/A")
                    val2 = summary2.get(key, "N/A")

                    key_aliases = {
                        ("Perimeter", "Circumference"),
                        ("Circumference", "Perimeter"),
                    }

                    if val1 == val2:
                        comparison = "Equal"
                    else:
                        try:
                            if key not in summary2 and any(key == k1 and k2 in summary2 for (k1, k2) in key_aliases):
                                val2 = summary2.get("Circumference", summary2.get("Perimeter", "N/A"))
                            if key not in summary1 and any(key == k2 and k1 in summary1 for (k1, k2) in key_aliases):
                                val1 = summary1.get("Circumference", summary1.get("Perimeter", "N/A"))

                            diff = abs(float(val1) - float(val2))
                            comparison = f"Difference: {diff:.4f}"
                        except:
                            comparison = "Different"

                    self.table.insert(tk.END, f"  {key}: {val1} vs {val2} → {comparison}\n")

                self.table.insert(tk.END, "\n")

    def save_comparison(self):
        try:
            from tkinter import filedialog
            file_path = filedialog.asksaveasfilename(
                defaultextension=".txt",
                filetypes=[("Text Files", "*.txt")],
                title="Save Comparison Summary"
            )
            if not file_path:
                return

            content = self.table.get("1.0", tk.END)
            with open(file_path, "w", encoding="utf-8") as f:
                f.write(content)

            messagebox.showinfo("Saved", f"Comparison saved to {file_path}")
        except Exception as e:
            messagebox.showerror("Save Error", str(e))

# App
class GeometryApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Geometry Breakdown")
        self.root.state('zoomed')  # Fullscreen by default on Windows
        try:
            root.state('zoomed')  # Windows
        except:
            root.attributes('-zoomed', True)  # Linux/Mac fallback
        self.root.protocol("WM_DELETE_WINDOW", self.quit_app)
        self.root.configure(bg="#cddbe8")

        self.canvas = tk.Canvas(root, bg="#cddbe8")
        self.scroll_x = tk.Scrollbar(root, orient=tk.HORIZONTAL, command=self.canvas.xview, width=30)
        self.canvas.configure(xscrollcommand=self.scroll_x.set)

        self.scroll_x.pack(fill=tk.X, side=tk.BOTTOM)
        self.canvas.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)

        self.shape_container = tk.Frame(self.canvas, bg="#f0f4f8")
        self.canvas.create_window((0, 0), window=self.shape_container, anchor='nw')

        # Main toolbar
        self.toolbar = tk.Frame(self.root, bg="#bccee0", pady=5)
        self.toolbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Add New Shape
        add_shape_btn = styled_button(self.toolbar, text="Add New Shape", command=self.open_new_shape_window)
        add_shape_btn.pack(pady=15, fill="x")

        # Compare Selected
        compare_btn = styled_button(self.toolbar, text="Compare Selected", command=self.compare_selected_shapes)
        compare_btn.pack(pady=15, fill="x")

        # Batch Export dropdown menu
        export_btn = tk.Menubutton(
            self.toolbar, text="Batch Export...",
            relief="raised", bd=2, font=BUTTON_FONT,
            bg=BUTTON_NORMAL, activebackground=BUTTON_HOVER
        )
        export_menu = tk.Menu(export_btn, tearoff=0)
        export_btn.config(menu=export_menu)

        #Hover style
        export_btn.bind("<Enter>", lambda e: export_btn.config(bg=BUTTON_HOVER))
        export_btn.bind("<Leave>", lambda e: export_btn.config(bg=BUTTON_NORMAL))

        export_menu.add_command(label="Export All", command=self.export_all_shapes)  # Placeholder
        export_menu.add_command(label="Export Selected", command=self.export_selected_shapes)  # Placeholder
        export_btn.pack(pady=15, fill="x")

        self.shape_container.bind("<Configure>", self.on_frame_configure)

        self.open_new_shape_window()

        self.root.update_idletasks()  # Force layout calculations
        self.root.state('zoomed')     # Maximize window after layout is stable

    def on_frame_configure(self, event):
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))

    def open_new_shape_window(self):
        ShapeWindow(self.shape_container)

    def quit_app(self):
        self.root.quit()
        self.root.destroy()

    def export_all_shapes(self):
        from PIL import Image, ImageDraw, ImageFont
        from io import BytesIO

        if not open_windows:
            messagebox.showinfo("No Shapes", "There are no shapes to export.")
            return

        base_folder = filedialog.askdirectory(title="Select Folder to Save All Summaries")
        if not base_folder:
            return

        timestamp = datetime.datetime.now().strftime("export_%Y-%m-%d_%H%M%S")
        folder = os.path.join(base_folder, timestamp)
        os.makedirs(folder, exist_ok=True)

        if not folder:
            return

        for window in open_windows:
            if hasattr(window.analyzer, 'shape'):
                try:
                    # Save graph + summary to memory

                    buf = BytesIO()
                    window.analyzer.fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
                    buf.seek(0)
                    graph_img = PIL.Image.open(buf).convert("RGB")

                    summary = display_summary(window.analyzer.shape.summary())
                    lines = summary.strip().split("\n")

                    try:
                        font = PIL.ImageFont.truetype("arial.ttf", size=32)
                    except IOError:
                        font = PIL.ImageFont.load_default()
                    dummy_img = PIL.Image.new("RGB", (10, 10))
                    draw = PIL.ImageDraw.Draw(dummy_img)
                    padding = 10
                    line_height = draw.textbbox((0, 0), "A", font=font)[3] + 2
                    text_width = max(draw.textbbox((0, 0), line, font=font)[2] for line in lines) + 2 * padding
                    text_height = len(lines) * line_height + 2 * padding

                    text_img = PIL.Image.new("RGB", (text_width, text_height), "white")
                    draw = PIL.ImageDraw.Draw(text_img)
                    for i, line in enumerate(lines):
                        draw.text((padding, padding + i * line_height), line, fill="black", font=font)

                    total_width = max(graph_img.width, text_img.width)
                    total_height = graph_img.height + text_img.height
                    combined = PIL.Image.new("RGB", (total_width, total_height), "white")
                    combined.paste(graph_img, (0, 0))
                    combined.paste(text_img, (0, graph_img.height))

                    file_path = os.path.join(folder, f"shape_{window.shape_id}_summary.png")
                    combined.save(file_path)

                except Exception as e:
                    messagebox.showerror("Export Error", f"Failed to export Shape {window.shape_id}:\n{e}")
        messagebox.showinfo("Export Complete", "All shapes have been exported.")

    def compare_selected_shapes(self):
        if not open_windows or len(open_windows) < 2:
            messagebox.showinfo("Compare Shapes", "At least two shapes must be open.")
            return

        selector_window = tk.Toplevel(self.root)
        selector_window.title("Select Shapes to Compare")
        selector_window.geometry("300x400")

        selection_vars = {}

        canvas = tk.Canvas(selector_window)
        scrollbar = tk.Scrollbar(selector_window, orient="vertical", command=canvas.yview)
        frame = tk.Frame(canvas)

        frame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )

        canvas.create_window((0, 0), window=frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        # Scrollable checklist setup
        canvas.pack(side="top", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        # Define the action BEFORE you create the button
        def apply_selection():
            selected_ids = [sid for sid, var in selection_vars.items() if var.get()]
            if len(selected_ids) < 2:
                messagebox.showwarning("Selection Error", "Select at least two shapes to compare.")
                return

            for window in open_windows:
                is_selected = window.shape_id in selected_ids
                window.compare_var.set(is_selected)
                window.update_comparison()

            selector_window.destroy()

        # Bottom frame with the Compare button
        bottom_frame = tk.Frame(selector_window)
        bottom_frame.pack(side="bottom", fill="x", pady=5)

        compare_button = styled_button(bottom_frame, text="Compare", command=apply_selection)
        compare_button.pack(pady=5)

        for window in open_windows:
            var = tk.BooleanVar(value=window.compare_var.get())
            selection_vars[window.shape_id] = var
            tk.Checkbutton(frame, text=f"Shape {window.shape_id}", variable=var).pack(anchor="w")

        bottom_frame = tk.Frame(selector_window)
        bottom_frame.pack(fill="x")

    def export_selected_shapes(self):
        if not open_windows:
            messagebox.showinfo("No Shapes", "There are no shapes to export.")
            return

        selector_window = tk.Toplevel(self.root)
        selector_window.title("Select Shapes to Export")
        selector_window.geometry("300x400")

        selection_vars = {}

        canvas = tk.Canvas(selector_window)
        scrollbar = tk.Scrollbar(selector_window, orient="vertical", command=canvas.yview)
        frame = tk.Frame(canvas)

        frame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )

        canvas.create_window((0, 0), window=frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="top", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        for window in open_windows:
            var = tk.BooleanVar(value=False)
            selection_vars[window.shape_id] = var
            tk.Checkbutton(frame, text=f"Shape {window.shape_id}", variable=var).pack(anchor="w")

        def apply_export():
            selected_windows = [w for w in open_windows if selection_vars[w.shape_id].get()]
            if not selected_windows:
                messagebox.showwarning("No Selection", "Please select at least one shape to export.")
                return

            folder = filedialog.askdirectory(title="Select Folder to Save Selected Summaries")
            if not folder:
                return

            try:
                # Create new folder for this batch
                timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")
                batch_folder = os.path.join(folder, f"export_batch_{timestamp}")
                os.makedirs(batch_folder, exist_ok=True)

                from PIL import Image, ImageDraw, ImageFont
                from io import BytesIO

                for window in selected_windows:
                    if not hasattr(window.analyzer, 'shape'):
                        continue

                    buf = BytesIO()
                    window.analyzer.fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
                    buf.seek(0)
                    graph_img = Image.open(buf).convert("RGB")

                    summary = display_summary(window.analyzer.shape.summary())
                    lines = summary.strip().split("\n")

                    try:
                        font = ImageFont.truetype("arial.ttf", 32)
                    except:
                        font = ImageFont.load_default()

                    dummy_img = Image.new("RGB", (10, 10))
                    draw = ImageDraw.Draw(dummy_img)
                    padding = 10
                    line_height = draw.textbbox((0, 0), "A", font=font)[3] + 2
                    text_width = max(draw.textbbox((0, 0), line, font=font)[2] for line in lines) + 2 * padding
                    text_height = len(lines) * line_height + 2 * padding

                    text_img = Image.new("RGB", (text_width, text_height), "white")
                    draw = ImageDraw.Draw(text_img)
                    for i, line in enumerate(lines):
                        draw.text((padding, padding + i * line_height), line, fill="black", font=font)

                    total_width = max(graph_img.width, text_img.width)
                    total_height = graph_img.height + text_img.height
                    combined = Image.new("RGB", (total_width, total_height), "white")
                    combined.paste(graph_img, (0, 0))
                    combined.paste(text_img, (0, graph_img.height))

                    filename = f"shape_{window.shape_id}_summary.png"
                    file_path = os.path.join(batch_folder, filename)
                    combined.save(file_path)

                messagebox.showinfo("Export Complete", f"{len(selected_windows)} shapes exported to:\n{batch_folder}")
                selector_window.destroy()

            except Exception as e:
                messagebox.showerror("Export Error", f"Something went wrong:\n{e}")

        # Bottom frame with Export Selected button
        bottom_frame = tk.Frame(selector_window)
        bottom_frame.pack(side="bottom", fill="x", pady=5)

        export_btn = styled_button(bottom_frame, text="Export Selected", command=apply_export)
        export_btn.pack(pady=5)



# Run
if __name__ == "__main__":
    root = tk.Tk()
    apply_global_styles(root)

    # Build the UI
    app = GeometryApp(root)

    # Force full screen after a short delay
    def force_maximize():
        root.state('zoomed')

    root.after(100, force_maximize)

    # Set background
    root.configure(bg="#f0f4f8")

    # Start main loop AFTER everything is initialized
    root.mainloop()

